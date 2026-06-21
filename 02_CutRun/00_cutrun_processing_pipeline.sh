#!/usr/bin/env bash
set -euo pipefail

# Process one paired-end CUT&RUN sample through QC, alignment, filtering,
# bigWig generation and MACS2 peak calling.

# Sample, reference and workflow configuration.
SAMPLE="${SAMPLE:?sample identifier}"
FASTQ_R1="${FASTQ_R1:?read 1 FASTQ.gz}"
FASTQ_R2="${FASTQ_R2:?read 2 FASTQ.gz}"
OUTDIR="${OUTDIR:?output directory}"
GENOME_INDEX_PREFIX="${GENOME_INDEX_PREFIX:?Bowtie2 index prefix}"
ADAPTER_REF="${ADAPTER_REF:?BBduk adapter FASTA}"
PICARD_JAR="${PICARD_JAR:?path to picard.jar}"

THREADS="${THREADS:-14}"
PRIMARY_GENOME_PREFIX="${PRIMARY_GENOME_PREFIX:-}"
SPIKEIN_GENOME_PREFIX="${SPIKEIN_GENOME_PREFIX:-}"
MACS_GENOME="${MACS_GENOME:-hs}"
MACS_QVALUE="${MACS_QVALUE:-0.05}"
BAMCOVERAGE_MODE="${BAMCOVERAGE_MODE:-RPKM_BIN1}"
REMOVE_DUPLICATES="${REMOVE_DUPLICATES:-false}"
USE_MODULES="${USE_MODULES:-false}"

if [[ "$USE_MODULES" == "true" ]]; then
  # Load representative tool versions recorded in source workflows.
  module load fastqc/0.11.5
  module load multiqc/1.8
  module load bowtie2/2.3.4.1
  module load samtools/1.9
  module load deeptools/3.5.0
  module load sambamba/0.6.7
  module load picard/3.0.0
  module load macs/2.2.7.1
fi

for command_name in fastqc multiqc bbduk.sh bowtie2 samtools sambamba bamCoverage macs2 java; do
  command -v "$command_name" >/dev/null || {
    echo "Required command not found: $command_name" >&2
    exit 1
  }
done

for input_file in "$FASTQ_R1" "$FASTQ_R2" "$ADAPTER_REF" "$PICARD_JAR"; do
  [[ -f "$input_file" ]] || {
    echo "Input file not found: $input_file" >&2
    exit 1
  }
done

mkdir -p \
  "$OUTDIR/qc/raw" "$OUTDIR/qc/trimmed" "$OUTDIR/qc/bam" \
  "$OUTDIR/trimmed" "$OUTDIR/alignment" "$OUTDIR/bam" \
  "$OUTDIR/bigwig" "$OUTDIR/macs2"

trimmed_r1="$OUTDIR/trimmed/${SAMPLE}_R1.trim.fastq.gz"
trimmed_r2="$OUTDIR/trimmed/${SAMPLE}_R2.trim.fastq.gz"
combined_bam="$OUTDIR/bam/${SAMPLE}.combined.sorted.bam"
primary_bam="$OUTDIR/bam/${SAMPLE}.primary.bam"
filtered_bam="$OUTDIR/bam/${SAMPLE}.filtered.bam"
final_bam="$OUTDIR/bam/${SAMPLE}.final.bam"

# Raw and post-trimming read QC.
fastqc -t "$THREADS" -o "$OUTDIR/qc/raw" "$FASTQ_R1" "$FASTQ_R2"
multiqc --force --filename multiqc_raw_fastqc.html --outdir "$OUTDIR/qc/raw" "$OUTDIR/qc/raw"

bbduk.sh \
  in1="$FASTQ_R1" in2="$FASTQ_R2" \
  out1="$trimmed_r1" out2="$trimmed_r2" \
  ref="$ADAPTER_REF" \
  ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=5 minlen=20 \
  stats="$OUTDIR/trimmed/${SAMPLE}_bbduk.trimstats.txt" \
  refstats="$OUTDIR/trimmed/${SAMPLE}_bbduk.refstats.txt"

fastqc -t "$THREADS" -o "$OUTDIR/qc/trimmed" "$trimmed_r1" "$trimmed_r2"
multiqc --force --filename multiqc_trimmed_fastqc.html --outdir "$OUTDIR/qc/trimmed" "$OUTDIR/qc/trimmed"

# Paired-end Bowtie2 alignment. The flags match the common CUT&RUN workflows.
bowtie2 \
  --dovetail --local --very-sensitive --no-mixed --no-discordant \
  --phred33 -I 10 -X 700 -p "$THREADS" \
  -x "$GENOME_INDEX_PREFIX" -1 "$trimmed_r1" -2 "$trimmed_r2" \
  2> "$OUTDIR/alignment/${SAMPLE}_bowtie2.txt" \
  | samtools view -h -b -@ "$THREADS" - \
  | samtools sort -@ "$THREADS" -o "$combined_bam" -
samtools index -@ "$THREADS" "$combined_bam"

# For hybrid indexes, select primary-genome contigs by prefix. For a normal
# single-genome index, leave PRIMARY_GENOME_PREFIX empty and retain all contigs.
if [[ -n "$PRIMARY_GENOME_PREFIX" ]]; then
  mapfile -t primary_contigs < <(
    samtools idxstats "$combined_bam" | cut -f 1 | grep "^${PRIMARY_GENOME_PREFIX}" || true
  )
  ((${#primary_contigs[@]} > 0)) || {
    echo "No contigs matched PRIMARY_GENOME_PREFIX=$PRIMARY_GENOME_PREFIX" >&2
    exit 1
  }
  samtools view -h -b -@ "$THREADS" "$combined_bam" "${primary_contigs[@]}" \
    | samtools view -h -@ "$THREADS" - \
    | sed "s/${PRIMARY_GENOME_PREFIX}//g" \
    | samtools view -b -@ "$THREADS" -o "$primary_bam" -
else
  cp "$combined_bam" "$primary_bam"
fi
samtools index -@ "$THREADS" "$primary_bam"

# Optionally exclude non-standard chromosomes or other unwanted contigs.
if [[ -n "${EXCLUDE_CHROM_REGEX:-}" ]]; then
  mapfile -t retained_contigs < <(
    samtools idxstats "$primary_bam" | cut -f 1 | grep -Ev "$EXCLUDE_CHROM_REGEX" || true
  )
  ((${#retained_contigs[@]} > 0)) || {
    echo "EXCLUDE_CHROM_REGEX removed every contig" >&2
    exit 1
  }
  samtools view -h -b -@ "$THREADS" "$primary_bam" "${retained_contigs[@]}" > "$filtered_bam"
else
  cp "$primary_bam" "$filtered_bam"
fi

sambamba sort -t "$THREADS" -o "$OUTDIR/bam/${SAMPLE}.filtered.sorted.bam" "$filtered_bam"

# Picard marks duplicates in the source workflow. Set REMOVE_DUPLICATES=true
# only if duplicate removal is intended and reported for the experiment.
java -Xmx8g -jar "$PICARD_JAR" MarkDuplicates \
  INPUT="$OUTDIR/bam/${SAMPLE}.filtered.sorted.bam" \
  OUTPUT="$OUTDIR/bam/${SAMPLE}.dup_processed.bam" \
  METRICS_FILE="$OUTDIR/qc/bam/${SAMPLE}_dup_metrics.txt" \
  REMOVE_DUPLICATES="$REMOVE_DUPLICATES" \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

samtools view -q 2 -F 0x04 -b -@ "$THREADS" \
  "$OUTDIR/bam/${SAMPLE}.dup_processed.bam" \
  | samtools sort -@ "$THREADS" -o "$final_bam" -
samtools index -@ "$THREADS" "$final_bam"

samtools idxstats "$final_bam" > "$OUTDIR/qc/bam/${SAMPLE}_idxstats.txt"
samtools flagstat -@ "$THREADS" "$final_bam" > "$OUTDIR/qc/bam/${SAMPLE}_flagstat.txt"
samtools stats -@ "$THREADS" "$final_bam" > "$OUTDIR/qc/bam/${SAMPLE}_stats.txt"
java -Xmx8g -jar "$PICARD_JAR" CollectInsertSizeMetrics \
  I="$final_bam" \
  O="$OUTDIR/qc/bam/${SAMPLE}_insert_size_metrics.txt" \
  H="$OUTDIR/qc/bam/${SAMPLE}_insert_size_histogram.pdf"
multiqc --force --filename multiqc_alignment.html --outdir "$OUTDIR/qc/bam" "$OUTDIR/qc/bam"

# Count reads assigned to a spike-in/normalization genome in a hybrid index.
if [[ -n "$SPIKEIN_GENOME_PREFIX" ]]; then
  mapfile -t spikein_contigs < <(
    samtools idxstats "$combined_bam" | cut -f 1 | grep "^${SPIKEIN_GENOME_PREFIX}" || true
  )
  if ((${#spikein_contigs[@]} > 0)); then
    samtools view -c -F 260 "$combined_bam" "${spikein_contigs[@]}" \
      > "$OUTDIR/qc/bam/${SAMPLE}_spikein_mapped_reads.txt"
  else
    echo "No contigs matched SPIKEIN_GENOME_PREFIX=$SPIKEIN_GENOME_PREFIX" >&2
    exit 1
  fi
fi

case "$BAMCOVERAGE_MODE" in
  RPKM_BIN1)
    bamCoverage -b "$final_bam" -o "$OUTDIR/bigwig/${SAMPLE}.bw" \
      --binSize 1 --normalizeUsing RPKM --numberOfProcessors 5
    ;;
  BPM_BIN20_SMOOTH60_EXTEND150)
    bamCoverage -b "$final_bam" -o "$OUTDIR/bigwig/${SAMPLE}.bw" \
      --binSize 20 --normalizeUsing BPM --smoothLength 60 \
      --extendReads 150 --numberOfProcessors 5
    ;;
  SCALE_FACTOR_BIN20_SMOOTH60)
    bamCoverage -b "$final_bam" -o "$OUTDIR/bigwig/${SAMPLE}.bw" \
      --binSize 20 --smoothLength 60 --extendReads \
      --scaleFactor "${SCALE_FACTOR:?required for scale-factor normalization}" \
      --numberOfProcessors "$THREADS"
    ;;
  *)
    echo "Unknown BAMCOVERAGE_MODE: $BAMCOVERAGE_MODE" >&2
    exit 1
    ;;
esac

macs_args=(
  callpeak -t "$final_bam" -f BAMPE -g "$MACS_GENOME"
  -q "$MACS_QVALUE" -n "$SAMPLE" --outdir "$OUTDIR/macs2"
)
if [[ -n "${CONTROL_BAM:-}" ]]; then
  [[ -f "$CONTROL_BAM" ]] || {
    echo "CONTROL_BAM not found: $CONTROL_BAM" >&2
    exit 1
  }
  macs_args+=(-c "$CONTROL_BAM")
fi
macs2 "${macs_args[@]}"
