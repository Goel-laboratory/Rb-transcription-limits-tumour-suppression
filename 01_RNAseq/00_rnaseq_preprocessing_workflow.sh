#!/usr/bin/env bash

# Process RNA-seq FASTQs through QC, trimming, STAR alignment, featureCounts
# and coverage-track generation. Inputs are configured through environment variables.

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:?project root}"
FASTQ_DIR="${FASTQ_DIR:?directory containing FASTQ files}"
OUTDIR="${RNASEQ_OUTPUT_DIR:-${PROJECT_DIR}/code_availability/RNAseq/output}"
TRIM_DIR="${PROJECT_DIR}/trim_data"
BAM_DIR="${OUTDIR}/bamfiles"
BIGWIG_DIR="${OUTDIR}/bigwig"
MULTIQC_DIR="${OUTDIR}/multiqc"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

GENOME_INDEX="${GENOME_INDEX:?STAR genome index}"
GTF="${GTF:?gene annotation GTF}"
BBDUK_DIR="${BBDUK_DIR:?BBMap installation directory}"
BBDUK_ADAPTERS="${BBDUK_DIR}/resources/truseq.fa.gz"

THREADS="${THREADS:-8}"
FASTQ_GLOB="${FASTQ_GLOB:-*.fastq.gz}"
READ_LAYOUT="${READ_LAYOUT:-single}"
USE_MODULES="${USE_MODULES:-false}"

case "$READ_LAYOUT" in
  single) SINGLE_END=true; PAIRED_END=false ;;
  paired) SINGLE_END=false; PAIRED_END=true ;;
  *) echo "READ_LAYOUT must be single or paired" >&2; exit 1 ;;
esac

for path in "$FASTQ_DIR" "$GENOME_INDEX" "$GTF" "$BBDUK_ADAPTERS"; do
  [[ -e "$path" ]] || { echo "Required input not found: $path" >&2; exit 1; }
done

# Create shared output directories before running the individual workflow stages.
mkdir -p "${OUTDIR}" "${TRIM_DIR}" "${BAM_DIR}" "${BIGWIG_DIR}" "${MULTIQC_DIR}" "${SCRIPTS_DIR}"

###############################################################################
# 1. Raw-read QC
###############################################################################

if [[ "$USE_MODULES" == true ]]; then
  module load fastqc/0.11.5
  module load multiqc/1.8
  module load java/1.8.0_312-jdk
  module load star/2.7.5b
  module load subread/2.0.6
  module load samtools/1.9
  module load deeptools/3.5.0
fi

mkdir -p "${PROJECT_DIR}/raw_data/fastqc"
cd "${FASTQ_DIR}"
fastqc -t "${THREADS}" -o "${PROJECT_DIR}/raw_data/fastqc" ${FASTQ_GLOB}
multiqc --filename multiqc_RNAseq_raw "${PROJECT_DIR}/raw_data/fastqc"

###############################################################################
# 2. Adapter and quality trimming
###############################################################################


mkdir -p "${TRIM_DIR}/fastqc_out_trim"
cd "${FASTQ_DIR}"

if [[ "${SINGLE_END}" == true ]]; then
  for read1 in ${FASTQ_GLOB}; do
    [[ -f "${read1}" ]] || continue
    sample=$(basename "${read1}" .fastq.gz)
    "${BBDUK_DIR}/bbduk.sh" \
      in="${read1}" \
      out="${TRIM_DIR}/${sample}.trim.fastq.gz" \
      ref="${BBDUK_ADAPTERS}" \
      ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=5 minlen=20 \
      stats="${TRIM_DIR}/${sample}_bbduk.trimstats.txt" \
      refstats="${TRIM_DIR}/${sample}_bbduk.refstats.txt" \
      &> "${TRIM_DIR}/${sample}_bbduk.stdout_stats.txt"
  done
fi

if [[ "${PAIRED_END}" == true ]]; then
  for read1 in *_R1*.fastq.gz; do
    [[ -f "${read1}" ]] || continue
    read2="${read1/_R1/_R2}"
    sample=$(basename "${read1}" _R1.fastq.gz)
    "${BBDUK_DIR}/bbduk.sh" \
      in1="${read1}" in2="${read2}" \
      out1="${TRIM_DIR}/${sample}_R1.trim.fastq.gz" \
      out2="${TRIM_DIR}/${sample}_R2.trim.fastq.gz" \
      ref="${BBDUK_ADAPTERS}" \
      ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=5 minlen=20 \
      stats="${TRIM_DIR}/${sample}_bbduk.trimstats.txt" \
      refstats="${TRIM_DIR}/${sample}_bbduk.refstats.txt" \
      &> "${TRIM_DIR}/${sample}_bbduk.stdout_stats.txt"
  done
fi

cd "${TRIM_DIR}"
fastqc -t "${THREADS}" -o fastqc_out_trim *.trim.fastq.gz
multiqc --filename fastqc_out_trim fastqc_out_trim

###############################################################################
# 3. STAR alignment
###############################################################################


cd "${TRIM_DIR}"
if [[ "${SINGLE_END}" == true ]]; then
  for read1 in *.trim.fastq.gz; do
    [[ -f "${read1}" ]] || continue
    sample=$(basename "${read1}" .trim.fastq.gz)
    STAR --runThreadN "${THREADS}" \
      --genomeDir "${GENOME_INDEX}" \
      --readFilesIn "${read1}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${BAM_DIR}/${sample}." \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMunmapped Within \
      --outSAMattributes All
  done
fi

if [[ "${PAIRED_END}" == true ]]; then
  for read1 in *_R1.trim.fastq.gz; do
    [[ -f "${read1}" ]] || continue
    read2="${read1/_R1.trim.fastq.gz/_R2.trim.fastq.gz}"
    sample=$(basename "${read1}" _R1.trim.fastq.gz)
    STAR --runThreadN "${THREADS}" \
      --genomeDir "${GENOME_INDEX}" \
      --readFilesIn "${read1}" "${read2}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${BAM_DIR}/${sample}." \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMunmapped Within \
      --outSAMattributes All
  done
fi

###############################################################################
# 4. Gene-level counts with featureCounts
###############################################################################
# Original workflows counted exons by gene_id for all strandedness modes.
# Use the summary files to choose the strandedness with the highest Assigned reads.
# most downstream analyses used featurecounts_stranded_1; one current workflow
# explicitly concluded reverse-stranded / strand 2 for its library. Do not assume:
# inspect the summary files for each dataset.

mkdir -p "${BAM_DIR}/featurecounts_unstranded" \
         "${BAM_DIR}/featurecounts_stranded_1" \
         "${BAM_DIR}/featurecounts_stranded_2"

featureCounts -T "${THREADS}" -s 0 -t exon -g gene_id \
  -a "${GTF}" \
  -o "${BAM_DIR}/featurecounts_unstranded/featurecounts_unstranded.txt" \
  "${BAM_DIR}"/*Aligned.sortedByCoord.out.bam

featureCounts -T "${THREADS}" -s 1 -t exon -g gene_id \
  -a "${GTF}" \
  -o "${BAM_DIR}/featurecounts_stranded_1/featurecounts_stranded_1.txt" \
  "${BAM_DIR}"/*Aligned.sortedByCoord.out.bam

featureCounts -T "${THREADS}" -s 2 -t exon -g gene_id \
  -a "${GTF}" \
  -o "${BAM_DIR}/featurecounts_stranded_2/featurecounts_stranded_2.txt" \
  "${BAM_DIR}"/*Aligned.sortedByCoord.out.bam

###############################################################################
# 5. BAM index and bigWig generation
###############################################################################
# Optional. bigWigs used deepTools bamCoverage
# with BPM normalization.

mkdir -p "${BIGWIG_DIR}"
for bam in "${BAM_DIR}"/*Aligned.sortedByCoord.out.bam; do
  [[ -f "${bam}" ]] || continue
  sample=$(basename "${bam}" Aligned.sortedByCoord.out.bam)
  samtools index -@ "${THREADS}" "${bam}"
  bamCoverage -b "${bam}" \
    -o "${BIGWIG_DIR}/${sample}.bw" \
    --normalizeUsing BPM \
    -p "${THREADS}"
done

###############################################################################
# 6. Final MultiQC report
###############################################################################

mkdir -p "${MULTIQC_DIR}"
multiqc --filename multiqc_summary \
  --outdir "${MULTIQC_DIR}" \
  "${TRIM_DIR}/fastqc_out_trim" \
  "${TRIM_DIR}"/*_stats.txt \
  "${BAM_DIR}"/*Log.final.out \
  "${BAM_DIR}"/featurecounts_*/*.summary

###############################################################################
# 7. Strand-choice helper
###############################################################################
# Inspect Assigned rows in:
#   featurecounts_unstranded/featurecounts_unstranded.txt.summary
#   featurecounts_stranded_1/featurecounts_stranded_1.txt.summary
#   featurecounts_stranded_2/featurecounts_stranded_2.txt.summary
# Choose the strandedness with the highest Assigned counts for downstream DESeq2.
