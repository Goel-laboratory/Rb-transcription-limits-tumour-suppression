#!/usr/bin/env bash
set -euo pipefail

# Generate deepTools matrix, heatmap and profile panels from one BED and
# matched bigWig/label lists.

BED_FILE="${BED_FILE:?BED file of regions to plot}"
BIGWIGS_CSV="${BIGWIGS_CSV:?comma-separated bigWig paths}"
SAMPLE_LABELS_CSV="${SAMPLE_LABELS_CSV:?comma-separated labels matching bigWigs}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-deeptools_panel}"
OUTPUT_DIR="${OUTPUT_DIR:-code_availability/Figures/output/08_deeptools}"
THREADS="${THREADS:-8}"
UPSTREAM="${UPSTREAM:-3000}"
DOWNSTREAM="${DOWNSTREAM:-3000}"
REGION_BODY_LENGTH="${REGION_BODY_LENGTH:-3000}"
COLOR_MAP="${COLOR_MAP:-Blues}"

IFS=',' read -r -a bigwigs <<< "$BIGWIGS_CSV"
IFS=',' read -r -a sample_labels <<< "$SAMPLE_LABELS_CSV"
if ((${#bigwigs[@]} != ${#sample_labels[@]})); then
  echo "BIGWIGS_CSV and SAMPLE_LABELS_CSV must have equal lengths." >&2
  exit 1
fi
[[ -f "$BED_FILE" ]] || { echo "BED_FILE not found: $BED_FILE" >&2; exit 1; }
for bigwig in "${bigwigs[@]}"; do
  [[ -f "$bigwig" ]] || { echo "bigWig not found: $bigwig" >&2; exit 1; }
done

mkdir -p "$OUTPUT_DIR"
matrix="$OUTPUT_DIR/${OUTPUT_PREFIX}_matrix.gz"

# Build one signal matrix used by both downstream plots.
computeMatrix scale-regions \
  -S "${bigwigs[@]}" \
  -R "$BED_FILE" \
  --beforeRegionStartLength "$UPSTREAM" \
  --regionBodyLength "$REGION_BODY_LENGTH" \
  --afterRegionStartLength "$DOWNSTREAM" \
  --skipZeros --missingDataAsZero \
  --numberOfProcessors "$THREADS" \
  -o "$matrix"

# Render the region-by-sample heatmap from the saved matrix.
plotHeatmap -m "$matrix" \
  --samplesLabel "${sample_labels[@]}" \
  --heatmapHeight 12 --heatmapWidth 4 \
  --colorMap "$COLOR_MAP" \
  -out "$OUTPUT_DIR/${OUTPUT_PREFIX}_heatmap.pdf"

# Summarize the same matrix as average signal profiles.
plotProfile -m "$matrix" \
  --samplesLabel "${sample_labels[@]}" \
  --perGroup \
  -out "$OUTPUT_DIR/${OUTPUT_PREFIX}_profile.pdf"
