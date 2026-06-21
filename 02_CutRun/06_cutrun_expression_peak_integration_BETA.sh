#!/usr/bin/env bash
set -euo pipefail

# Integrate one CUT&RUN peak set with differential expression using BETA.
# Inputs are a peak BED, BETA-format expression table and RefSeq reference.

PEAK_FILE="${PEAK_FILE:?peak BED file}"
EXPR_FILE="${EXPR_FILE:?expression table in BETA-compatible format}"
REFERENCE_REFSEQ="${REFERENCE_REFSEQ:?BETA RefSeq reference file}"
OUTPUT_DIR="${OUTPUT_DIR:?output directory}"

BETA_BIN="${BETA_BIN:-BETA}"
BETA_MODE="${BETA_MODE:-basic}"
BETA_INFO="${BETA_INFO:-11,3,7}"
GENOME="${GENOME:-hg38}"
PEAK_NAME="${PEAK_NAME:-$(basename "$PEAK_FILE" .bed)}"
DIFF_FDR="${DIFF_FDR:-0.05}"
PEAK_NUMBER="${PEAK_NUMBER:-100000}"

for input_file in "$PEAK_FILE" "$EXPR_FILE" "$REFERENCE_REFSEQ"; do
  [[ -f "$input_file" ]] || {
    echo "Input file not found: $input_file" >&2
    exit 1
  }
done
command -v "$BETA_BIN" >/dev/null || {
  echo "BETA executable not found: $BETA_BIN" >&2
  exit 1
}

# Use source promoter/non-promoter distance rules unless overridden.
if [[ -n "${DISTANCE:-}" ]]; then
  beta_distance="$DISTANCE"
elif [[ "$PEAK_NAME" == *"_nonpromoter_"* || "$PEAK_NAME" == *"_nonpromoter" ]]; then
  beta_distance=100000
elif [[ "$PEAK_NAME" == *"_promoter_"* || "$PEAK_NAME" == *"_promoter" ]]; then
  beta_distance=2000
else
  echo "Set DISTANCE explicitly when PEAK_NAME is not labelled promoter/nonpromoter." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR/$PEAK_NAME"

"$BETA_BIN" "$BETA_MODE" \
  -p "$PEAK_FILE" \
  -e "$EXPR_FILE" \
  -k O \
  --info "$BETA_INFO" \
  -g "$GENOME" \
  -r "$REFERENCE_REFSEQ" \
  -n "$PEAK_NAME" \
  --output "$OUTPUT_DIR/$PEAK_NAME" \
  --gname2 \
  --da 1 \
  --df "$DIFF_FDR" \
  -c 1 \
  -d "$beta_distance" \
  --pn "$PEAK_NUMBER"

echo "BETA analysis complete: $OUTPUT_DIR/$PEAK_NAME"
