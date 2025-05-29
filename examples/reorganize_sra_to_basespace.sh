#!/usr/bin/env bash
# reorganize_sra_to_basespace.sh
# Usage: ./reorganize_sra_to_basespace.sh <id_list.txt> <src_dir> <dest_root>

set -euo pipefail

ID_FILE="$1"
SRC_DIR="$2"
DEST_ROOT="$3"

# Make sure destination root exists
mkdir -p "$DEST_ROOT"

counter=1
while IFS= read -r ID; do
  # Create a sample-specific folder (<ID>-ds) under the BaseSpace run directory
  DS_DIR="${DEST_ROOT}/${ID}-ds"
  mkdir -p "$DS_DIR"

  # For each read pair (R1/R2), look for .fastq.gz or .fastq, gzip if needed, and rename to BaseSpace convention
  for RP in 1 2; do
    for EXT in fastq.gz fastq; do
      SRC_FASTQ="${SRC_DIR}/${ID}_${RP}.${EXT}"
      if [ -f "$SRC_FASTQ" ]; then
        if [ "$EXT" = "fastq" ]; then
          gzip -c "$SRC_FASTQ" > "${DS_DIR}/${ID}_S${counter}_L001_R${RP}_001.fastq.gz"
        else
          cp "$SRC_FASTQ" "${DS_DIR}/${ID}_S${counter}_L001_R${RP}_001.fastq.gz"
        fi
        break
      fi
    done
  done

  counter=$((counter + 1))
done < "$ID_FILE"

