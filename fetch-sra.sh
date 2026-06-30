#!/usr/bin/env bash
#
# fetch-sra.sh — download paired-end FASTQs from the NCBI SRA for this pipeline.
#
# Usage:
#   ./fetch-sra.sh <accession_list.txt> [output_directory]
#
#   <accession_list.txt>  one SRA run accession per line (e.g. SRR26513986),
#                         such as examples/canine-h3n2-sra/ids.txt
#   [output_directory]    where to save the FASTQs (default: ./sra-data)
#
# Requires the NCBI SRA Toolkit (fasterq-dump). It is NOT installed by conda for
# this pipeline — install the official build and configure it once:
#   https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
#   https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration

set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <accession_list.txt> [output_directory]" >&2
  exit 1
fi

accession_list="$1"
output_dir="${2:-sra-data}"

if ! command -v fasterq-dump >/dev/null 2>&1; then
  echo "Error: fasterq-dump not found. Install the NCBI SRA Toolkit:" >&2
  echo "  https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit" >&2
  exit 1
fi

mkdir -p "$output_dir"

# Download each accession listed in the text file, one line at a time.
while IFS= read -r accession; do
  [ -z "$accession" ] && continue   # skip blank lines

  if [ -f "$output_dir/${accession}_1.fastq.gz" ]; then
    echo "Already have $accession, skipping."
    continue
  fi

  echo "Downloading $accession ..."
  fasterq-dump --split-files --outdir "$output_dir" "$accession"
  gzip "$output_dir/${accession}_1.fastq" "$output_dir/${accession}_2.fastq"
done < "$accession_list"

# Write a starter metadata file with one row per accession. You fill in the
# SampleId and Replicate columns by hand. We never overwrite an existing one.
metadata="data/metadata.tsv"
if [ -f "$metadata" ]; then
  echo "Note: $metadata already exists, leaving it untouched."
else
  mkdir -p data
  printf 'SequencingId\tSampleId\tReplicate\n' > "$metadata"
  while IFS= read -r accession; do
    [ -z "$accession" ] && continue
    printf '%s\t\t\n' "$accession" >> "$metadata"
  done < "$accession_list"
  echo "Wrote a starter metadata file: $metadata"
fi

abs_output_dir="$(cd "$output_dir" && pwd)"
echo
echo "Done. FASTQs are in: $abs_output_dir"
echo
echo "Next steps:"
echo "  1. Set this line in config.yml:"
echo "         data_root_directory: \"$abs_output_dir\""
echo "  2. Fill in the SampleId and Replicate columns in $metadata."
echo "     For SRA, each accession is usually its own sample, single replicate."
echo "  3. Run:  python mlip/dataflow.py flow --sra-mode"
