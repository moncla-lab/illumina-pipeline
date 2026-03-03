# Build-DB: Annotated Influenza Reference Database Pipeline

Self-contained Snakemake pipeline that builds per-segment VAPOR reference databases
from NCBI Virus sequences. The output is consumed directly by the parent illumina
pipeline's VAPOR-based reference selection step.

## Quick Start

```bash
cd build-db
snakemake --cores 8
```

## Prerequisites

- **Conda environment**: `mlip` (same as the parent pipeline)
- **CLI tools**: `efetch` (EDirect), `mafft`, `tn93`, `tn93-cluster`
- **Python packages**: BioPython (included in `mlip` env)
- **Input data**: `sequences.fasta` — downloaded from NCBI Virus
  (H5NX, currently 217K sequences from Feb 2026)

## Input

| File | Description |
|------|-------------|
| `sequences.fasta` | Multi-FASTA from NCBI Virus. Each record's header contains segment keywords and strain name. |
| `config.yml` | Pipeline configuration (NC_ accessions, thresholds, clustering parameters). |

## Output

For each of the 8 influenza segments (`pb2`, `pb1`, `pa`, `ha`, `np`, `na`, `mp`, `ns`):

```
output/db/{segment}/
  all.fasta           # Multi-FASTA for VAPOR reference selection
  metadata.gb         # NC_ reference GenBank (default fallback annotation)
  sequence.fasta      # NC_ reference FASTA
  genbanks/
    {accession}.gb    # Per-sequence annotated GenBank files
```

This matches the contract expected by the parent pipeline when `use_vapor: true`.

## Pipeline Steps

### 1. Fetch NC_ References (`fetch_nc_references`)

Fetches 8 GenBank records from NCBI for the A/goose/Guangdong/1/1996 reference
genome (NC_007357–NC_007364). These have complete, curated CDS annotations
including M1/M2 splicing, PA-X frameshift, and PB1-F2. Requests are serialized
with a 1-second delay to avoid NCBI rate limiting. FASTA sequences are extracted
from the GenBanks locally.

### 2. Parse & Filter Sequences (`parse_sequences`)

Parses all 217K input sequences in a single pass:

- **Segment identification**: Keyword voting system matches description text
  against segment-specific keywords (e.g., "hemagglutinin" → HA, "matrix protein"
  → MP). Handles ambiguous descriptions by picking the segment with the most
  keyword matches.
- **Strain extraction**: Parses `A/{host}/{location}/{id}/{year}` from headers.
- **Uni12/Uni13 detection**: Checks for conserved influenza terminal sequences
  via Hamming distance (threshold configurable, default 1 mismatch). This is
  recorded as metadata but **not used as a filter** — truncated sequences are
  kept because alignment against the NC_ reference handles missing termini.
- **Complete genomes**: Identifies strains with all 8 segments present.
- **Deduplication**: When a strain has multiple sequences for the same segment,
  keeps only the longest.

QC files are written to `intermediate/qc/` (ambiguous segments, unidentified
segments, missing strain names).

### 3. Align to References (`align_to_reference`)

Aligns all sequences for each segment to the corresponding NC_ reference using
MAFFT `--addfragments`. This produces an MSA with the reference as the first
sequence and ~23K queries aligned to it. The alignment provides the coordinate
mapping needed for annotation transfer.

Runs 8 times in parallel (one per segment), each using 4 threads.

### 4. Create GenBank Files (`create_genbank_files`)

Transfers CDS annotations from the NC_ GenBank to each query sequence using the
MAFFT alignment from step 3:

1. Build a coordinate map from reference ungapped positions to alignment columns
2. Build a reverse map from alignment columns to query ungapped positions
3. For each CDS feature (including CompoundLocation for spliced genes like M1/M2),
   map the coordinates through: `ref_pos → alignment_col → query_pos`
4. Truncated sequences get partial feature annotations (`BeforePosition`/
   `AfterPosition`)

Output: one `.gb` file per sequence in `intermediate/genbanks/{segment}/`.

### 5. Parse Strain Metadata (`parse_strain_metadata`)

Extracts host, country, and year from strain names. Countries are normalized
(e.g., US state names → "USA", Chinese provinces → "China"). Years are binned
into configurable windows (default 2 years: 2024-2025, 2022-2023, etc.).

If GenoFLU lineage results are available, per-segment lineage assignments are
incorporated into the grouping.

Groups are defined as `(country, year_bin)` or `(country, year_bin, lineage)`.

### 6. TN93 Clustering (`tn93_cluster_by_group`)

Performs diversity-aware clustering within each metadata group using TN93
(Tamura-Nei 93) pairwise distances:

- Groups with ≤ `min_seqs_per_group` sequences are kept as-is
- Larger groups are clustered at the configured threshold (default 1.5%)
- Cluster centroids (representative sequences) are kept
- All group centroids are concatenated into the clustered FASTA

This ensures the VAPOR database has good spatial, temporal, and genetic diversity
without being oversized. Heavily-sampled lineages/regions don't dominate.

### 7. Assemble Database (`assemble_database`)

Combines the clustered sequences with the NC_ reference into the final output
structure. Copies the relevant GenBank files from the annotation step. The NC_
reference is always included in `all.fasta`.

## DAG

```
sequences.fasta           NC_ accessions (config)
       |                        |
       v                        v
  parse_sequences         fetch_nc_references
       |                    |           |
       +----+----+          |           |
       |    |    |          |           |
       v    |    v          v           |
  parse_    |  align_to_reference       |
  strain_   |        |                  |
  metadata  |        v                  |
       |    |  create_genbank_files     |
       v    |        |                  |
  tn93_     |        |                  |
  cluster_  |        |                  |
  by_group  |        |                  |
       |    |        |                  |
       +----+--------+------------------+
            |
            v
     assemble_database
            |
            v
    output/db/{segment}/
```

Two independent branches run in parallel:
- **Clustering branch**: `parse_sequences` → `parse_strain_metadata` → `tn93_cluster_by_group`
- **Annotation branch**: `parse_sequences` + `fetch_nc_references` → `align_to_reference` → `create_genbank_files`

Both converge at `assemble_database`.

## Configuration

`config.yml` parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nc_references` | NC_007357–64 | GenBank accessions for A/goose/Guangdong/1/1996 |
| `segments` | all 8 | Segment names (lowercase) |
| `uni12_variants` | 2 variants | Conserved 5' terminal sequences |
| `uni13_rc` | `ccttgtttctact` | Conserved 3' terminal reverse complement |
| `hamming_distance_threshold` | 1 | Max mismatches for Uni12/13 detection |
| `input_fasta` | `sequences.fasta` | Path to NCBI Virus download |
| `clustering.tn93_threshold` | 0.015 | TN93 distance threshold (1.5%) |
| `clustering.year_bin_size` | 2 | Years per time bin |
| `clustering.min_seqs_per_group` | 1 | Minimum sequences kept per metadata group |
| `clustering.genoflu_results` | null | Optional path to GenoFLU results.tsv |

## Integration with Parent Pipeline

Copy the output into an analysis directory's reference folder:

```bash
cp -r output/db/* /path/to/analysis/reference/
```

Then set `use_vapor: true` in the parent pipeline's `config.yml`. The parent
pipeline expects:
- `reference/{segment}/all.fasta` — VAPOR database
- `reference/{segment}/metadata.gb` — default GenBank annotation
- `reference/{segment}/sequence.fasta` — default reference FASTA
- `reference/{segment}/genbanks/{accession}.gb` — per-accession GenBanks
