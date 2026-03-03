"""
build_db.py - Helper functions for the build-db Snakemake pipeline.

Parses NCBI Virus FASTA downloads, identifies segments, filters for complete
genomes, transfers CDS annotations via alignment, and performs metadata-aware
TN93 clustering to build per-segment VAPOR reference databases.
"""

import csv
import json
import os
import re
import subprocess
import tempfile
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    FeatureLocation,
    SeqFeature,
)
from Bio.SeqRecord import SeqRecord


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SEGMENT_KEYWORDS = {
    "pb2": ["polymerase PB2", "segment 1", "PB2 gene", "PB2)"],
    "pb1": ["polymerase PB1", "segment 2", "PB1 gene", "PB1)", "PB1-F2"],
    "pa": ["polymerase PA", "segment 3", "PA gene", "PA)", "PA-X"],
    "ha": ["hemagglutinin", "segment 4", "HA gene", "HA)"],
    "np": [
        "nucleocapsid protein",
        "nucleoprotein",
        "segment 5",
        "NP gene",
        "NP)",
    ],
    "na": ["neuraminidase", "segment 6", "NA gene", "NA)"],
    "mp": [
        "matrix protein",
        "segment 7",
        "M gene",
        "M1 gene",
        "M2 gene",
        "matrix protein 1",
        "matrix protein 2",
    ],
    "ns": [
        "nonstructural protein",
        "segment 8",
        "NS1",
        "NS2",
        "NEP",
        "NS gene",
        "nuclear export protein",
    ],
}

REQUIRED_SEGMENTS = {"pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"}


# ---------------------------------------------------------------------------
# Segment identification
# ---------------------------------------------------------------------------


def identify_segment(description):
    """Identify influenza segment from a FASTA description via keyword voting.

    Returns (segment, votes, is_ambiguous, is_critical).
    - segment: winning segment name (lowercase) or None if tied/no match
    - votes: dict of {segment: vote_count}
    - is_ambiguous: matched keywords from more than one segment
    - is_critical: tie at max votes (cannot resolve)
    """
    desc_lower = description.lower()
    votes = {}
    for segment, keywords in SEGMENT_KEYWORDS.items():
        segment_votes = sum(
            1 for kw in keywords if kw.lower() in desc_lower
        )
        if segment_votes > 0:
            votes[segment] = segment_votes

    if not votes:
        return None, votes, False, False

    max_votes = max(votes.values())
    winners = [seg for seg, v in votes.items() if v == max_votes]
    is_ambiguous = len(votes) > 1
    is_critical = len(winners) > 1
    segment = winners[0] if len(winners) == 1 else None
    return segment, votes, is_ambiguous, is_critical


# ---------------------------------------------------------------------------
# Strain name extraction
# ---------------------------------------------------------------------------


def extract_strain_name(description):
    """Extract A/host/location/.../year from an NCBI FASTA header."""
    match = re.search(r"A/.+?/\d{4}", description)
    if not match:
        match = re.search(r"A/.+?/\d{2}(?!\d)", description)
    if not match:
        return None
    return match.group(0)


# ---------------------------------------------------------------------------
# Uni12 / Uni13 primer checking
# ---------------------------------------------------------------------------


def within_hamming_distance(s1, s2, max_distance=1):
    """Check whether two equal-length strings are within Hamming distance."""
    if len(s1) != len(s2):
        return False
    return sum(c1 != c2 for c1, c2 in zip(s1, s2)) <= max_distance


def check_uni12_uni13(seq_str, uni12_variants, uni13_rc, hamming_threshold=1):
    """Check whether a sequence has Uni12 at 5' and Uni13-RC at 3'.

    Returns (has_uni12, has_uni13).
    """
    seq_lc = seq_str.lower()
    has_uni12 = any(
        within_hamming_distance(seq_lc[: len(v)], v, hamming_threshold)
        for v in uni12_variants
    )
    has_uni13 = within_hamming_distance(
        seq_lc[-len(uni13_rc) :], uni13_rc, hamming_threshold
    )
    return has_uni12, has_uni13


# ---------------------------------------------------------------------------
# Main parsing pass
# ---------------------------------------------------------------------------


def parse_all_sequences(
    input_fasta,
    uni12_variants,
    uni13_rc,
    hamming_threshold=1,
    qc_dir=None,
):
    """Parse an NCBI Virus FASTA, identify segments, strains, and primers.

    Returns:
        strain_segment_seqs: {strain: {segment: [(record, seq_len, has_both)]}}
        strain_segments:     {strain: set_of_segments}
        stats:               dict of parsing statistics
    """
    strain_segment_seqs = defaultdict(lambda: defaultdict(list))
    strain_segments = defaultdict(set)

    total = 0
    no_segment = 0
    critical_ambiguous = 0
    no_strain = 0

    ambiguous_records = []
    no_segment_records = []
    no_strain_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        total += 1
        desc = record.description
        seq_str = str(record.seq)

        segment, votes, is_ambiguous, is_critical = identify_segment(desc)
        if segment is None:
            if is_critical:
                critical_ambiguous += 1
                ambiguous_records.append(
                    (record.id, desc, dict(votes))
                )
            else:
                no_segment += 1
                no_segment_records.append((record.id, desc))
            continue

        strain = extract_strain_name(desc)
        if strain is None:
            no_strain += 1
            no_strain_records.append((record.id, desc))
            continue

        has_uni12, has_uni13 = check_uni12_uni13(
            seq_str, uni12_variants, uni13_rc, hamming_threshold
        )
        has_both = has_uni12 and has_uni13

        strain_segments[strain].add(segment)
        strain_segment_seqs[strain][segment].append(
            (record, len(seq_str), has_both)
        )

    # Write QC files
    if qc_dir:
        os.makedirs(qc_dir, exist_ok=True)
        _write_qc_files(
            qc_dir,
            ambiguous_records,
            no_segment_records,
            no_strain_records,
        )

    stats = {
        "total_sequences": total,
        "no_segment_identified": no_segment,
        "critical_ambiguous": critical_ambiguous,
        "no_strain_identified": no_strain,
        "total_strains": len(strain_segments),
    }
    return dict(strain_segment_seqs), dict(strain_segments), stats


def _write_qc_files(qc_dir, ambiguous, no_segment, no_strain):
    with open(os.path.join(qc_dir, "ambiguous_segments.tsv"), "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "description", "votes"])
        for acc, desc, votes in ambiguous:
            writer.writerow([acc, desc, json.dumps(votes)])

    with open(os.path.join(qc_dir, "no_segment.tsv"), "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "description"])
        for acc, desc in no_segment:
            writer.writerow([acc, desc])

    with open(os.path.join(qc_dir, "no_strain.tsv"), "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "description"])
        for acc, desc in no_strain:
            writer.writerow([acc, desc])


# ---------------------------------------------------------------------------
# Complete genomes & deduplication
# ---------------------------------------------------------------------------


def find_complete_genomes(strain_segments):
    """Return set of strain names that have all 8 segments."""
    return {
        strain
        for strain, segs in strain_segments.items()
        if segs >= REQUIRED_SEGMENTS
    }


def deduplicate_sequences(strain_segment_seqs, complete_strains):
    """Keep the longest sequence per segment per strain (complete genomes only).

    Returns:
        segment_records: {segment: [SeqRecord, ...]}
        strain_segment_map: {strain: {segment: accession}}
    """
    segment_records = defaultdict(list)
    strain_segment_map = {}

    for strain in sorted(complete_strains):
        strain_segment_map[strain] = {}
        for segment in REQUIRED_SEGMENTS:
            seqs = strain_segment_seqs.get(strain, {}).get(segment, [])
            if not seqs:
                continue
            best_record, best_len, best_primers = max(
                seqs, key=lambda x: x[1]
            )
            segment_records[segment].append(best_record)
            strain_segment_map[strain][segment] = best_record.id

    return dict(segment_records), strain_segment_map


def write_segment_fastas(segment_records, output_dir):
    """Write per-segment FASTA files from deduplicated records."""
    os.makedirs(output_dir, exist_ok=True)
    for segment, records in segment_records.items():
        out_path = os.path.join(output_dir, f"{segment}.fasta")
        SeqIO.write(records, out_path, "fasta")


# ---------------------------------------------------------------------------
# Coordinate mapping (ported from mlip/dataflow.py)
# ---------------------------------------------------------------------------


def build_coordinate_map(aligned_seq_str):
    """Map 1-based ungapped position to 1-based alignment column."""
    coord_map = {}
    ungapped_pos = 0
    for col_idx, base in enumerate(aligned_seq_str):
        if base != "-":
            ungapped_pos += 1
            coord_map[ungapped_pos] = col_idx + 1
    return coord_map


def build_reverse_coordinate_map(aligned_seq_str):
    """Map 1-based alignment column to 1-based ungapped position."""
    rev_map = {}
    ungapped_pos = 0
    for col_idx, base in enumerate(aligned_seq_str):
        if base != "-":
            ungapped_pos += 1
            rev_map[col_idx + 1] = ungapped_pos
    return rev_map


# ---------------------------------------------------------------------------
# Annotation transfer
# ---------------------------------------------------------------------------


def transfer_cds_features(ref_genbank, ref_aligned_str, query_aligned_str, query_seq_str):
    """Transfer CDS features from an NC_ GenBank to a query sequence via alignment.

    Args:
        ref_genbank: path to NC_ GenBank file
        ref_aligned_str: aligned reference sequence string (with gaps)
        query_aligned_str: aligned query sequence string (with gaps)
        query_seq_str: ungapped query sequence string

    Returns:
        list of SeqFeature objects with transferred coordinates
    """
    ref_record = SeqIO.read(ref_genbank, "genbank")
    ref_coord_map = build_coordinate_map(ref_aligned_str)
    query_rev_map = build_reverse_coordinate_map(query_aligned_str)

    query_ungapped_len = len(query_seq_str)
    features = []

    for feature in ref_record.features:
        if feature.type != "CDS":
            continue

        parts = (
            feature.location.parts
            if isinstance(feature.location, CompoundLocation)
            else [feature.location]
        )

        transferred_parts = []
        partial_start = False
        partial_end = False

        for part in parts:
            ref_start = int(part.start) + 1  # BioPython 0-based -> 1-based
            ref_end = int(part.end)           # BioPython end is exclusive, so this is the last position

            # Map reference positions through alignment to query positions
            aln_col_start = ref_coord_map.get(ref_start)
            aln_col_end = ref_coord_map.get(ref_end)

            if aln_col_start is None or aln_col_end is None:
                continue

            query_start = query_rev_map.get(aln_col_start)
            query_end = query_rev_map.get(aln_col_end)

            # Handle positions that fall in gaps: scan for nearest mapped position
            if query_start is None:
                query_start = _find_nearest_mapped_position(
                    query_rev_map, aln_col_start, direction=1
                )
                partial_start = True
            if query_end is None:
                query_end = _find_nearest_mapped_position(
                    query_rev_map, aln_col_end, direction=-1
                )
                partial_end = True

            if query_start is None or query_end is None:
                continue

            # Check for truncation at sequence edges
            if query_start <= 1:
                partial_start = True
            if query_end >= query_ungapped_len:
                partial_end = True

            # Build FeatureLocation (back to 0-based for BioPython)
            start_pos = (
                BeforePosition(query_start - 1)
                if partial_start
                else query_start - 1
            )
            end_pos = (
                AfterPosition(query_end) if partial_end else query_end
            )
            transferred_parts.append(
                FeatureLocation(start_pos, end_pos, strand=part.strand)
            )

        if not transferred_parts:
            continue

        if len(transferred_parts) == 1:
            location = transferred_parts[0]
        else:
            location = CompoundLocation(transferred_parts)

        new_feature = SeqFeature(
            location=location,
            type="CDS",
            qualifiers=dict(feature.qualifiers),
        )
        features.append(new_feature)

    return features


def _find_nearest_mapped_position(rev_map, col, direction=1, max_search=50):
    """Scan alignment columns in `direction` to find a mapped query position."""
    for offset in range(1, max_search + 1):
        candidate = col + (offset * direction)
        if candidate in rev_map:
            return rev_map[candidate]
    return None


# ---------------------------------------------------------------------------
# GenBank file creation
# ---------------------------------------------------------------------------


def create_genbank_record(accession, seq_str, features, description=""):
    """Create a BioPython SeqRecord with GenBank features."""
    record = SeqRecord(
        Seq(seq_str),
        id=accession,
        name=accession,
        description=description,
        annotations={"molecule_type": "DNA"},
    )
    record.features = features
    return record


def create_genbank_files_for_segment(
    alignment_fasta,
    ref_genbank,
    ref_accession,
    output_dir,
):
    """Create annotated GenBank files for all sequences in a segment alignment.

    Args:
        alignment_fasta: path to MAFFT alignment (ref is first sequence)
        ref_genbank: path to NC_ GenBank file
        ref_accession: accession of the NC_ reference
        output_dir: directory for output .gb files
    """
    os.makedirs(output_dir, exist_ok=True)

    records = list(SeqIO.parse(alignment_fasta, "fasta"))
    if not records:
        return []

    # Reference is the first record in the alignment
    ref_record = records[0]
    ref_aligned_str = str(ref_record.seq)

    created_files = []
    for record in records:
        accession = record.id
        aligned_str = str(record.seq)
        ungapped_str = aligned_str.replace("-", "")

        if accession == ref_accession:
            # For the reference itself, just copy the original GenBank
            import shutil
            src = ref_genbank
            dst = os.path.join(output_dir, f"{accession}.gb")
            shutil.copy2(src, dst)
            created_files.append(dst)
            continue

        features = transfer_cds_features(
            ref_genbank, ref_aligned_str, aligned_str, ungapped_str
        )
        gb_record = create_genbank_record(
            accession, ungapped_str, features, description=record.description
        )

        out_path = os.path.join(output_dir, f"{accession}.gb")
        SeqIO.write(gb_record, out_path, "genbank")
        created_files.append(out_path)

    return created_files


# ---------------------------------------------------------------------------
# Metadata parsing
# ---------------------------------------------------------------------------

# Common country/region normalization for influenza strain names
COUNTRY_ALIASES = {
    "hong kong": "Hong Kong",
    "inner mongolia": "China",
    "guangdong": "China",
    "hunan": "China",
    "yunnan": "China",
    "jiangsu": "China",
    "hebei": "China",
    "qinghai": "China",
    "shandong": "China",
    "texas": "USA",
    "california": "USA",
    "ohio": "USA",
    "michigan": "USA",
    "colorado": "USA",
    "kansas": "USA",
    "iowa": "USA",
    "minnesota": "USA",
    "wisconsin": "USA",
    "south dakota": "USA",
    "north dakota": "USA",
    "idaho": "USA",
    "wyoming": "USA",
    "montana": "USA",
    "oregon": "USA",
    "washington": "USA",
    "alaska": "USA",
    "pennsylvania": "USA",
    "new york": "USA",
    "north carolina": "USA",
    "virginia": "USA",
    "england": "UK",
    "scotland": "UK",
}


def parse_strain_metadata(strain_name):
    """Extract host, country, and year from a strain name.

    Strain format: A/{host}/{location}/{id}/{year}
    Returns dict with keys: host, country, year (int or None).
    """
    parts = strain_name.split("/")
    if len(parts) < 4:
        return {"host": None, "country": None, "year": None}

    host = parts[1] if len(parts) > 1 else None
    location = parts[2] if len(parts) > 2 else None

    # Year is typically the last numeric component
    year = None
    year_str = parts[-1]
    if re.match(r"^\d{4}$", year_str):
        year = int(year_str)
    elif re.match(r"^\d{2}$", year_str):
        y = int(year_str)
        year = 2000 + y if y < 50 else 1900 + y

    # Normalize location to country
    country = location
    if location:
        loc_lower = location.lower().strip()
        if loc_lower in COUNTRY_ALIASES:
            country = COUNTRY_ALIASES[loc_lower]

    return {"host": host, "country": country, "year": year}


def assign_year_bin(year, bin_size=2):
    """Assign a year to a bin. E.g., bin_size=2: 2024 -> '2024-2025'."""
    if year is None:
        return "unknown"
    bin_start = year - (year % bin_size)
    bin_end = bin_start + bin_size - 1
    return f"{bin_start}-{bin_end}"


def build_metadata_groups(
    segment_fasta,
    year_bin_size=2,
    genoflu_results=None,
    segment_name=None,
):
    """Group sequences by (country, year_bin) with optional lineage.

    Args:
        segment_fasta: path to per-segment FASTA
        year_bin_size: number of years per bin
        genoflu_results: path to GenoFLU results.tsv (optional)
        segment_name: segment name for lineage lookup (optional)

    Returns:
        groups: {group_key: [accession_id, ...]}
        metadata: {accession_id: {host, country, year, year_bin, ...}}
    """
    # Load GenoFLU lineage data if available
    lineage_map = {}
    if genoflu_results and os.path.exists(genoflu_results):
        lineage_map = _load_genoflu_lineages(genoflu_results, segment_name)

    groups = defaultdict(list)
    metadata = {}

    for record in SeqIO.parse(segment_fasta, "fasta"):
        acc = record.id
        strain = extract_strain_name(record.description)
        if strain:
            meta = parse_strain_metadata(strain)
        else:
            meta = {"host": None, "country": None, "year": None}

        meta["year_bin"] = assign_year_bin(meta["year"], year_bin_size)
        meta["lineage"] = lineage_map.get(acc)
        metadata[acc] = meta

        # Build group key
        country = meta["country"] or "unknown"
        year_bin = meta["year_bin"]
        if meta["lineage"]:
            group_key = (country, year_bin, meta["lineage"])
        else:
            group_key = (country, year_bin)

        groups[group_key].append(acc)

    return dict(groups), metadata


def _load_genoflu_lineages(results_tsv, segment_name=None):
    """Load per-accession lineage assignments from GenoFLU results."""
    lineage_map = {}
    with open(results_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = row.get("name", row.get("sample", ""))
            # GenoFLU has per-segment lineage columns
            if segment_name:
                lineage = row.get(segment_name, row.get("genotype", ""))
            else:
                lineage = row.get("genotype", "")
            if acc and lineage:
                lineage_map[acc] = lineage
    return lineage_map


# ---------------------------------------------------------------------------
# TN93 clustering
# ---------------------------------------------------------------------------


def run_tn93_cluster(
    input_fasta,
    output_fasta,
    threshold=0.015,
    work_dir=None,
):
    """Run TN93 distance calculation and clustering on a FASTA file.

    Keeps cluster centroids (representative sequences).
    Returns list of centroid accession IDs.
    """
    if work_dir is None:
        work_dir = tempfile.mkdtemp()
    os.makedirs(work_dir, exist_ok=True)

    distances_csv = os.path.join(work_dir, "distances.csv")
    cluster_json = os.path.join(work_dir, "clusters.json")

    # Count sequences
    n_seqs = sum(1 for _ in SeqIO.parse(input_fasta, "fasta"))
    if n_seqs <= 1:
        # Single sequence: just copy
        records = list(SeqIO.parse(input_fasta, "fasta"))
        if records:
            SeqIO.write(records, output_fasta, "fasta")
        return [r.id for r in records]

    # Compute pairwise TN93 distances
    subprocess.run(
        [
            "tn93",
            "-t", str(threshold),
            "-o", distances_csv,
            input_fasta,
        ],
        check=True,
        capture_output=True,
    )

    # Cluster using tn93-cluster
    with open(distances_csv) as infile, open(cluster_json, "w") as outfile:
        subprocess.run(
            ["tn93-cluster", "-t", str(threshold)],
            stdin=infile,
            stdout=outfile,
            check=True,
            capture_output=False,
        )

    # Parse cluster output to find centroids
    centroid_ids = _parse_tn93_cluster_centroids(cluster_json, input_fasta)

    # Extract centroid sequences
    centroid_set = set(centroid_ids)
    records = [
        r for r in SeqIO.parse(input_fasta, "fasta") if r.id in centroid_set
    ]
    SeqIO.write(records, output_fasta, "fasta")
    return centroid_ids


def _parse_tn93_cluster_centroids(cluster_json, input_fasta):
    """Parse tn93-cluster output to identify cluster centroids.

    tn93-cluster outputs a JSON mapping of sequence -> cluster_id.
    The centroid of each cluster is the first sequence assigned to that cluster.
    """
    with open(cluster_json) as f:
        content = f.read().strip()

    if not content:
        # Fallback: return all sequences
        return [r.id for r in SeqIO.parse(input_fasta, "fasta")]

    cluster_assignments = json.loads(content)

    # Group sequences by cluster
    clusters = defaultdict(list)
    for seq_id, cluster_id in cluster_assignments.items():
        clusters[cluster_id].append(seq_id)

    # Centroid = first member of each cluster (alphabetical as tie-breaker)
    centroids = []
    for cluster_id in sorted(clusters.keys()):
        members = clusters[cluster_id]
        centroids.append(members[0])

    return centroids


def cluster_by_metadata_groups(
    segment_fasta,
    groups,
    output_fasta,
    tn93_threshold=0.015,
    min_seqs_per_group=1,
    work_dir=None,
):
    """Run TN93 clustering within each metadata group, then combine centroids.

    Args:
        segment_fasta: path to full per-segment FASTA
        groups: {group_key: [accession_id, ...]}
        output_fasta: path to write combined centroid FASTA
        tn93_threshold: TN93 distance threshold
        min_seqs_per_group: keep at least this many seqs per group
        work_dir: scratch directory for intermediate files
    """
    if work_dir is None:
        work_dir = tempfile.mkdtemp()
    os.makedirs(work_dir, exist_ok=True)

    # Index all sequences
    seq_index = SeqIO.index(segment_fasta, "fasta")

    all_centroids = []
    for group_idx, (group_key, accessions) in enumerate(sorted(groups.items())):
        group_dir = os.path.join(work_dir, f"group_{group_idx}")
        os.makedirs(group_dir, exist_ok=True)

        # Write group FASTA
        group_fasta = os.path.join(group_dir, "group.fasta")
        group_records = [
            seq_index[acc] for acc in accessions if acc in seq_index
        ]
        if not group_records:
            continue

        SeqIO.write(group_records, group_fasta, "fasta")

        if len(group_records) <= min_seqs_per_group:
            # Keep all sequences in small groups
            all_centroids.extend(group_records)
        else:
            # Cluster and keep centroids
            clustered_fasta = os.path.join(group_dir, "centroids.fasta")
            try:
                run_tn93_cluster(
                    group_fasta,
                    clustered_fasta,
                    threshold=tn93_threshold,
                    work_dir=group_dir,
                )
                centroid_records = list(
                    SeqIO.parse(clustered_fasta, "fasta")
                )
                # Ensure minimum representation
                if len(centroid_records) < min_seqs_per_group:
                    all_centroids.extend(group_records[:min_seqs_per_group])
                else:
                    all_centroids.extend(centroid_records)
            except (subprocess.CalledProcessError, FileNotFoundError):
                # If tn93 fails, keep all sequences from this group
                all_centroids.extend(group_records)

    seq_index.close()

    SeqIO.write(all_centroids, output_fasta, "fasta")
    return len(all_centroids)


# ---------------------------------------------------------------------------
# Database assembly
# ---------------------------------------------------------------------------


def assemble_segment_database(
    clustered_fasta,
    nc_fasta,
    nc_genbank,
    nc_accession,
    genbanks_dir,
    output_dir,
):
    """Assemble final database for one segment matching the parent pipeline contract.

    Output structure:
        output_dir/
            all.fasta          - all sequences for VAPOR
            metadata.gb        - NC_ reference GenBank (default fallback)
            sequence.fasta     - NC_ reference FASTA
            genbanks/          - per-sequence annotated GenBanks
                {accession}.gb
    """
    os.makedirs(output_dir, exist_ok=True)
    genbanks_out = os.path.join(output_dir, "genbanks")
    os.makedirs(genbanks_out, exist_ok=True)

    # Build all.fasta: clustered sequences + NC_ reference
    all_records = []
    seen_ids = set()

    # Add NC_ reference first
    for record in SeqIO.parse(nc_fasta, "fasta"):
        all_records.append(record)
        seen_ids.add(record.id)

    # Add clustered sequences
    for record in SeqIO.parse(clustered_fasta, "fasta"):
        if record.id not in seen_ids:
            all_records.append(record)
            seen_ids.add(record.id)

    all_fasta_path = os.path.join(output_dir, "all.fasta")
    SeqIO.write(all_records, all_fasta_path, "fasta")

    # Copy NC_ reference files
    import shutil
    shutil.copy2(nc_genbank, os.path.join(output_dir, "metadata.gb"))
    shutil.copy2(nc_fasta, os.path.join(output_dir, "sequence.fasta"))

    # Copy GenBank files for sequences in all.fasta
    for record in all_records:
        src_gb = os.path.join(genbanks_dir, f"{record.id}.gb")
        dst_gb = os.path.join(genbanks_out, f"{record.id}.gb")
        if os.path.exists(src_gb) and not os.path.exists(dst_gb):
            shutil.copy2(src_gb, dst_gb)

    return len(all_records)
