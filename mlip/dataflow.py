# Native Python module imports for environment check
import os
import sys
import shutil


def perform_initial_environment_check():
    """
    Checks for essential CLI tools and Python modules.
    Prints a simple success message if all are found.
    Prints detailed errors and exits script if any dependency is missing.
    """

    error_messages = []
    found_all_cli_tools = True
    found_all_python_modules = True

    required_cli_tools = [
        "bedtools",
        "bcftools",
        "bowtie2",
        "multiqc",
        "samtools",
        "trimmomatic",
        "varscan",
        "efetch",
        "seqkit",
        "sed",
        "csvcut",
        "perbase",
        "vapor.py",
        "mafft",
        "ivar",
    ]
    for tool in required_cli_tools:
        tool_path = shutil.which(tool)
        if not tool_path:
            error_messages.append(f"  ❌ ERROR: CLI tool '{tool}' not found in PATH.")
            found_all_cli_tools = False

    required_python_modules = {
        "Bio": "BioPython",
        "pandas": "Pandas",
        "yaml": "PyYAML",  # Essential for config loading
        "snakemake": "Snakemake",  # Crucial for the main pipeline
        "altair": "Altair",  # For visualization
    }
    for module_name, display_name in required_python_modules.items():
        try:
            __import__(module_name)
        except ModuleNotFoundError:
            error_messages.append(
                f"  ❌ ERROR: Python module {display_name} ({module_name}) not found."
            )
            found_all_python_modules = False
        except ImportError as e:
            # Handle cases where the module is found, but fails during import
            error_messages.append(
                f"  ❌ ERROR: Failed importing Python module {display_name} ({module_name})."
            )
            error_messages.append(
                f"      └─> This often means an issue with its dependencies or environment setup."
            )
            error_messages.append(
                f"      └─> Original Error: {e}"
            )  # Include the original error message
            found_all_python_modules = False
        except Exception as e:
            # Catch any other totally unexpected error during import attempt
            error_messages.append(
                f"  ❌ ERROR: Unexpected error while trying to import {display_name} ({module_name})."
            )
            error_messages.append(f"      └─> Original Error: {e}")
            found_all_python_modules = False

    # Print detailed errors and exit OR print simple success message
    if not found_all_cli_tools or not found_all_python_modules:
        # Print errors only if there are errors
        print("\n--- Initial Environment Check FAILED ---")  # Clear failure indication
        print("Required dependencies missing:")
        for msg in error_messages:
            print(msg)

        # Print the detailed alert message
        print("\n---------------------------------------------------------------")
        print("❌ CRITICAL ERROR: Missing one or more essential dependencies.")
        print(
            "   Your environment does not appear to be correctly installed or activated."
        )
        print("\n   Please ensure you have followed the installation instructions:")
        print("   and activated the Conda environment, e.g.:\n")
        print("         conda activate mlip\n")
        print("     Please also make sure that when you run the command:\n")
        print("         which python\n")
        print("     you see something like:\n")
        print("         /path/to/anaconda/envs/mlip/python\n")
        print(
            "     to ensure that your shell is utilizing python from the MLIP environment."
        )
        print("\n   If the issue persists, you may need to:")
        print("    1) reinstall the environment using the command from the README.")
        print("    2) assess the way your shell is currently configured.\n")
        print("Option 2 can be technically challenging. Know what shell you are")
        print("using (e.g., zsh vs bash). Know what startup scripts run (e.g")
        print(".zshrc), how they modify your environment, and why. It is possible")
        print("that this error occurs due personal customizations of your")
        print("environment that are difficult/impossible to account for in advance.")
        print("---------------------------------------------------------------")
        sys.exit(1)  # Exit the script
    else:
        # All checks passed, print the single, simple success message
        print("✅ MLIP environment detected.")
        # No verbose output needed on success


perform_initial_environment_check()

import csv
import gzip
import re
import json
from pathlib import Path
import argparse
import subprocess
from collections import defaultdict
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation
import pandas as pd
import numpy as np
import yaml


def load_reference_dictionary(reference, using_zip=False):
    if using_zip:
        return {}
    reference_dictionary = {}
    with open("references.tsv", "r", newline="") as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter="\t")
        for row in reader:
            if row["reference"] == reference:
                reference_dictionary[row["segment_key"]] = row
    return reference_dictionary


def load_manifest(analysis_dir):
    with open(f"{analysis_dir}/file_manifest.json") as manifest_file:
        manifest = json.load(manifest_file)
    return manifest


def load_segments(config, analysis_dir):
    ref = config.get("reference", "")
    if ref.endswith(".zip"):
        # After unzip, expect one folder per segment under analysis_dir/reference/
        base = Path(analysis_dir) / "reference"
        return [
            p.name for p in base.iterdir() if p.is_dir() and not p.name.startswith(".")
        ]
    else:
        # Treat as a named key—load the dictionary here
        ref_dict = load_reference_dictionary(ref)
        return list(ref_dict.keys())


def extract_genes(input_gtf, output_gene_list):
    pattern = re.compile(r'gene_id\s+"([^"]+)"')
    with open(input_gtf) as input_file, open(output_gene_list, "w") as output_file:
        for line in input_file:
            match = pattern.search(line)
            if match:
                output_file.write(match.group(1) + "\n")


def check_duplicates(lines):
    seen = {}
    saw_dupes = False
    for line_num, line in enumerate(lines, 1):
        id_ = line.strip()
        if id_ in seen:
            saw_dupes = True
            print(f"Duplicate ID '{id_}' found on lines {seen[id_]} and {line_num}")
        else:
            seen[id_] = line_num
    if saw_dupes:
        print(
            "You have duplicate sequencing IDs! MLIP cannot proceed until this is rectified."
        )
        sys.exit(1)


def preprocess(id_filepath, seq_key="Seq"):
    with open(id_filepath, "r") as f:
        lines = f.read().splitlines()
    check_duplicates(lines)

    # Load config to get analysis directory
    config = load_mlip_config()
    analysis_dir = config.get('analysis', '')
    if not analysis_dir:
        print("ERROR: 'analysis' key not found in config.yml.")
        print("Please add an analysis name to your config.yml.")
        sys.exit(1)

    sorted_seq_ids = sorted(lines, key=lambda x: x.lower())
    seq_key_pattern = re.compile(rf"_{seq_key}(\d+)", re.IGNORECASE)
    key_hash = {}

    fieldnames = ["SequencingId", "SampleId", "Replicate"]
    os.makedirs(analysis_dir, exist_ok=True)
    f = open(f"{analysis_dir}/metadata.tsv", "w", newline="")
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    # process keys
    for seq_id in sorted_seq_ids:
        match = seq_key_pattern.search(seq_id)
        found_match = False
        if match:
            sample_id = seq_id.replace(match.group(0), "").lower()
            found_match = True
        else:
            sample_id = ""
        writer.writerow(
            {"SequencingId": seq_id, "SampleId": sample_id, "Replicate": ""}
        )

    f.close()
    print(f"Metadata spreadsheet written to {analysis_dir}/metadata.tsv.")
    print("Please edit, then run the flow step.")
    return


def preprocess_cli(args):
    if args.key:
        preprocess(args.file, args.key)
    else:
        preprocess(args.file)


def tokenize(identifier):
    return identifier.lower().replace("_", "").replace("-", "")


def is_forward_fastq_path(fastq_path):
    return fastq_path.split("_")[-2] == "R1"


def convert_forward_to_reverse(fastq_path):
    dirname, filename = os.path.split(fastq_path)
    new_filename = filename[::-1].replace("1R", "2R", 1)[::-1]  # Reverse swap
    return os.path.join(dirname, new_filename)


def fastq_is_low_coverage(filepath, min_reads=50):
    try:
        with gzip.open(filepath, "rt") as f:
            lines = 0
            for line in f:
                lines += 1
                if lines > min_reads * 4:
                    return False
            return (lines // 4) < min_reads
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return True


def flow(args):
    # manifest_samples_data will store: {sample_id: {replicate_num: [exp_details_list]}}
    manifest_samples_data = defaultdict(lambda: defaultdict(list))
    config = load_mlip_config()

    analysis_dir = config.get('analysis', '')
    if not analysis_dir:
        print("ERROR: 'analysis' key not found in config.yml.")
        sys.exit(1)

    with open(f"{analysis_dir}/metadata.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    data_root = Path(config["data_root_directory"]).expanduser()

    fastq_paths = [
        str(fastq_path)
        for fastq_path in data_root.rglob("*.fastq.gz")
        if is_forward_fastq_path(str(fastq_path))
    ]
    fastq_tokens = [
        tokenize(os.path.basename(fastq_path)) for fastq_path in fastq_paths
    ]
    fastq_hash = {
        token: {
            "path": path,
            "reverse": convert_forward_to_reverse(path),
            "seen": False,
        }
        for token, path in zip(fastq_tokens, fastq_paths)
    }

    for row in rows:
        sample_id = row["SampleId"]
        sequencing_id_from_meta = row["SequencingId"]
        replicate_from_meta = row["Replicate"]

        if not sample_id or not sequencing_id_from_meta or not replicate_from_meta:
            continue

        sequencing_token_from_meta = tokenize(sequencing_id_from_meta)
        token_pattern = re.compile(rf"^{sequencing_token_from_meta}s\d+l\d+")

        found_match_for_row = False
        for fastq_file_token, file_data in fastq_hash.items():
            if not file_data["seen"] and token_pattern.search(fastq_file_token):
                old_forward_path = file_data["path"]
                old_reverse_path = file_data["reverse"]

                if (
                    not Path(old_forward_path).exists()
                    or not Path(old_reverse_path).exists()
                ):
                    continue

                is_low_cov = fastq_is_low_coverage(
                    old_forward_path
                ) or fastq_is_low_coverage(old_reverse_path)

                experiment_entry = {
                    "sequencing_id_from_metadata": sequencing_id_from_meta,
                    "source_forward_path": str(Path(old_forward_path).resolve()),
                    "source_reverse_path": str(Path(old_reverse_path).resolve()),
                    "source_gzipped": True,
                    "is_low_coverage": is_low_cov,
                }
                manifest_samples_data[sample_id][replicate_from_meta].append(
                    experiment_entry
                )

                file_data["seen"] = True
                found_match_for_row = True
                break

    # Convert defaultdicts to dicts for cleaner JSON output
    final_samples_data = {
        s_id: dict(replicates) for s_id, replicates in manifest_samples_data.items()
    }

    final_manifest = {
        "manifest_details": {
            "manifest_mode_used": "basespace",
            "data_root_queried": str(data_root.resolve()),
        },
        "samples": final_samples_data,
    }

    if not final_samples_data:  # Check if any samples had processable experiments
        sys.stderr.write(
            "\nERROR: No FASTQ files were found or matched for any samples.\n"
            f"       Please check '{analysis_dir}/metadata.tsv', 'data_root_directory' in 'config.yml' ({data_root}),\n"
            f"       and ensure files follow the BaseSpace convention.\n\n"
            f"       Also, you are in default mode... did you mean to run in SRA mode?"
        )
        sys.exit(1)

    os.makedirs(analysis_dir, exist_ok=True)
    with open(f"{analysis_dir}/file_manifest.json", "w") as f_out:
        json.dump(final_manifest, f_out, indent=2)


def fastq_is_low_coverage_sra_generic(filepath_str, min_reads=50):
    line_count = 0
    try:
        with gzip.open(filepath_str, "rt") as f_gz:
            for line in f_gz:
                line_count += 1
                if line_count > min_reads * 4:
                    return False
            return (line_count // 4) < min_reads
    except gzip.BadGzipFile:  # Not a gzip file, try as plain text
        try:
            with open(filepath_str, "r") as f_plain:
                for line in f_plain:
                    line_count += 1
                    if line_count > min_reads * 4:
                        return False
                return (line_count // 4) < min_reads
        except Exception:
            return True  # Error reading plain, assume low
    except FileNotFoundError:
        return True  # File not found, assume low
    except Exception:
        return True  # Other gzip errors, assume low


def sra_flow(args):
    manifest_samples_data = defaultdict(lambda: defaultdict(list))
    config = load_mlip_config()

    analysis_dir = config.get('analysis', '')
    if not analysis_dir:
        print("ERROR: 'analysis' key not found in config.yml.")
        sys.exit(1)

    data_root = Path(config["data_root_directory"]).expanduser()

    with open(f"{analysis_dir}/metadata.tsv", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    for row in rows:
        sample_id = row["SampleId"]
        # For this mode, SequencingId from metadata IS the base filename part
        sequencing_id_base = row["SequencingId"]
        replicate_from_meta = row["Replicate"]

        if not sample_id or not sequencing_id_base or not replicate_from_meta:
            continue

        fwd_path_to_use = None
        rev_path_to_use = None
        is_gzipped_source = False

        # Check for gzipped pair first, then plain fastq
        fwd_gz = data_root / f"{sequencing_id_base}_1.fastq.gz"
        rev_gz = data_root / f"{sequencing_id_base}_2.fastq.gz"
        fwd_plain = data_root / f"{sequencing_id_base}_1.fastq"
        rev_plain = data_root / f"{sequencing_id_base}_2.fastq"

        if fwd_gz.exists() and rev_gz.exists():
            fwd_path_to_use = str(fwd_gz.resolve())
            rev_path_to_use = str(rev_gz.resolve())
            is_gzipped_source = True
        elif fwd_plain.exists() and rev_plain.exists():
            fwd_path_to_use = str(fwd_plain.resolve())
            rev_path_to_use = str(rev_plain.resolve())
            is_gzipped_source = False
        else:
            # print(f"SRA/Generic Mode: File pair for {sequencing_id_base} not found.")
            continue  # Skip if pair not found

        is_low_cov = fastq_is_low_coverage_sra_generic(
            fwd_path_to_use
        ) or fastq_is_low_coverage_sra_generic(rev_path_to_use)

        experiment_entry = {
            "sequencing_id_from_metadata": sequencing_id_base,  # Storing the ID from metadata
            "source_forward_path": fwd_path_to_use,
            "source_reverse_path": rev_path_to_use,
            "source_gzipped": is_gzipped_source,
            "is_low_coverage": is_low_cov,
        }
        manifest_samples_data[sample_id][replicate_from_meta].append(experiment_entry)

    final_samples_data = {
        s_id: dict(replicates) for s_id, replicates in manifest_samples_data.items()
    }

    final_manifest = {
        "manifest_details": {
            "manifest_mode_used": "sra_generic",  # Indicate mode used
            "data_root_queried": str(data_root.resolve()),
        },
        "samples": final_samples_data,
    }

    if not final_samples_data:  # Check if any samples had processable experiments
        sys.stderr.write(
            "\nERROR: No FASTQ files were found or matched for any samples.\n"
            f"       Please check '{analysis_dir}/metadata.tsv', 'data_root_directory' in 'config.yml' ({data_root}),\n"
            f"       and ensure files follow the '{Path('{SequencingId}')}_1.fastq[.gz]' convention.\n\n"
        )
        sys.exit(1)

    os.makedirs(analysis_dir, exist_ok=True)
    with open(f"{analysis_dir}/file_manifest.json", "w") as f_out:
        json.dump(final_manifest, f_out, indent=2)

    total_experiments_in_manifest = sum(
        len(experiments)
        for reps in final_samples_data.values()
        for experiments in reps.values()
    )
    print(
        f"SRA/Generic mode: File manifest generated at {analysis_dir}/file_manifest.json with {total_experiments_in_manifest} total experiments for {len(final_samples_data)} samples."
    )


def _prepare_reference_from_zip_if_needed(config, analysis_dir):
    ref = Path(config.get("reference", "")).expanduser()
    marker = Path(analysis_dir) / "reference" / ".unzipped"

    if not ref.suffix == ".zip":
        print("Reference will be downloaded from Genbank")
        return
    if not ref.is_file():
        sys.exit(f"ERROR: Reference ZIP not found: {ref}")

    ref_dir = Path(analysis_dir) / "reference"
    shutil.rmtree(ref_dir, ignore_errors=True)
    subprocess.run(["unzip", "-o", str(ref), "-d", analysis_dir], check=True)
    marker.touch()


def flow_cli(args):
    config = load_mlip_config()
    analysis_dir = config.get('analysis', '')
    if not analysis_dir:
        print("ERROR: 'analysis' key not found in config.yml.")
        sys.exit(1)
    _prepare_reference_from_zip_if_needed(config, analysis_dir)
    if args.sra_mode:
        print("Initiating SRA/Generic mode manifest generation...")
        sra_flow(args)
    else:
        print("Initiating BaseSpace mode manifest generation...")
        flow(args)  # Call the renamed BaseSpace function


def concatenate_replicates_from_manifest_py(
    manifest_filepath,
    sample_id,
    replicate_num_str,
    output_forward_fastq_path,
    output_reverse_fastq_path,
):
    with open(manifest_filepath, "r") as f_manifest:
        manifest_data = json.load(f_manifest)

    sample_specific_data = manifest_data.get("samples", {}).get(sample_id, {})
    experiments_for_replicate = sample_specific_data.get(replicate_num_str, [])

    valid_experiments = [
        exp
        for exp in experiments_for_replicate
        if not exp.get("is_low_coverage", False)
    ]

    # Ensure output directory exists (minimal check)
    os.makedirs(os.path.dirname(output_forward_fastq_path), exist_ok=True)

    if not valid_experiments:
        Path(output_forward_fastq_path).touch()
        Path(output_reverse_fastq_path).touch()
        return

    with open(output_forward_fastq_path, "wb") as fwd_out_handle, open(
        output_reverse_fastq_path, "wb"
    ) as rev_out_handle:

        for exp_data in valid_experiments:
            fwd_src_path = exp_data["source_forward_path"]
            rev_src_path = exp_data["source_reverse_path"]
            is_gzipped = exp_data["source_gzipped"]

            if is_gzipped:
                with gzip.open(fwd_src_path, "rb") as f_in:
                    shutil.copyfileobj(f_in, fwd_out_handle)
            else:
                with open(fwd_src_path, "rb") as f_in:
                    shutil.copyfileobj(f_in, fwd_out_handle)

            if is_gzipped:
                with gzip.open(rev_src_path, "rb") as f_in:
                    shutil.copyfileobj(f_in, rev_out_handle)
            else:
                with open(rev_src_path, "rb") as f_in:
                    shutil.copyfileobj(f_in, rev_out_handle)


def load_mlip_config():
    config_filename = "config.yml"
    config_path = Path(config_filename)  # Using pathlib

    if not config_path.exists():
        raise FileNotFoundError()
    try:
        with open(config_path, "r") as f:
            config_data = yaml.safe_load(f)

        return config_data
    except:
        raise Exception("Could not load configuration file!")


def command_line_interface():
    # --- Main Parser Setup ---
    parser = argparse.ArgumentParser(
        # Description appears at the top, explaining the script's overall purpose.
        description=(
            "=====================================================\n"
            " MLIP Dataflow & Setup Tool\n"
            "=====================================================\n"
            "This script manages the initial data setup steps for\n"
            "this viral deep sequencing pipeline. It helps you:\n"
            "  - Check your environment and configuration.\n"
            "  - Create a metadata file from sequencing IDs.\n"
            "  - Arrange your FASTQ files for ingestion by the pipeline."
        ),
        # Epilog appears at the very bottom, after all arguments.
        epilog=(
            "-----------------------------------------------------\n"
            "Typical Workflow Steps:\n"
            " 1. `check`: Verify your setup (run this as much as you want!).\n"
            "    Usage: python mlip/dataflow.py check\n"
            " 2. `preprocess`: Create metadata sheet from sequence IDs.\n"
            "    Usage: python mlip/dataflow.py preprocess -f <id_file>\n"
            " 3. *Manually Edit* the generated `data/metadata.tsv`.\n"
            " 4. `flow`: Move FASTQ files based on completed metadata.\n"
            "    Usage: python mlip/dataflow.py flow\n\n"
            "For detailed help on any command and its specific options:\n"
            "  python mlip/dataflow.py <command> -h \n"
            "e.g.:\n"
            "  python mlip/dataflow.py preprocess -h\n"
            "-----------------------------------------------------"
        ),
        # This formatter preserves your line breaks in description and epilog.
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- Subparsers Setup ---
    # Give the group of commands a title and description
    subparsers = parser.add_subparsers(
        dest="command",  # Stores the chosen command name
        required=True,  # A command MUST be provided
        title="Available Commands",
        description="Choose one of the following commands to perform a specific task:",
        metavar="<command>",  # How the placeholder appears in the usage line
    )

    # --- Check Subcommand Parser ---
    check_parser = subparsers.add_parser(
        "check",
        help="✅ Verify environment, config, and data setup status.",  # Concise help for the list
        description=(  # More detailed help shown with 'check -h'
            "Performs checks on your MLIP setup:\n"
            " - Verifies required software and Python packages are installed.\n"
            " - Checks if `config.yml` exists and is valid.\n"
            " - Looks for `data/metadata.tsv` and assesses its status.\n"
            " - Checks if data appears to have been moved by the `flow` command.\n"
            "Run this first or if you encounter problems."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    check_parser.set_defaults(
        func=check_cli
    )  # No arguments specific to 'check' needed yet

    # --- Preprocess Subcommand Parser ---
    pre_parser = subparsers.add_parser(
        "preprocess",
        help="📄 Create initial metadata sheet from sequencing IDs.",
        description=(  # More detailed help shown with 'preprocess -h'
            "Reads a simple text file containing one Sequencing ID per line \n"
            "(e.g., 'SampleA_Seq1', matching your FASTQ file names from BaseSpace/SRA).\n"
            "It automatically creates a template spreadsheet at `data/metadata.tsv`,\n"
            "attempting to guess 'SampleId' based on common patterns.\n\n"
            "--> IMPORTANT: You MUST manually open and edit `data/metadata.tsv` \n"
            "    after running this command to verify/correct 'SampleId' and \n"
            "    fill in the 'Replicate' column before running the `flow` command."
        ),
        epilog="Example: python mlip/dataflow.py preprocess -f ./my_sequence_ids.txt",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    pre_parser.add_argument(
        "-f",
        "--file",
        required=True,
        type=str,
        metavar="<path/to/id_list.txt>",  # More descriptive placeholder
        help="REQUIRED: Path to the input text file containing Sequencing IDs (one per line).",
    )
    pre_parser.add_argument(
        "-k",
        "--key",
        required=False,
        type=str,
        default="Seq",
        metavar="<key>",
        help="Keyword in Sequencing ID that denotes the sequencing run (default: 'Seq'). Used for auto-generating SampleId (e.g., 'SampleA_Seq1' -> SampleId 'SampleA').",
    )
    pre_parser.set_defaults(func=preprocess_cli)

    # --- Flow Subcommand Parser ---
    flow_parser = subparsers.add_parser(
        "flow",
        help="🚚 Move FASTQ files based on the completed metadata sheet.",
        description=(  # More detailed help shown with 'flow -h'
            "Reads your completed `data/metadata.tsv` spreadsheet AND your `config.yml` file.\n"
            "It finds the corresponding FASTQ files (R1/R2 pairs) within the \n"
            "`data_root_directory` specified in your config, and copies them into \n"
            "an organized structure within the `data/` directory (e.g., data/SampleA/sequencing-1/).\n"
            "This prepares the data for the main Snakemake pipeline.\n\n"
            "--> Ensure `config.yml` points to the correct `data_root_directory` \n"
            "    and `data/metadata.tsv` is fully edited and saved before running."
        ),
        epilog="Example: python mlip/dataflow.py flow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    flow_parser.add_argument(
        "--sra-mode", action="store_true", help="Enable SRA mode for data flowing."
    )
    flow_parser.set_defaults(func=flow_cli)

    # --- Argument Parsing ---
    # If the script is called with no arguments other than the script name itself
    # (e.g., "python3 mlip/dataflow.py" with nothing after it)
    if len(sys.argv) == 1:  # sys.argv[0] is the script name
        parser.print_help()  # Display the full help message
        sys.exit(0)  # Exit gracefully

    # Only parse args if the script is run directly (standard practice)
    if __name__ == "__main__":
        args = parser.parse_args()
        # Execute the function associated with the chosen subcommand
        if hasattr(args, "func"):
            args.func(args)
        else:
            # Should not happen if 'required=True' is set for subparsers,
            # but good practice to handle case where no command is given
            # if required=False were used.
            parser.print_help()


# --- Helper functions for printing status ---
def print_section_header(title):
    print(f"\n--- {title} ---")


def print_status_item(message, level="info"):
    # Simple print; can be enhanced with colors if desired using external libs
    prefix_map = {
        "info": "ℹ️ ",
        "success": "✅ ",
        "warning": "⚠️ ",
        "error": "❌ ",
        "bell": "🔔 ",
    }
    print(f"  {prefix_map.get(level.lower(), 'ℹ️ ')}{message}")


def print_guidance(message):
    print(f"       {message}")


# --- CLI handler for 'check' and its worker ---
def check_cli(args):
    """
    Handles the 'check' subcommand.
    Orchestrates the detailed pipeline status reporting.
    """
    print("\n🔎 Running MLIP Setup & Data Status Check...")

    all_checks_clear = report_pipeline_status()

    if all_checks_clear:
        print(
            "\n👍 All critical setup items appear to be in order for the next likely step."
        )
        print(
            "   Please review any ℹ️ informational or ⚠️ warning messages above for further context."
        )
        print("   If you haven't already, run the following command:.")
        print("")
        print("      python mlip/dataflow.py flow")
        print("")
        print("   to move data out of your data root directory and into the pipeline.")
        print("   You should then be ready to run Snakemake commands.")
    else:
        print(
            "\n⚠️ Some setup items need attention. Please review the ❌ errors and ⚠️ warnings above."
        )
    print("\n--- Check Complete ---")


def report_pipeline_status():
    """
    Performs and reports on various stages of the pipeline setup.
    Returns True if no critical errors are found, False otherwise.
    """
    overall_status_ok = True
    config = None
    analysis_dir = None

    # --- 1. Configuration File (`config.yml`) ---
    print_section_header("1. Configuration (`config.yml`)")
    try:
        config = load_mlip_config()
        print_status_item("`config.yml` found and loaded successfully.", "success")

        essential_keys = [
            "analysis",
            "reference",
            "data_root_directory",
            "consensus_minimum_coverage",
            "variants_minimum_coverage",
        ]  # Add more as needed
        missing_keys = [k for k in essential_keys if k not in config]
        if missing_keys:
            print_status_item(
                f"`config.yml` is missing essential keys: {', '.join(missing_keys)}",
                "error",
            )
            overall_status_ok = False
        else:
            print_status_item("Essential configuration keys are present.", "success")

        # Get analysis directory
        analysis_dir = config.get('analysis', '')
        if not analysis_dir:
            print_status_item("`analysis` key not specified in `config.yml`.", "error")
            overall_status_ok = False

        data_root_val = config.get("data_root_directory")
        if not data_root_val:
            print_status_item(
                "`data_root_directory` not specified in `config.yml`.", "error"
            )
            overall_status_ok = False
        else:
            data_root = Path(data_root_val).expanduser()
            if not data_root.exists():
                print_status_item(
                    f"`data_root_directory` ('{data_root}') from `config.yml` does not exist.",
                    "error",
                )
                print_guidance(
                    "Please create this directory or correct the path in `config.yml`."
                )
                # This is a warning, might not block all subsequent checks but flow will fail.
            else:
                print_status_item(
                    f"`data_root_directory` ('{data_root}') exists.", "success"
                )

    except FileNotFoundError as e:
        print_status_item(f"`config.yml` not found.", "error")
        print_guidance(
            f"Please copy `config.yml.template` to `config.yml` and edit it."
        )
        overall_status_ok = False
    except (ValueError, yaml.YAMLError) as e:  # Errors from load_mlip_config
        print_status_item(f"Error loading `config.yml`: {e}", "error")
        print_guidance(
            "Please ensure `config.yml` is correctly formatted and not empty."
        )
        overall_status_ok = False
    except Exception as e:
        print_status_item(f"Unexpected error processing `config.yml`: {e}", "error")
        overall_status_ok = False

    # --- 2. Reference Genome Setup ---
    print_section_header("2. Reference Genome Setup")
    if config and "reference" in config and analysis_dir:
        reference_config_val = str(config["reference"])
        data_ref_path = Path(analysis_dir) / "reference"
        unzipped_marker_path = (
            data_ref_path / ".unzipped"
        )  # Matches marker from flow step

        is_zip_ref = reference_config_val.endswith(".zip")

        if is_zip_ref:
            configured_zip_path = Path(reference_config_val).expanduser()
            print_status_item(
                f"Reference mode: ZIP file ('{configured_zip_path}').", "info"
            )

            if not configured_zip_path.is_file():
                print_status_item(
                    f"Configured reference ZIP not found: '{configured_zip_path}'.",
                    "error",
                )
                overall_status_ok = False
            elif not unzipped_marker_path.exists():
                print_status_item(
                    f"Reference ZIP not yet processed by a 'flow' command.", "warning"
                )
                print_guidance(
                    f"Run 'python mlip/dataflow.py flow' to prepare from '{configured_zip_path}'."
                )
            else:  # Marker exists
                print_status_item(
                    f"Reference ZIP appears processed by 'flow'.", "success"
                )
                if data_ref_path.is_dir():
                    segment_dirs = [
                        d.name
                        for d in data_ref_path.iterdir()
                        if d.is_dir() and not d.name.startswith(".")
                    ]
                    if not segment_dirs:
                        print_status_item(
                            f"No segment subdirs in '{data_ref_path}' after 'flow' processing.",
                            "warning",
                        )
                    else:
                        print_status_item(
                            f"Found segment dirs from processed ZIP: {', '.join(segment_dirs)}.",
                            "info",
                        )
                else:  # Should not happen if marker exists and flow worked
                    print_status_item(
                        f"Directory '{data_ref_path}' missing despite unzipped marker.",
                        "error",
                    )
                    overall_status_ok = False

        else:  # Named reference from references.tsv
            print_status_item(
                f"Reference mode: Named ('{reference_config_val}' from 'references.tsv').",
                "info",
            )
            try:
                ref_details = load_reference_dictionary(reference_config_val)
                if not ref_details:
                    print_status_item(
                        f"Named reference '{reference_config_val}' not found in 'references.tsv'.",
                        "error",
                    )
                    overall_status_ok = False
                else:
                    print_status_item(
                        f"Named ref '{reference_config_val}' configured. Segments: {', '.join(sorted(ref_details.keys()))}.",
                        "success",
                    )
            except FileNotFoundError:
                print_status_item("'references.tsv' not found.", "error")
                overall_status_ok = False
            except Exception as e:
                print_status_item(f"Error loading from 'references.tsv': {e}", "error")
                overall_status_ok = False
    elif not config:
        print_status_item("Skipping reference check: `config.yml` not loaded.", "info")
    elif "reference" not in config:
        print_status_item("`reference` key missing in `config.yml`.", "error")
        overall_status_ok = False

    # --- 3. Metadata File ---
    print_section_header(f"3. Metadata File ('{analysis_dir}/metadata.tsv')" if analysis_dir else "3. Metadata File")
    if not analysis_dir:
        print_status_item("Cannot check metadata: analysis directory not configured.", "error")
        overall_status_ok = False
    else:
        metadata_path = Path(analysis_dir) / "metadata.tsv"
        if not metadata_path.exists():
            print_status_item(f"'{analysis_dir}/metadata.tsv' not found.", "error")
            print_guidance(
                "This file is generated by: `python mlip/dataflow.py preprocess -f <your_id_list.txt>`."
            )
            print("")
            print_guidance(
                "Alternatively, you may want to use an existing metadata spreadsheet."
            )
            print_guidance(f"This would need to be manually placed at '{analysis_dir}/metadata.tsv'.")
            overall_status_ok = False  # Treat missing metadata as blocking the next step
        else:
            print_status_item(f"'{analysis_dir}/metadata.tsv' found.", "success")
        try:
            df = pd.read_csv(metadata_path, sep="\t", dtype=str).fillna(
                ""
            )  # Read all as str, fill NaN with empty str
                if df.empty and not list(
                    df.columns
                ):  # No columns means truly empty file or just whitespace
                    print_status_item(
                        f"'{analysis_dir}/metadata.tsv' is an empty file (no headers, no data).",
                        "warning",
                    )
                    print_guidance(
                        "Run `preprocess` if you haven't, or ensure it's correctly formatted."
                    )
                elif (
                    df.iloc[:, 0].eq("").all() and len(df.columns) <= 1
                ):  # Heuristic for empty content with just headers
                    print_status_item(
                        f"'{analysis_dir}/metadata.tsv' appears to have headers but no data rows.",
                        "warning",
                    )
                    print_guidance("Please populate the file with your sample information.")
                else:
                    required_cols = ["SequencingId", "SampleId", "Replicate"]
                    missing_cols = [col for col in required_cols if col not in df.columns]
                    if missing_cols:
                        print_status_item(
                            f"'{analysis_dir}/metadata.tsv' is missing required columns: {', '.join(missing_cols)}.",
                            "error",
                        )
                        overall_status_ok = False
                    else:
                        print_status_item(
                            "Required columns found in metadata file.", "success"
                        )
                        # Check for non-empty SampleId and Replicate columns in data rows
                        if df.empty:  # Only headers, no data rows
                            print_status_item(
                                f"'{analysis_dir}/metadata.tsv' has headers but no data rows.",
                                "warning",
                            )
                            print_guidance(
                                "Please populate it with your sample information."
                            )
                            overall_status_ok = False
                        elif df["SampleId"].eq("").any() or df["Replicate"].eq("").any():
                            print_status_item(
                                f"'{analysis_dir}/metadata.tsv' has some empty values in 'SampleId' or 'Replicate'.",
                                "error",
                            )
                            print_guidance(
                                "Ensure these are fully populated for all intended samples."
                            )
                            overall_status_ok = False
                        else:
                            print_status_item(
                                "Metadata file appears populated.", "success"
                            )
            except pd.errors.EmptyDataError:
                print_status_item(
                    f"'{analysis_dir}/metadata.tsv' exists but is completely empty (cannot be parsed).",
                    "error",
                )
                print_guidance(
                    "Run `python mlip/dataflow.py preprocess -h` and follow instructions therein."
                )
                overall_status_ok = False
            except Exception as e:
                print_status_item(
                    f"Error reading or parsing '{analysis_dir}/metadata.tsv': {e}", "error"
                )
                overall_status_ok = False

    return overall_status_ok


def load_metadata_dictionary(analysis_dir):
    f = open(f"{analysis_dir}/metadata.tsv", "r")
    reader = csv.DictReader(f, delimiter="\t")
    md_dict = defaultdict(lambda: defaultdict(list))
    counter = Counter()
    for row in reader:
        sample_id = row["SampleId"]
        replicate = row["Replicate"]
        counter[sample_id] += 1
        md_dict[sample_id][replicate].append(counter[sample_id])
    f.close()
    return md_dict


def samples_to_analyze(analysis_dir):
    manifest_filepath = Path(analysis_dir) / "file_manifest.json"
    samples_with_valid_data = set()

    with open(manifest_filepath, "r") as f:
        manifest_data = json.load(f)

    for sample_id, replicates_data in manifest_data.get("samples", {}).items():
        for _, experiments_list in replicates_data.items():
            for exp in experiments_list:
                if not exp.get("is_low_coverage", True):
                    samples_with_valid_data.add(sample_id)
                    break
            if sample_id in samples_with_valid_data:
                break
    return sorted(list(samples_with_valid_data))


def genbank_to_gtf(gbk_file, gtf_file, segment):
    with open(gtf_file, "w") as gtf_out:
        for record in SeqIO.parse(gbk_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    # Extract required fields
                    chrom = segment
                    source = "GenBank"
                    feature_type = "CDS"
                    start = int(feature.location.start) + 1  # GTF uses 1-based indexing
                    end = int(feature.location.end)  # Already 1-based in GenBank
                    score = "."
                    strand = "+" if feature.location.strand == 1 else "-"
                    frame = str(feature.qualifiers.get("codon_start", ["0"])[0])
                    gene_name = feature.qualifiers.get("gene", ["unknown"])[0]
                    transcript_id = feature.qualifiers.get("protein_id", ["unknown"])[0]

                    # Format the GTF line
                    gtf_line = f'{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tgene_id "{gene_name}"; transcript_id "{transcript_id}";\n'
                    gtf_out.write(gtf_line)


def clean_varscan(varscan_df):
    varscan_df["Frequency"] = [
        float(format_.split(":")[6][:-1]) / 100 for format_ in varscan_df["Sample1"]
    ]
    return varscan_df.loc[:, ["#CHROM", "POS", "REF", "ALT", "Frequency"]]


def merge_varscan(clean_varscan_dfs):
    variant_dictionary = {}
    for idx, clean_varscan_df in enumerate(clean_varscan_dfs):
        replicate_key = "replicate_%d" % (idx + 1)
        for _, row in clean_varscan_df.iterrows():
            if row["gene"] in ("PA-X", "PB1-F2"):
                continue
            variant_key = (
                row["segment"],
                row["reference_position"],
                row["reference_allele"],
            )
            variant_data = {
                "Frequency": row["frequency"],
                "CodingRegionChange": row["coding_region_change"],
                "Gene": row["gene"],
            }
            if variant_key in variant_dictionary:
                variant_dictionary[variant_key][replicate_key] = variant_data
            else:
                variant_dictionary[variant_key] = {replicate_key: variant_data}

    variant_list = []
    for variant_key, variant_value in variant_dictionary.items():
        coding_region_change = None
        gene = None
        for idx in range(len(clean_varscan_dfs)):
            replicate_key = "replicate_%d" % (idx + 1)
            if not replicate_key in variant_value:
                variant_value[replicate_key] = {"Frequency": 0}

        flattened = {}
        for idx in range(len(clean_varscan_dfs)):
            replicate_key = "replicate_%d" % (idx + 1)
            replicate_attributes = variant_value[replicate_key]
            for attribute_key, attribute_value in replicate_attributes.items():
                if attribute_key == "CodingRegionChange":
                    coding_region_change = attribute_value
                elif attribute_key == "Gene":
                    gene = attribute_value
                else:
                    flattened[f"{attribute_key}_{idx+1}"] = attribute_value
        variant = {
            "segment": variant_key[0],
            "position": variant_key[1],
            "allele": variant_key[2],
            "coding_region_change": coding_region_change,
            "gene": gene,
            **flattened,
        }
        variant_list.append(variant)

    if not variant_list:
        # If no variants were found, create an empty DataFrame with the expected columns.
        num_reps = len(clean_varscan_dfs)
        freq_cols = [f"Frequency_{i+1}" for i in range(num_reps)]
        final_cols = [
            "segment",
            "position",
            "allele",
            "coding_region_change",
            "gene",
        ] + freq_cols
        return pd.DataFrame(columns=final_cols)

    return pd.DataFrame(variant_list).sort_values(
        by=["segment", "position"], ascending=[True, True]
    )


def merge_varscan_io(input_tsv_filepaths, output_tsv_filepath):
    dfs = [
        pd.read_csv(input_tsv_filepath, sep="\t", na_filter=False)
        for input_tsv_filepath in input_tsv_filepaths
    ]
    merged_df = merge_varscan(dfs)
    merged_df.to_csv(output_tsv_filepath, sep="\t", index=False)


def get_coverage_bucket_labels(config):
    consensus_coverage = config["consensus_minimum_coverage"]
    variant_coverage = config["variants_minimum_coverage"]
    coverage_bucket_labels = [
        "0x",
        f"1-{consensus_coverage}x",
        f"{consensus_coverage}-{variant_coverage}x",
        f"{variant_coverage}-1000x",
        "1000x+",
    ]
    return coverage_bucket_labels


def assign_coverage_bucket(coverage):
    config = load_mlip_config()
    coverage_bucket_labels = get_coverage_bucket_labels(config)
    consensus_coverage = config["consensus_minimum_coverage"]
    variant_coverage = config["variants_minimum_coverage"]
    if coverage == 0:
        return coverage_bucket_labels[0]
    elif 1 <= coverage <= consensus_coverage:
        return coverage_bucket_labels[1]
    elif consensus_coverage <= coverage <= variant_coverage:
        return coverage_bucket_labels[2]
    elif variant_coverage <= coverage <= 1000:
        return coverage_bucket_labels[3]
    else:
        return coverage_bucket_labels[4]


def compute_coverage_categories(df):
    config = load_mlip_config()
    coverage_bucket_labels = get_coverage_bucket_labels(config)

    df["length"] = df["end"] - df["start"]

    def calculate_stats(group):
        total_length = group["length"].sum()
        if total_length == 0:
            return pd.Series({"mean_coverage": 0, "median_coverage": 0})

        weighted_sum = (group["coverage"] * group["length"]).sum()
        mean_cov = weighted_sum / total_length

        # Efficiently calculate median by repeating coverage value by its region length
        expanded_coverage = np.repeat(group["coverage"], group["length"])
        median_cov = np.median(expanded_coverage)

        return pd.Series({"mean_coverage": mean_cov, "median_coverage": median_cov})

    if not df.empty:
        stats_df = df.groupby("segment").apply(calculate_stats).reset_index()
    else:
        stats_df = pd.DataFrame(columns=["segment", "mean_coverage", "median_coverage"])

    if not df.empty:
        df["coverage_bucket"] = df["coverage"].apply(assign_coverage_bucket)
        bucket_summary = (
            df.groupby(["segment", "coverage_bucket"])["length"]
            .sum()
            .reset_index(name="count")
        )
    else:
        bucket_summary = pd.DataFrame(columns=["segment", "coverage_bucket", "count"])

    all_buckets = pd.DataFrame({"coverage_bucket": coverage_bucket_labels})
    segments = df["segment"].unique() if not df.empty else []

    if len(segments) > 0:
        all_combinations = pd.MultiIndex.from_product(
            [segments, all_buckets["coverage_bucket"]],
            names=["segment", "coverage_bucket"],
        ).to_frame(index=False)
        bucket_summary = pd.merge(
            all_combinations,
            bucket_summary,
            on=["segment", "coverage_bucket"],
            how="left",
        ).fillna(0)
        desired_structure = bucket_summary.pivot_table(
            index="segment", columns="coverage_bucket", values="count", fill_value=0
        ).reset_index()[["segment"] + coverage_bucket_labels]
        desired_structure.columns.name = None
    else:
        desired_structure = pd.DataFrame(columns=["segment"] + coverage_bucket_labels)

    final_summary_df = pd.merge(
        desired_structure, stats_df, on="segment", how="outer"
    ).fillna(0)

    final_cols = [
        "segment",
        "mean_coverage",
        "median_coverage",
    ] + coverage_bucket_labels
    final_summary_df = final_summary_df.reindex(columns=final_cols, fill_value=0)

    return final_summary_df


def compute_coverage_categories_io(input_coverage, output_summary):
    coverage_df = pd.read_csv(input_coverage, sep="\t")
    coverage_summary_df = compute_coverage_categories(coverage_df)
    coverage_summary_df.to_csv(output_summary, sep="\t", index=False)


def extract_cds_with_coords(genbank_file, consensus_record):
    """
    Extracts CDS data using coordinates from a GenBank file but sequence
    from a given consensus sequence record.
    """
    cds_data = []
    gb_record = SeqIO.read(genbank_file, "genbank")

    for feature in gb_record.features:
        if feature.type != "CDS":
            continue

        coding_seq = str(feature.location.extract(consensus_record.seq))

        # Rebuild the coordinate list as before
        coord_parts = []
        parts = feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]
        for part in parts:
            start, end = int(part.start), int(part.end)
            if part.strand == 1:
                coord_parts.extend(range(start, end))
            else:
                coord_parts.extend(reversed(range(start, end)))

        gene_name = feature.qualifiers.get("gene", [gb_record.id])[0]
        cds_data.append(
            {"gene": gene_name, "coding": coding_seq, "positions": coord_parts}
        )

    return cds_data


def extract_coding_regions(segments, replicate_consensus_fasta):
    """
    Loads consensus sequences and orchestrates CDS extraction for each segment.
    """
    # Load all consensus sequences into a dictionary for quick access
    consensus_sequences = SeqIO.to_dict(SeqIO.parse(replicate_consensus_fasta, "fasta"))

    coding_regions = {}
    for segment in segments:
        # Proceed only if the segment exists in the consensus file
        if segment in consensus_sequences:
            genbank_filepath = f"data/reference/{segment}/metadata.gb"
            # Pass the specific consensus sequence record to the extraction function
            consensus_record = consensus_sequences[segment]
            coding_regions[segment] = extract_cds_with_coords(genbank_filepath, consensus_record)

    return coding_regions


def extract_coding_regions_io(segments, replicate_consensus_fasta, output_json):
    """
    I/O wrapper for extracting coding regions for a replicate.
    """
    with open(output_json, "w") as json_file:
        coding_regions = extract_coding_regions(segments, replicate_consensus_fasta)
        json.dump(coding_regions, json_file)


def define_coding_regions(gtf_file):
    with open(gtf_file, "r") as gtf:
        coding_regions = {}

        for line in gtf:
            if (
                line.strip("\n") != ""
            ):  # ignore blank lines (otherwise throws an index error)
                sequence_name = line.split("\t")[0]
                annotation_type = line.split("\t")[2]
                start = (
                    int(line.split("\t")[3]) - 1
                )  # adding the -1 here for 0 indexing
                stop = int(line.split("\t")[4]) - 1  # adding the -1 here for 0 indexing
                gene_name = line.split("\t")[8]
                gene_name = gene_name.split(";")[0]
                gene_name = gene_name.replace("gene_id ", "")
                gene_name = gene_name.replace('"', "")

                if annotation_type.lower() == "cds":
                    if sequence_name not in coding_regions:
                        coding_regions[sequence_name] = {}
                        coding_regions[sequence_name][gene_name] = [start, stop]
                    elif (
                        sequence_name in coding_regions
                        and gene_name not in coding_regions[sequence_name]
                    ):
                        coding_regions[sequence_name][gene_name] = [start, stop]
                    elif gene_name in coding_regions[sequence_name]:
                        coding_regions[sequence_name][gene_name].extend([start, stop])

            # sort coding region coordinates so that they are always in the correct order
            for sequence_name in coding_regions:
                for gene in coding_regions[sequence_name]:
                    coding_regions[sequence_name][gene] = sorted(
                        coding_regions[sequence_name][gene]
                    )

    return coding_regions


amino_acid_abbreviations = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "O": "Pyl",
    "S": "Ser",
    "U": "Sec",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "B": "Asx",
    "Z": "Glx",
    "U": "Sec",
    "X": "Xaa",
    "J": "Xle",
    "*": "Stop",
}


def lookup_site(arr, site):
    # return the location of a site in an array, or None if it is not present
    return next((i for i, val in enumerate(arr) if val == site), None)


def slice_fastas(coding_regions, reference_sequence_path):
    reference_sequence = {}

    for seq in SeqIO.parse(reference_sequence_path, "fasta"):
        sequence_name = str(seq.id)
        sequence = str(seq.seq).lower()
        reference_sequence[sequence_name] = sequence

    transcripts = {}
    stop_codons = ["taa", "tag", "tga"]
    refseqname = "A sequence"

    for c in coding_regions:
        for gene in coding_regions[c]:
            transcripts[gene] = ""
            coordinates = coding_regions[c][
                gene
            ]  # define the coding regions for each gene

            for i in range(
                0, int(len(coordinates)), 2
            ):  # go along the coordinates in chunks of 2 at a time
                sequence_chunk = reference_sequence[c][
                    coordinates[i] : coordinates[i + 1] + 1
                ]
                transcripts[
                    gene
                ] += sequence_chunk  # append each piece of the transcript together

    # loop through each transcript to make sure that it begins with a start codon and ends with a stop codon
    for t in transcripts:
        if transcripts[t][0:3] != "atg":
            print(
                "WARNING! " + refseqname + " " + t + "does not contain a start codon!"
            )
        if transcripts[t][-3:] not in stop_codons:
            print(
                "WARNING! "
                + refseqname
                + " "
                + t
                + " does not contain a stop codon! These are the last 3 nucleotides: "
                + transcripts[t][-3:]
            )

    return transcripts


def annotate_amino_acid_changes(coding_regions, vcf, outfilename):
    csvfile = open(vcf, "r")
    outfile = open(outfilename, "w")

    reader = csv.reader(csvfile, delimiter="\t")
    writer = csv.writer(outfile, delimiter="\t")

    header = [
        "segment",
        "gene",
        "reference_position",
        "amino_acid_position",
        "reference_allele",
        "variant_allele",
        "reference_codon",
        "variant_codon",
        "reference_amino_acid",
        "variant_amino_acid",
        "coding_region_change",
        "synonymous/nonsynonymous",
        "frequency(%)",
        "frequency",
    ]
    writer.writerow(header)

    for row in reader:
        # ignore comment lines
        if "##" in row[0] or "#CHROM" in row[0]:
            continue

        segment = row[0]
        site = int(row[1]) - 1  # to make this 0 indexed
        reference_allele = row[3].lower()
        alternative_allele = row[4].lower()

        # pull out the frequency using a string search
        SearchStr = r".+\:([0-9]{1,2}\.{0,1}[0-9]{0,2}\%)"
        result = re.search(SearchStr, row[9])
        if result:
            frequency = "\t".join(result.groups())
            frequency_decimal = frequency.replace("%", "")
            frequency_decimal = (float(frequency_decimal)) / 100
            frequency_decimal = (
                "%.4f" % frequency_decimal
            )  # only include 4 numbers after decimal
        else:
            frequency = "none reported"

        for gene in coding_regions[segment]:
            gene_name = gene["gene"]
            cds = gene["coding"]
            positions = gene["positions"]
            position = lookup_site(positions, site)
            if position is None:
                continue
            # how far back should we go for the start of the codon?
            # zero-based indexing and modular arithmetic pair well here
            offset = position % 3
            start = position - offset
            stop = start + 3
            ref_codon = cds[start:stop]
            ref_aa = Seq(ref_codon).translate()

            # we convert to a list to dynamically modify a base, then join to a string
            variant_codon_list = list(ref_codon)
            variant_codon_list[offset] = alternative_allele
            variant_codon = "".join(variant_codon_list)
            variant_aa = Seq(variant_codon).translate()

            # return the amino acid changes, with single letters converted to 3-letter aa codes
            if ref_aa == variant_aa:
                syn_nonsyn = "synonymous"
            elif ref_aa != variant_aa and variant_aa == "*":
                syn_nonsyn = "stop_gained"
            elif ref_aa != variant_aa:
                syn_nonsyn = "nonsynonymous"

            # now determine whether the site is in the 1st, 2nd, or 3rd codon position
            # if SNP is in 1st position in codon:
            aa_site = int(position / 3) + 1

            amino_acid_change = (
                amino_acid_abbreviations[ref_aa]
                + str(aa_site)
                + amino_acid_abbreviations[variant_aa]
            )
            output = [
                segment,
                gene_name,
                str(site + 1),
                str(aa_site),
                reference_allele.upper(),
                alternative_allele.upper(),
                ref_codon.upper(),
                variant_codon.upper(),
                str(ref_aa),
                str(variant_aa),
                amino_acid_change,
                syn_nonsyn,
                frequency,
                frequency_decimal,
            ]
            writer.writerow(output)

    csvfile.close()
    outfile.close()
    return


def coverage_summary(input_tsvs, output_tsv):
    config = load_mlip_config()
    coverage_bucket_labels = get_coverage_bucket_labels(config)
    dfs = []
    for tsv_path in input_tsvs:
        split_path = tsv_path.split("/")
        sample_id = split_path[1]
        replicate = split_path[2]
        df = pd.read_csv(tsv_path, sep="\t")
        df["sample_id"] = f"{sample_id}-{replicate}"
        dfs.append(df)
    full_df = pd.concat(dfs, ignore_index=True)
    full_df["total_coverage"] = (
        full_df[coverage_bucket_labels[0]]
        + full_df[coverage_bucket_labels[1]]
        + full_df[coverage_bucket_labels[2]]
        + full_df[coverage_bucket_labels[3]]
        + full_df[coverage_bucket_labels[4]]
    )

    full_df.to_csv(output_tsv, sep="\t", index=False)


def get_bases_and_qualities(raw_bases, qual_str):
    """Parses a pileup raw base string and its quality string, returning paired lists."""
    bases, qualities = [], []
    i, q_index = 0, 0
    while i < len(raw_bases):
        c = raw_bases[i]
        if c == "^":
            i += 2  # skip the read-start marker and its mapping quality
            continue
        if c == "$":
            i += 1
            continue
        if c in "+-":
            i += 1
            num = ""
            while i < len(raw_bases) and raw_bases[i].isdigit():
                num += raw_bases[i]
                i += 1
            if num:
                i += int(num)
            continue
        bases.append(c)
        if q_index < len(qual_str):
            qualities.append(ord(qual_str[q_index]) - 33)
            q_index += 1
        else:
            qualities.append(0)
        i += 1
    return bases, qualities


def parse_pileup(pileup_file, qual_threshold):
    """
    Parses a samtools pileup file and, for each position,
    computes the majority base and effective coverage using only bases with quality >= qual_threshold.

    Returns:
        dict: {segment: {position: (majority_base, effective_coverage)}}
    """
    pileup_data = {}
    with open(pileup_file) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 5:
                continue
            segment = fields[0]
            pos = int(fields[1])
            ref_base = fields[2].upper()
            raw_bases = fields[4]
            qual_str = fields[5] if len(fields) >= 6 else ""

            bases, quals = get_bases_and_qualities(raw_bases, qual_str)
            filtered = []
            for b, q in zip(bases, quals):
                if q >= qual_threshold:
                    filtered.append(ref_base if b in [".", ","] else b.upper())
            effective_cov = len(filtered)
            majority = Counter(filtered).most_common(1)[0][0] if filtered else "N"

            if segment not in pileup_data:
                pileup_data[segment] = {}
            pileup_data[segment][pos] = (majority, effective_cov)
    return pileup_data


def check_consensus(fasta_file, pileup_data, coverage_threshold):
    """
    Compares a consensus FASTA with pileup-derived majority bases (from high-quality bases only).
    Flags positions where:
      - A base is present when effective coverage is below threshold (should be 'N').
      - An 'N' is present despite sufficient effective coverage.
      - The consensus base does not match the majority call.
    """
    errors = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        segment = record.id
        seq = record.seq.upper()
        if segment not in pileup_data:
            print(f"Warning: Segment {segment} not found in pileup.")
            continue
        for pos in range(1, len(seq) + 1):
            base = seq[pos - 1]
            majority, cov = pileup_data[segment].get(pos, ("N", 0))
            expected = "N" if cov < coverage_threshold else majority
            if base != expected:
                if base == "N" and cov >= coverage_threshold:
                    err_type = "Unexpected N"
                elif base != "N" and cov < coverage_threshold:
                    err_type = "N expected but not present"
                else:
                    err_type = "Consensus base mismatch"
                errors.append((segment, pos, base, expected, cov, err_type))
    if errors:
        print("Discrepancies found:")
        print(
            f"{'Segment':<12}{'Pos':<8}{'Consensus':<12}{'Expected':<12}{'Coverage':<10}{'Issue'}"
        )
        for seg, pos, base, exp, cov, issue in errors:
            print(f"{seg:<12}{pos:<8}{base:<12}{exp:<12}{cov:<10}{issue}")
    else:
        print("No discrepancies found.")
    return errors


def check_consensus_io(input_consensus, input_pileup, output_tsv, sample, replicate):
    config = load_mlip_config()
    coverage_threshold = config["consensus_minimum_coverage"]
    quality_threshold = config["minimum_quality_score"]
    pileup_data = parse_pileup(input_pileup, 0)
    output = check_consensus(input_consensus, pileup_data, coverage_threshold)
    headers = [
        "Sample",
        "Replicate",
        "Segment",
        "Position",
        "Consensus",
        "Expected",
        "Coverage",
        "Issue",
    ]
    with open(output_tsv, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        writer.writerows([(sample, replicate) + row for row in output])


def merge_variant_calls(input, output):
    if len(input) == 0:
        Path(output).touch()
        return
    dfs = []
    for fp in input:
        parts = fp.split("/")
        sample = parts[1]
        df = pd.read_csv(fp, sep="\t", na_filter=False)
        df["sample"] = sample
        dfs.append(df)

    # Concatenate all DataFrames and write to CSV
    merged = pd.concat(dfs, ignore_index=True)
    merged.to_csv(output, index=False, sep="\t")


unambig = set("ATCG")


def call_sample_consensus(input_replicates, output_sample):

    if len(input_replicates) == 1:
        shutil.copy(input_replicates[0], output_sample)
        return

    # Use Biopython's built-in function to create a dictionary (assumes unique IDs)
    dicts = [SeqIO.to_dict(SeqIO.parse(f, "fasta")) for f in input_replicates]
    output_records = []

    for rec_id in dicts[0]:
        if rec_id not in dicts[1]:
            continue

        r1 = dicts[0][rec_id]
        r2 = dicts[1][rec_id]
        s1, s2 = str(r1.seq), str(r2.seq)

        if not s1 or not s2:
            if not s1 and s2:
                print(
                    f"Warning: record {rec_id} is empty in first replicate. Using second replicate."
                )
                output_records.append(r2)
            elif not s2 and s1:
                print(
                    f"Warning: record {rec_id} is empty in second replicate. Using first replicate."
                )
                output_records.append(r1)
            else:
                print(f"Warning: record {rec_id} is empty in both replicates.")
                output_records.append(r1)
            continue

        merged = []
        for pos, (b1, b2) in enumerate(zip(s1, s2)):
            if b1 == b2:
                merged.append(b1)
            elif b1 not in unambig and b2 in unambig:
                merged.append(b2)
            elif b2 not in unambig and b1 in unambig:
                merged.append(b1)
            else:
                merged.append("N")
                print(
                    f"Warning: conflict in record {rec_id} at position {pos+1}: {b1} vs {b2}"
                )

        consensus_record = SeqRecord(Seq("".join(merged)), id=rec_id, description="")
        output_records.append(consensus_record)

    SeqIO.write(output_records, output_sample, "fasta")


def translate_consensus_genes(consensus_fasta, output_dir, sample):
    """
    For each consensus record (whose id corresponds to a segment/genbank file),
    translate each CDS gene using the consensus sequence (aligned to the GenBank).
    The protein sequences are written to output_dir/{gene}.fasta.

    Parameters:
      consensus_fasta : str
          Path to consensus FASTA file.
      genbank_dir : str
          Directory containing GenBank files named by segment id (e.g., SEGID.gb).
      output_dir : str
          Directory where the translation files will be written.
      seg_start : int, optional
          Offset of the segment relative to the GenBank sequence (default 0).
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load consensus records into a dict keyed by id.
    consensus_dict = SeqIO.to_dict(SeqIO.parse(consensus_fasta, "fasta"))

    for seg_id, consensus_record in consensus_dict.items():
        gb_file = os.path.join("data", "reference", seg_id, "metadata.gb")
        try:
            gb_record = SeqIO.read(gb_file, "genbank")
        except Exception as e:
            print(f"Error reading GenBank file {gb_file}: {e}")
            continue

        cons_seq = consensus_record.seq

        # Process each CDS feature in the GenBank record.
        for feature in gb_record.features:
            if feature.type != "CDS":
                continue
            gene_start = int(feature.location.start)
            gene_end = int(feature.location.end)

            # Map gene coordinates to consensus.
            gene_seq = cons_seq[gene_start:gene_end]

            # Get codon table; default to standard table 1.
            table = int(feature.qualifiers.get("transl_table", ["1"])[0])
            protein = gene_seq.translate(table=table, to_stop=True)

            gene_id = feature.qualifiers.get("gene", [f"{gene_start}_{gene_end}"])[0]
            output_file = os.path.join(output_dir, f"{gene_id}.fasta")

            # Write the translation. The record id is set to the segment id.
            protein_record = SeqRecord(protein, id=sample, description=gene_id)
            SeqIO.write(protein_record, output_file, "fasta")


def fill(input_alignment, output_sequence):
    fasta = SeqIO.parse(input_alignment, "fasta")
    background = next(fasta)
    foreground = next(fasta)
    output_bases = []
    for background_base, foreground_base in zip(background.seq, foreground.seq):
        if background_base != "-":
            if not foreground_base in ["-", "N"]:
                output_bases.append(foreground_base)
            else:
                output_bases.append(background_base)
    record = SeqRecord(
        Seq("".join(output_bases)), id=foreground.id, description="hybrid"
    )
    SeqIO.write(record, output_sequence, "fasta")


def get_duplicate_samples(metadata_dictionary):
    duplicate_samples = [
        sample for sample, metadata in metadata_dictionary.items() if len(metadata) == 2
    ]
    return duplicate_samples


def compare_remappings_io(
    penultimate_fasta_path, final_fasta_path, output_tsv_path, sample, replicate
):
    """
    Compare consensus sequences between penultimate and final remapping stages.

    Parameters:
        penultimate_fasta_path: Path to the penultimate remapping consensus FASTA
        final_fasta_path: Path to the final remapping consensus FASTA
        output_tsv_path: Path to write the consensus summary report
        sample: Sample ID for the comparison
        replicate: Replicate number for the comparison
    """
    try:
        penultimate_dict = SeqIO.to_dict(SeqIO.parse(penultimate_fasta_path, "fasta"))
        final_dict = SeqIO.to_dict(SeqIO.parse(final_fasta_path, "fasta"))

        differences = []

        for segment_id in penultimate_dict:
            if segment_id not in final_dict:
                continue

            penultimate_seq = str(penultimate_dict[segment_id].seq)
            final_seq = str(final_dict[segment_id].seq)

            if penultimate_seq == final_seq:
                continue

            min_length = min(len(penultimate_seq), len(final_seq))
            for pos in range(min_length):
                if penultimate_seq[pos] != final_seq[pos]:
                    differences.append(
                        {
                            "sample": sample,
                            "replicate": replicate,
                            "segment": segment_id,
                            "position": pos + 1,  # 1-based positioning
                            "penultimate_base": penultimate_seq[pos],
                            "final_base": final_seq[pos],
                        }
                    )

        if differences:
            df = pd.DataFrame(differences)
            df.to_csv(output_tsv_path, sep="\t", index=False)
        else:
            # Create empty file with header if no differences found
            df = pd.DataFrame(
                columns=[
                    "sample",
                    "replicate",
                    "segment",
                    "position",
                    "penultimate_base",
                    "final_base",
                ]
            )
            df.to_csv(output_tsv_path, sep="\t", index=False)

    except Exception as e:
        print(f"Error comparing remappings for {sample} replicate {replicate}: {e}")
        Path(output_tsv_path).touch()


def aggregate_consensus_summaries_io(input_files, output_file):
    """
    Aggregate consensus summary TSV files, handling empty files gracefully.

    Parameters:
        input_files: List of input TSV file paths
        output_file: Path to write aggregated output
    """
    # Expected columns for empty case
    columns = [
        "sample",
        "replicate",
        "segment",
        "position",
        "penultimate_base",
        "final_base",
    ]

    # Read all non-empty files
    dfs = []
    for file_path in input_files:
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            try:
                df = pd.read_csv(file_path, sep="\t")
                if not df.empty:
                    dfs.append(df)
            except pd.errors.EmptyDataError:
                continue

    # Combine or create empty DataFrame
    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
    else:
        combined_df = pd.DataFrame(columns=columns)

    # Write output
    combined_df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    command_line_interface()
