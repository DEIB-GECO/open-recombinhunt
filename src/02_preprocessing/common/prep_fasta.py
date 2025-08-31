# src/02_preprocessing/ncbi/prep_fasta_ncbi.py

import argparse
from pathlib import Path
import sys
import logging
import yaml
from typing import Set

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
    from Bio import SeqIO
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    print("Please ensure Biopython is installed ('pip install biopython').")
    sys.exit(1)

def load_ids_from_file(filepath: Path) -> Set[str]:
    """Loads sequence IDs from a text file into a set, ignoring empty lines."""
    logging.info(f"Loading IDs from: {filepath}")
    try:
        with open(filepath, 'r') as f:
            # Use a set comprehension for efficiency
            ids_set = {line.strip() for line in f if line.strip()}
        logging.info(f"Loaded {len(ids_set)} unique IDs.")
        return ids_set
    except FileNotFoundError:
        logging.error(f"Required ID file not found at {filepath}")
        # We can continue with an empty set, but the filtering will be incomplete.
        # It's better to exit if a required file is missing.
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading ID file {filepath}: {e}")
        sys.exit(1)

def filter_fasta_file(
    input_fasta_path: Path,
    output_fasta_path: Path,
    ids_to_drop_set: Set[str]
):
    """
    Filters an input FASTA file, writing sequences to the output path
    only if their ID is NOT in the ids_to_drop_set.
    """
    logging.info(f"Filtering FASTA file: {input_fasta_path}...")
    logging.info(f"Removing sequences whose ID is in the combined dropped list ({len(ids_to_drop_set)} total IDs).")
    
    kept_records = []
    original_count = 0
    kept_count = 0

    try:
        with open(input_fasta_path, "r") as handle_in:
            for record in SeqIO.parse(handle_in, "fasta"):
                original_count += 1
                # The record.id from Biopython typically gives the main identifier before any spaces.
                # For NCBI FASTA files, this is usually the accession number.
                sequence_id = record.id
                
                if sequence_id not in ids_to_drop_set:
                    kept_records.append(record)
                    kept_count += 1
        
        # Write all kept records in one go for efficiency
        with open(output_fasta_path, "w") as handle_out:
            SeqIO.write(kept_records, handle_out, "fasta")

    except FileNotFoundError:
        logging.error(f"Input FASTA file not found at {input_fasta_path}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error occurred during FASTA processing: {e}", exc_info=True)
        sys.exit(1)

    logging.info(f"FASTA Filtering Summary:")
    logging.info(f"  Original sequences read: {original_count}")
    logging.info(f"  Sequences dropped      : {original_count - kept_count}")
    logging.info(f"  Sequences kept         : {kept_count}")
    logging.info(f"Filtered FASTA saved to  : {output_fasta_path}")

def main():
    """Main function to run the FASTA preprocessing."""
    parser = argparse.ArgumentParser(description="Filter a raw FASTA file based on dropped ID lists.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    config = yaml.safe_load(Path(args.config).read_text())
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_02.2_prep_fasta")

    paths_config = config.get(PATHS, {})
    
    raw_dir = Path(paths_config.get(RAW_DATA, 'data/raw')) / args.virus
    processed_dir = Path(paths_config.get(PROCESSED_DATA, 'data/processed')) / args.virus
    
    # Define all file paths
    input_fasta = raw_dir / "raw_sequences.fasta"
    filtered_ids_file = processed_dir / "filtered_sequence_ids.txt"
    non_classified_ids_file = processed_dir / "non_classified_sequence_ids.txt"
    output_fasta = processed_dir / "sequences.fasta"

    # 1. Load both sets of IDs to be dropped
    logging.info("Loading IDs of sequences to be dropped...")
    filtered_ids = load_ids_from_file(filtered_ids_file)
    non_classified_ids = load_ids_from_file(non_classified_ids_file)
    logging.info(f"NUmber of Filtered IDs: {len(filtered_ids)}")
    logging.info(f"Number of Non-Classified IDs: {len(non_classified_ids)}")

    # 2. Combine them into a single master set
    master_drop_set = filtered_ids.union(non_classified_ids)
    logging.info(f"Total unique IDs to drop from FASTA file: {len(master_drop_set)}")

    # 3. Filter the FASTA file
    filter_fasta_file(input_fasta, output_fasta, master_drop_set)

    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()
