# src/02_preprocessing/ftp/fetch_ftp_sequences_by_id.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import yaml
import shutil
import zipfile

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging, run_command, find_and_rename
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    sys.exit(1)

def fetch_sequences_step(virus_name: str, config: dict):
    """
    Orchestrates the fetching of FASTA sequences based on a list of accession IDs.
    """
    logging.info(f"--- Starting Step 2.2: Fetch Sequences by ID for '{virus_name}' ---")

    # --- 1. Get Configuration and Paths ---
    paths_config = config.get(PATHS, {})
    
    processed_dir = Path(paths_config.get(PROCESSED_DATA, 'data/processed')) / virus_name
    
    accession_ids_filepath = processed_dir / f"{virus_name}-accession-ids.txt"
    output_fasta_filepath = processed_dir / "sequences.fasta"

    # --- 2. Check for Accession IDs File ---
    if not accession_ids_filepath.exists():
        logging.error(f"Accession ID file not found at: {accession_ids_filepath}")
        logging.error("Please run the 'prep_metadata_ftp.py' script first.")
        sys.exit(1)
        
    # Check if the file is empty
    if accession_ids_filepath.stat().st_size == 0:
        logging.warning("The accession IDs file is empty. No sequences will be downloaded.")
        # Create an empty sequences.fasta file for pipeline consistency
        output_fasta_filepath.touch()
        logging.info(f"Created empty sequences file at: {output_fasta_filepath}")
        return

    # --- 3. Create Temporary Working Directory ---
    temp_work_dir = processed_dir / "temp_fetch_sequences"
    temp_work_dir.mkdir(exist_ok=True)
    logging.info(f"Created temporary working directory: {temp_work_dir}")

    try:
        # --- 4. Prepare for CLI Command ---
        # The NCBI tool expects the input file to be in the current working directory.
        # So, we copy our accession ID list into the temp directory.
        temp_accession_ids_path = temp_work_dir / f"{virus_name}-accession-ids.txt"
        shutil.copy(accession_ids_filepath, temp_accession_ids_path)

        # Get the command template from the config
        sequence_commands = config.get(FTP_CLI, {}).get(SEQUENCES, {})
        if not sequence_commands:
            logging.error(f"No sequence download commands found for '{virus_name}' in config under 'ftp-cli'.")
            raise ValueError("Missing CLI command in config.")

        # --- 5. Execute Download Command ---
        command_template = sequence_commands[0]
        command = command_template.format(virus_name=virus_name)
        
        run_command(command, working_dir=temp_work_dir)

        # --- 6. Unzip and Finalize Output ---
        zip_filename = f"{virus_name}-sequences.zip"
        zip_filepath = temp_work_dir / zip_filename
        
        logging.info(f"Unzipping {zip_filepath}...")
        with zipfile.ZipFile(zip_filepath, 'r') as zip_ref:
            zip_ref.extractall(temp_work_dir)
            
        # Find the .fna file, rename it to sequences.fasta, and move it to the final destination
        find_and_rename(temp_work_dir, "*.fna", "sequences.fasta", processed_dir)

    finally:        
        # --- 7. Cleanup ---
        logging.info("Cleaning up temporary working directory...")
        if temp_work_dir.exists():
            shutil.rmtree(temp_work_dir)
        logging.info("Cleanup complete.")

def main():
    """Main function to run the sequence fetching step."""
    parser = argparse.ArgumentParser(description="Fetch FASTA sequences for a given list of accession IDs.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process (e.g., 'influenza').")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    try:
        config = yaml.safe_load(Path(args.config).read_text())
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_02.2_fetch_sequences_by_id")

    try:
        fetch_sequences_step(args.virus, config)
        logging.info("Sequence fetching step finished successfully.")
    except Exception as e:
        logging.critical(f"Sequence fetching step failed with a critical error: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
