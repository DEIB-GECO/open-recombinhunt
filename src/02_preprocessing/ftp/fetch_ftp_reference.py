# src/02_preprocessing/ftp/fetch_ftp_reference.py

import argparse
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

def fetch_reference_step(virus_name: str, config: dict):
    """
    Orchestrates the fetching of the reference FASTA sequence from NCBI.
    """
    logging.info(f"--- Starting Step 2.3: Fetch Reference Sequence for '{virus_name}' ---")

    # --- 1. Get Configuration and Paths ---
    paths_config = config.get(PATHS, {})
    virus_config = config.get(VIRUSES, {}).get(virus_name, {})
    
    processed_dir = Path(paths_config.get(PROCESSED_DATA, 'data/processed')) / virus_name
    
    # --- 2. Get Reference Accession ID ---
    ref_accession = virus_config.get(REFERENCE, {}).get(ACCESSION_ID)
    if not ref_accession:
        logging.error(f"No reference accession_id found in config for '{virus_name}'. Cannot proceed.")
        sys.exit(1)

    logging.info(f"Found reference accession ID: {ref_accession}")

    # --- 3. Create Temporary Working Directory ---
    temp_work_dir = processed_dir / "temp_fetch_reference"
    temp_work_dir.mkdir(exist_ok=True)
    logging.info(f"Created temporary working directory: {temp_work_dir}")

    try:
        # --- 4. Get and Execute Download Command ---
        ftp_ref_commands = config.get(FTP_CLI, {}).get(REFERENCE, [])
        if not ftp_ref_commands:
            logging.error("No reference download commands found in config under 'ftp-cli'.")
            raise ValueError("Missing CLI command in config.")

        command_template = ftp_ref_commands[0]
        command = command_template.format(
            accession_id=ref_accession,
            virus_name=virus_name
        )
        
        run_command(command, working_dir=temp_work_dir)

        # --- 5. Unzip and Finalize Output ---
        zip_filename = f"{virus_name}-reference.zip"
        zip_filepath = temp_work_dir / zip_filename
        
        logging.info(f"Unzipping {zip_filepath}...")
        with zipfile.ZipFile(zip_filepath, 'r') as zip_ref:
            zip_ref.extractall(temp_work_dir)
            
        # Find the .fna file, rename it to reference.fasta, and move it to the final destination
        unzipped_data_path = temp_work_dir / "ncbi_dataset" / "data"
        find_and_rename(unzipped_data_path, "*.fna", "reference.fasta", processed_dir)

    finally:
        # --- 6. Cleanup ---
        logging.info("Cleaning up temporary working directory...")
        if temp_work_dir.exists():
            shutil.rmtree(temp_work_dir)
        logging.info("Cleanup complete.")


def main():
    """Main function to run the reference sequence fetching step."""
    parser = argparse.ArgumentParser(description="Fetch the reference FASTA sequence for a virus.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process (e.g., 'influenza').")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    try:
        config = yaml.safe_load(Path(args.config).read_text())
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_02.3_fetch_reference")

    try:
        fetch_reference_step(args.virus, config)
        logging.info("Reference sequence fetching step finished successfully.")
    except Exception as e:
        logging.critical(f"Reference sequence fetching step failed with a critical error: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
