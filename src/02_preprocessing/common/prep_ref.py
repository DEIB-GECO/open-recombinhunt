# src/02_preprocessing/ncbi/prep_ref_ncbi.py

import argparse
from pathlib import Path
import sys
import logging
import yaml
import shutil

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    print("Please ensure 'utils.py' and 'constants.py' are in the 'src/utils' directory.")
    sys.exit(1)

def copy_reference_sequence(source_path: Path, destination_path: Path):
    """
    Copies the reference sequence from the raw data directory to the processed directory.
    """
    logging.info(f"Copying reference sequence...")
    logging.info(f"  Source: {source_path}")
    logging.info(f"  Destination: {destination_path}")
    
    try:
        shutil.copy(source_path, destination_path)
        logging.info("Reference sequence successfully copied to processed directory.")
    except FileNotFoundError:
        logging.error(f"Source reference file not found at '{source_path}'. Cannot copy.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error occurred while copying the reference file: {e}", exc_info=True)
        sys.exit(1)

def main():
    """Main function to run the reference sequence preparation."""
    parser = argparse.ArgumentParser(description="Copy the reference FASTA file from the raw to the processed directory.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    try:
        config = yaml.safe_load(Path(args.config).read_text())
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"CRITICAL ERROR: Error parsing YAML config file: {e}", file=sys.stderr)
        sys.exit(1)
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_02.3_prep_ref")

    paths_config = config.get(PATHS, {})
    
    raw_dir = Path(paths_config.get(RAW_DATA, 'data/raw')) / args.virus
    processed_dir = Path(paths_config.get(PROCESSED_DATA, 'data/processed')) / args.virus
    
    # Define file paths
    source_ref_file = raw_dir / "reference.fasta"
    dest_ref_file = processed_dir / "reference.fasta"

    # Run the copy function
    copy_reference_sequence(source_ref_file, dest_ref_file)

    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()
