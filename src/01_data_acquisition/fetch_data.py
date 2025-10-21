# src/01_data_acquisition/fetch_data.py

import argparse
import yaml
import requests
import zipfile
import lzma
import gzip
import shutil
from pathlib import Path
import sys
import logging
import zstandard as zstd
from datetime import datetime

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))

try:
    # Import from utils.py which is in the same 'src' directory.
    from utils.utils import setup_logging, run_command, find_and_rename
    # Import all constants from the constants file
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    print("Please ensure 'src/utils.py' and 'src/constants.py' exist and the script is run from the project root.")
    sys.exit(1)

# --- Main Fetching Logic Functions ---

def update_download_date(config_path: str, virus_name: str):
    """Update the download_date for a virus in the config.yaml file."""
    try:
        # Read the current config
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Update the download_date for the specific virus
        if VIRUSES in config and virus_name in config[VIRUSES]:
            current_date = datetime.now().strftime("%Y-%m-%d")
            config[VIRUSES][virus_name]["download_date"] = current_date
            
            # Write the updated config back to file
            with open(config_path, 'w') as f:
                yaml.dump(config, f, default_flow_style=False, sort_keys=False)
            
            logging.info(f"Updated download_date for {virus_name} to {current_date}")
        else:
            logging.warning(f"Could not find virus '{virus_name}' in config to update download_date")
            
    except Exception as e:
        logging.error(f"Failed to update download_date in config: {e}")

def download_and_decompress(url: str, final_filename: str, virus_raw_dir: Path):
        """Helper to download a file and decompress it based on its extension."""
        if not url:
            logging.warning(f"No URL provided for {final_filename}. Skipping.")
            return

        logging.info(f"Downloading from {url}")
        
        # Determine temporary filename and decompression method
        if url.endswith(".gz"):
            temp_path = virus_raw_dir / f"temp_{final_filename}.gz"
            decompressor = gzip.open
        elif url.endswith(".zst"):
            temp_path = virus_raw_dir / f"temp_{final_filename}.zst"
            decompressor = zstd.open
        elif url.endswith(".xz"):
            temp_path = virus_raw_dir / f"temp_{final_filename}.xz"
            decompressor = lzma.open
        else:
            logging.warning(f"Unrecognized compression format for URL: {url}. Assuming no compression.")
            temp_path = virus_raw_dir / final_filename
            decompressor = None

        final_path = virus_raw_dir / final_filename

        try:
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(temp_path, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)

            if decompressor:
                logging.info(f"Decompressing {temp_path}...")
                with decompressor(temp_path, 'rb') as f_in:
                    with open(final_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                temp_path.unlink() # Delete the compressed file
            
            logging.info(f"Successfully created {final_path}")

        except Exception as e:
            logging.error(f"Failed to download or process {url}: {e}")

def fetch_reference_sequence(virus_config: dict, global_config: dict, virus_raw_dir: Path):
    """
    Downloads and processes the reference sequence for a virus from NCBI.
    This is a self-contained helper function.
    """
    logging.info("Fetching reference sequence from NCBI...")
    ref_accession = virus_config.get(REFERENCE, {}).get(ACCESSION_ID)
    if not ref_accession:
        logging.warning("No reference accession_id found in config. Skipping reference download.")
        return

    # Create a dedicated temporary directory for the reference download
    temp_ref_dir = virus_raw_dir / "temp_ref_download"
    temp_ref_dir.mkdir(exist_ok=True)

    try:
        # Run download command in the dedicated 'ref' subdirectory
        ncbi_ref_commands = global_config.get(NCBI_CLI, {}).get(REFERENCE, [])
        for command_template in ncbi_ref_commands:
            command = command_template.format(
                accession_id=ref_accession,
                virus_name=virus_config[NAME]
            )
            run_command(command, temp_ref_dir)
            
        ref_zip_path = temp_ref_dir / f"{virus_config[NAME]}-reference.zip"
        logging.info(f"Unzipping {ref_zip_path}...")
        with zipfile.ZipFile(ref_zip_path, 'r') as zip_ref:
            zip_ref.extractall(temp_ref_dir)
            
        # Find and rename the reference file from the 'ref' subdirectory
        unzipped_data_path = temp_ref_dir / "ncbi_dataset" / "data"
        find_and_rename(unzipped_data_path, "genomic.fna", "reference.fasta", virus_raw_dir)

    finally:
        # Clean up the temporary reference directory
        logging.info("Cleaning up temporary reference sequence files...")
        if temp_ref_dir.exists():
            shutil.rmtree(temp_ref_dir)

def fetch_from_ncbi(virus_config: dict, global_config: dict, virus_raw_dir: Path):
    """Handles the entire data acquisition process for NCBI-sourced viruses."""
    logging.info("Starting NCBI data acquisition via CLI...")
    
    # Create a temporary directory for the main dataset download
    temp_complete_dir = virus_raw_dir / "temp_complete_download"
    temp_complete_dir.mkdir(exist_ok=True)

    try:
        # 1. Download and process main sequences and metadata
        logging.info("Step 1: Downloading main dataset...")
        ncbi_seq_commands = global_config.get(NCBI_CLI, {}).get(SEQUENCES, [])
        
        if len(ncbi_seq_commands) >= 2:
            # Step 1.1: Run the download command in the 'complete' subdirectory
            download_command_template = ncbi_seq_commands[0]
            download_command = download_command_template.format(
                taxon_id=virus_config[TAXON_ID],
                virus_name=virus_config[NAME]
            )
            run_command(download_command, temp_complete_dir)
            
            # Step 1.2: Unzip the downloaded file
            main_zip_path = temp_complete_dir / f"{virus_config[NAME]}-download.zip"
            logging.info(f"Unzipping {main_zip_path}...")
            with zipfile.ZipFile(main_zip_path, 'r') as zip_ref:
                zip_ref.extractall(temp_complete_dir)

            # Step 1.3: Run the dataformat command
            format_command_template = ncbi_seq_commands[1]
            format_command = format_command_template.format(
                virus_name=virus_config[NAME]
            )
            run_command(format_command, temp_complete_dir)
        else:
            logging.error("NCBI CLI configuration for sequences is incomplete. Expected at least 2 commands.")
            return

        # Move and rename the final standardized files from the 'complete' subdirectory
        unzipped_data_path = temp_complete_dir / "ncbi_dataset" / "data"
        find_and_rename(unzipped_data_path, "genomic.fna", "raw_sequences.fasta", virus_raw_dir)
        find_and_rename(temp_complete_dir, f"{virus_config[NAME]}-raw-metadata.tsv", "raw_metadata.tsv", virus_raw_dir)

    finally:
        logging.info("Step 3: Cleaning up temporary files...")
        if temp_complete_dir.exists():
            shutil.rmtree(temp_complete_dir)
        logging.info("Cleanup complete.")

    # 2. Download and process reference sequence by calling the new helper function
    logging.info("Step 2: Downloading reference sequence...")
    fetch_reference_sequence(virus_config, global_config, virus_raw_dir)

def fetch_from_nextstrain(virus_config: dict, global_config: dict, virus_raw_dir: Path):
    """Handles data acquisition for Nextstrain-sourced viruses via URL."""
    logging.info("Starting Nextstrain data acquisition via URL...")

    # Get the virus name from the config
    virus_name = virus_config.get(NAME)
    if not virus_name:
        logging.error("Virus name not found in the configuration. Cannot proceed with Nextstrain data acquisition.")
        return
    
    # Get the URL block for the specific virus from the global config
    virus_urls = global_config.get(NEXTSTRAIN_URL, {}).get(virus_name, {})
    if not virus_urls:
        logging.error(f"No URLs found for '{virus_name}' under the '{NEXTSTRAIN_URL}' key in the config file.")
        return

    # # 1. Download and decompress metadata
    # metadata_url = virus_urls.get(METADATA)
    # if metadata_url:
    #     logging.info(f"Step 1: Downloading metadata from {metadata_url}")
    #     gz_path = virus_raw_dir / "temp_metadata.tsv.gz"
    #     final_path = virus_raw_dir / "raw_metadata.tsv"
    #     try:
    #         with requests.get(metadata_url, stream=True) as r:
    #             r.raise_for_status()
    #             with open(gz_path, 'wb') as f:
    #                 shutil.copyfileobj(r.raw, f)
            
    #         logging.info(f"Decompressing {gz_path}...")
    #         with gzip.open(gz_path, 'rb') as f_in:
    #             with open(final_path, 'wb') as f_out:
    #                 shutil.copyfileobj(f_in, f_out)
    #         gz_path.unlink()
    #     except Exception as e:
    #         logging.error(f"Failed to download or process metadata: {e}")

    # # 2. Download and decompress sequences
    # sequences_url = virus_urls.get(SEQUENCES)
    # if sequences_url:
    #     logging.info(f"Step 2: Downloading sequences from {sequences_url}")
    #     if sequences_url.endswith(".xz"):
    #         compressed_path = virus_raw_dir / "temp_sequences.fasta.xz"
    #         decompress_func = lzma.open
    #     else:
    #         logging.warning(f"Unrecognized compression format for URL: {sequences_url}. Assuming no compression.")
    #         compressed_path = None

    #     final_path = virus_raw_dir / "raw_sequences.fasta"
    #     try:
    #         with requests.get(sequences_url, stream=True) as r:
    #             r.raise_for_status()
    #             if compressed_path:
    #                 with open(compressed_path, 'wb') as f:
    #                     shutil.copyfileobj(r.raw, f)
                    
    #                 logging.info(f"Decompressing {compressed_path}...")
    #                 with decompress_func(compressed_path) as f_in:
    #                     with open(final_path, 'wb') as f_out:
    #                         shutil.copyfileobj(f_in, f_out)
    #                 compressed_path.unlink()
    #             else:
    #                 with open(final_path, 'wb') as f:
    #                     shutil.copyfileobj(r.raw, f)
    #     except Exception as e:
    #         logging.error(f"Failed to download or process sequences: {e}")

    # 1. Download and decompress metadata
    logging.info("Step 1: Downloading metadata...")
    metadata_url = virus_urls.get(METADATA)
    download_and_decompress(metadata_url, "raw_metadata.tsv", virus_raw_dir)

    if not virus_name == "sars-cov-2":
        # 2. Download and decompress sequences
        logging.info("Step 2: Downloading sequences...")
        sequences_url = virus_urls.get(SEQUENCES)
        download_and_decompress(sequences_url, "raw_sequences.fasta", virus_raw_dir)

        # 3. Download reference sequence from NCBI using its accession ID
        logging.info("Step 3: Fetching reference sequence from NCBI...")
        fetch_reference_sequence(virus_config, global_config, virus_raw_dir)
        
def fetch_from_ftp(virus_config: dict, global_config: dict, virus_raw_dir: Path):
    """Handles data acquisition for FTP-sourced viruses."""
    logging.info("Starting FTP data acquisition...")

    # Get the virus name from the config
    virus_name = virus_config.get(NAME)
    if not virus_name:
        logging.error("Virus name not found in the configuration. Cannot proceed with FTP data acquisition.")
        return

    # Get the FTP URL block for the specific virus from the global config
    ftp_urls = global_config.get(FTP_URL, {}).get(virus_name, {})
    if not ftp_urls:
        logging.error(f"No FTP URLs found for '{virus_name}' under the '{FTP_URL}' key in the config file.")
        return

    # 1. Download and decompress metadata
    metadata_url = ftp_urls.get(METADATA)
    if metadata_url:
        logging.info(f"Step 1: Downloading metadata from {metadata_url}")
        download_and_decompress(metadata_url, "raw_metadata.csv", virus_raw_dir)

def main():
    """Main function to orchestrate data fetching based on config."""
    parser = argparse.ArgumentParser(description="Fetch data for a specified virus based on the pipeline configuration.")
    parser.add_argument("--virus", required=True, help="The name of the virus to fetch data for (e.g., 'yellow-fever').")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    # Load main config
    try:
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"CRITICAL ERROR: Error parsing YAML config file: {e}", file=sys.stderr)
        sys.exit(1)

    # Get config for the specified virus
    virus_name = args.virus
    virus_config = config.get(VIRUSES, {}).get(virus_name)
    if not virus_config:
        print(f"CRITICAL ERROR: No configuration found for virus '{virus_name}' in '{args.config}'", file=sys.stderr)
        sys.exit(1)

    # Setup Logging
    log_base_path = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_base_path, log_name_prefix=f"{virus_name}_01_fetch_data")

    # Determine paths
    raw_data_base_path = Path(config.get(PATHS, {}).get(RAW_DATA, 'data/raw'))
    virus_raw_dir = raw_data_base_path / virus_name
    virus_raw_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Ensuring raw data directory exists: {virus_raw_dir}")

    # Route to the correct fetching function
    try:
        source = virus_config.get(SOURCE).lower()
        if source == NCBI:
            fetch_from_ncbi(virus_config, config, virus_raw_dir)
        elif source == NEXTSTRAIN:
            fetch_from_nextstrain(virus_config, config, virus_raw_dir)
        elif source == FTP:
            fetch_from_ftp(virus_config, config, virus_raw_dir)  
        else:
            logging.error(f"Unknown data acquisition method '{source}' for virus '{virus_name}'.")
            sys.exit(1)
        
        logging.info(f"Data fetching process for '{virus_name}' complete.")
        
        # Update the download_date in config.yaml after successful data acquisition
        logging.info("Updating download_date in config.yaml...")
        update_download_date(args.config, virus_name)
        
    except Exception as e:
        logging.critical(f"A critical error occurred during the fetching process: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
