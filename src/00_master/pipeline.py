# src/00_master/pipeline.py

import argparse
import yaml
from pathlib import Path
import sys
import logging
import subprocess

# Add the project's 'src' directory to the Python path.
# This allows us to import modules like 'utils' and 'constants'
SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    print("Please ensure 'utils.py' and 'constants.py' exist in the 'src/utils' directory.")
    sys.exit(1)

def run_pipeline_step(command: list):
    """
    Executes a pipeline step as a subprocess and handles errors.
    The command should be a list of arguments, e.g., ['python3', 'script.py', '--arg', 'value']
    """
    command_str = " ".join(command)
    logging.info(f"--- Executing: {command_str} ---")
    try:
        # Using capture_output=True to get stdout/stderr for logging
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True
        )
        if result.stdout:
            logging.info(f"Stdout from step:\n{result.stdout}")
        if result.stderr:
            logging.warning(f"Stderr from step:\n{result.stderr}")
        logging.info(f"--- Step Succeeded: {command_str} ---")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"--- Step FAILED: {command_str} ---")
        logging.error(f"Return Code: {e.returncode}")
        logging.error(f"Stderr:\n{e.stderr}")
        logging.error(f"Stdout:\n{e.stdout}")
        return False

def main():
    """Main orchestrator for the OpenRecombinHunt pipeline."""
    parser = argparse.ArgumentParser(description="Main orchestrator for the OpenRecombinHunt pipeline.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process, or 'all' to process all viruses in the config.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    # Load main config
    try:
        config_path = Path(args.config)
        config = yaml.safe_load(config_path.read_text())
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"CRITICAL ERROR: Error parsing YAML config file: {e}", file=sys.stderr)
        sys.exit(1)

    # Determine which viruses to process
    if args.virus.lower() == ALL:
        viruses_to_process = list(config.get(VIRUSES, {}).keys())
        if not viruses_to_process:
            print("No viruses found in the configuration file.")
            sys.exit(0)
    else:
        if args.virus not in config.get(VIRUSES, {}):
            print(f"Error: Virus '{args.virus}' not found in the configuration file.")
            sys.exit(1)
        viruses_to_process = [args.virus]

    # --- Main Loop: Process each virus ---
    for virus_name in viruses_to_process:
        # Setup Logging for this entire virus run
        log_dir = Path(config.get(PATHS, {}).get(LOGS))
        setup_logging(log_dir=log_dir, log_name_prefix=f"{virus_name}_pipeline_run")
        
        logging.info(f"================== STARTING PIPELINE FOR: {virus_name.upper()} ==================")
        
        virus_config = config.get(VIRUSES, {}).get(virus_name, {})
        source = virus_config.get(SOURCE)

        # --- Define the sequence of scripts to run ---
        # The command is built as a list for subprocess.run
        python_executable = sys.executable # Use the same python that is running the pipeline
        
        # Step 1: Data Acquisition
        step1_fetch = [python_executable, "src/01_data_acquisition/fetch_data.py", "--virus", virus_name, "--config", args.config]
        
        # Step 2: Preprocessing (Conditional based on source)
        if source == "ncbi":
            step2_prep_meta = [python_executable, "src/02_preprocessing/ncbi/prep_metadata_ncbi.py", "--virus", virus_name, "--config", args.config]
            step2_prep_fasta = [python_executable, "src/02_preprocessing/common/prep_fasta.py", "--virus", virus_name, "--config", args.config]
            step2_prep_ref = [python_executable, "src/02_preprocessing/common/prep_ref.py", "--virus", virus_name, "--config", args.config]
        elif source == "nextstrain":
            step2_prep_meta = [python_executable, "src/02_preprocessing/nextstrain/prep_metadata_nextstrain.py", "--virus", virus_name, "--config", args.config]
            step2_prep_fasta = [python_executable, "src/02_preprocessing/common/prep_fasta.py", "--virus", virus_name, "--config", args.config]
            step2_prep_ref = [python_executable, "src/02_preprocessing/common/prep_ref.py", "--virus", virus_name, "--config", args.config]
        elif source == "ftp":
            step2_prep_meta = [python_executable, "src/02_preprocessing/ftp/prep_metadata_ftp.py", "--virus", virus_name, "--config", args.config]
            step2_prep_fasta = [python_executable, "src/02_preprocessing/ftp/fetch_ftp_sequences_by_id.py", "--virus", virus_name, "--config", args.config]
            step2_prep_ref = [python_executable, "src/02_preprocessing/ftp/fetch_ftp_reference.py", "--virus", virus_name, "--config", args.config]
        else:
            logging.error(f"Unknown source '{source}' for virus '{virus_name}'. Cannot determine preprocessing script.")
            continue # Skip to the next virus

        # Step 3: Run HaploCoV
        step3_run_haplocov = [python_executable, "src/03_haplocov/run_haplocov.py", "--virus", virus_name, "--config", args.config]

        # Step 4: Post-processing
        step4_postprocess_haplocov = [python_executable, "src/04_postprocessing/format_haplocov_variations.py", "--virus", virus_name, "--config", args.config]
        step4_postprocess_covid = [python_executable, "src/04_postprocessing/format_covid_variations.py", "--virus", virus_name, "--config", args.config]

        # Step 5: Prepare for RecombinHunt
        step5_create_env = [python_executable, "src/05_prepare_recombinhunt/create_environment.py", "--virus", virus_name, "--config", args.config]
        step5_create_samples = [python_executable, "src/05_prepare_recombinhunt/create_samples.py", "--virus", virus_name, "--config", args.config]

        # Step 6: Run RecombinHunt
        step6_run_recombinhunt = [python_executable, "src/06_recombinhunt/run_recombinhunt.py", "--virus", virus_name, "--config", args.config]

        # --- Execute the pipeline sequentially ---
        if virus_name.lower() == "sars-cov-2":
            pipeline_steps = [
                #step1_fetch,
                #step2_prep_meta,
                #step4_postprocess_covid,
                #step5_create_env,
                step5_create_samples,
                step6_run_recombinhunt
            ]
        else:
            pipeline_steps = [
                step1_fetch,
                step2_prep_meta,
                step2_prep_fasta,
                step2_prep_ref,
                step3_run_haplocov,
                step4_postprocess_haplocov,
                step5_create_env,
                step5_create_samples,
                step6_run_recombinhunt
            ]
            
        for step_command in pipeline_steps:
            success = run_pipeline_step(step_command)
            if not success:
                logging.critical(f"Pipeline for '{virus_name}' failed at step: {' '.join(step_command)}")
                logging.critical("Aborting pipeline for this virus.")
                break # Stop processing this virus and move to the next one
        
        logging.info(f"================== PIPELINE FOR {virus_name.upper()} FINISHED ==================\n\n")

if __name__ == "__main__":
    main()
# This script is the main entry point for the OpenRecombinHunt pipeline.
# It orchestrates the execution of various steps based on the provided configuration.