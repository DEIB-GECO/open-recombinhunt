# src/03_haplocov/run_haplocov.py

import argparse
import yaml
from pathlib import Path
import sys
import logging
import os
import shutil

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging, run_command
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    print("Please ensure 'utils.py' and 'constants.py' exist in the 'src/utils' directory.")
    sys.exit(1)

def run_haplocov_step(virus_name: str, config: dict):
    """
    Orchestrates the execution of the HaploCoV Perl scripts in an isolated environment.
    """
    logging.info(f"--- Starting Step 3: Run HaploCoV for '{virus_name}' ---")
    
    # --- 1. Get Configuration and Paths ---
    paths_config = config.get(PATHS, {})
    virus_config = config.get(VIRUSES, {}).get(virus_name, {})
    haplocov_params = virus_config.get(PARAMETERS, {}).get(HAPLOCOV, {})
    
    dist = haplocov_params.get(DIST)
    size = haplocov_params.get(SIZE)
    designation_mode = haplocov_params.get(DESIGNATION_MODE, BOTH)

    if designation_mode in [NOMENCLATURE]:
        logging.info("Designation mode is set to 'nomenclature'. Setting dist and size to 0.")
        dist = 0
        size = 0
    
    # Define main paths
    processed_dir = Path(paths_config.get(PROCESSED_DATA)) / virus_name
    results_base_dir = Path(paths_config.get(RESULTS))
    haplocov_tool_path = Path(config.get(PATHS).get(LIBS)) / HAPLOCOV

    # Define the specific output directory for this run
    param_string = f"{DIST}{dist}{SIZE}{size}"
    haplocov_output_dir = results_base_dir / "haplocov_output" / virus_name / param_string
    
    # --- 2. Create Temporary Working Directory ---
    # We create it inside the final output directory to keep related files together.
    temp_work_dir = haplocov_output_dir / "temp_work"
    haplocov_output_dir.mkdir(parents=True, exist_ok=True)
    temp_work_dir.mkdir(exist_ok=True)
    logging.info(f"Created temporary working directory: {temp_work_dir}")

    try:
        # --- 3. Prepare Inputs using Symbolic Links ---
        logging.info("Creating symbolic links for input files in temporary directory...")
        
        input_files = {
            "metadata.tsv": processed_dir / "metadata.tsv",
            "sequences.fasta": processed_dir / "sequences.fasta",
            "reference.fasta": processed_dir / "reference.fasta",
            "areaFile": haplocov_tool_path / "areaFile",
            "align.pl": haplocov_tool_path / "align.pl"
        }
        
        for dest_name, src_path in input_files.items():
            dest_path = temp_work_dir / dest_name
            if not src_path.exists():
                logging.error(f"Required input file not found: {src_path}")
                raise FileNotFoundError(f"Missing input for HaploCoV: {src_path}")
            
            # Create symlink if it doesn't already exist
            if not dest_path.exists():
                os.symlink(src_path.resolve(), dest_path)

        # --- 4. Define and Execute Perl Commands based on Designation Mode ---
        logging.info(f"Executing HaploCoV with designation_mode: '{designation_mode}'")
        
        # Define the relative path to the Perl scripts from the temp_work_dir
        # This makes the command robust regardless of where the project is located.
        relative_perl_path = os.path.relpath(haplocov_tool_path, temp_work_dir)
        
        # Define the sequence of commands
        commands = [
            f"perl {relative_perl_path}/addToTableNCBI.pl --metadata metadata.tsv --ref reference.fasta --seq sequences.fasta --outfile out.HaploCoV --nproc 16",
            f"perl {relative_perl_path}/computeDefining.pl out.HaploCoV",
            f"perl {relative_perl_path}/assign.pl --dfile out.HaploCoV.definingVariants.txt --infile out.HaploCoV --outfile out.HaploCoV.assigned",
            f"perl {relative_perl_path}/augmentClusters.pl --metafile out.HaploCoV.assigned --deffile out.HaploCoV.definingVariants.txt --posFile out.HaploCoV.frequentVariants.txt --outfile out.HaploCoV.definingVariantsNew.txt --dist {dist} --size {size}",
            f"perl {relative_perl_path}/assign.pl --dfile out.HaploCoV.definingVariantsNew.txt --infile out.HaploCoV --outfile out.HaploCoV.assignedNew"
        ]

        if designation_mode in [NOMENCLATURE]:
            commands_to_run = commands[:1].copy() 
            final_output_file = "out.HaploCoV"

        elif designation_mode in [HAPLOCOV, BOTH]:
            commands_to_run = commands.copy()
            final_output_file = "out.HaploCoV.assignedNew"

        commands_to_run.append(
                f"cut -f 10 {final_output_file} | sort | uniq -c > {virus_name}-report.txt"
            )

        for command in commands_to_run:
            run_command(command, working_dir=temp_work_dir)
            
        # --- 5. Move Final Output ---
        final_output_path = temp_work_dir / final_output_file
        if final_output_path.exists():
            destination_path = haplocov_output_dir / "haplocov_assigned.tsv"
            logging.info(f"Moving final output file '{final_output_file}' to '{destination_path}'")
            shutil.move(final_output_path, destination_path)
        else:
            logging.error(f"Expected final output file '{final_output_file}' was not created by HaploCoV.")

    finally:
        logging.info("HaploCoV step completed.")
        # --- 6. Cleanup ---
        # logging.info("Cleaning up temporary working directory...")
        # if temp_work_dir.exists():
        #     shutil.rmtree(temp_work_dir)
        # logging.info("Cleanup complete.")


def main():
    """Main function to run the HaploCoV step."""
    parser = argparse.ArgumentParser(description="Run the HaploCoV tool on preprocessed data.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    try:
        config = yaml.safe_load(Path(args.config).read_text())
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_03_run_haplocov")

    try:
        run_haplocov_step(args.virus, config)
        logging.info("HaploCoV step finished successfully.")
    except Exception as e:
        logging.critical(f"HaploCoV step failed with a critical error: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
