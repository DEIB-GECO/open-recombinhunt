# src/05_prepare_recombinhunt/create_samples.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import yaml
import json

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    print("Please ensure 'utils.py' and 'constants.py' exist in the 'src/utils' directory.")
    sys.exit(1)

def is_recombinant(lineage_name):
    """Helper function to identify top-level recombinant lineages."""
    if not isinstance(lineage_name, str):
        return False
    lineage_name = lineage_name.upper()
    return lineage_name.startswith('X') and not '.' in lineage_name

def create_samples_step(virus_name: str, config: dict):
    """
    Orchestrates the creation of per-lineage sample files and a JSON summary.
    """
    logging.info(f"--- Starting Step 5.2: Create Samples for '{virus_name}' ---")

    # --- 1. Determine Paths and Parameters ---
    paths_config = config.get(PATHS, {})
    virus_config = config.get(VIRUSES, {}).get(virus_name, {})
    
    # Determine input directory based on HaploCov parameters from config
    haplocov_params = virus_config.get(PARAMETERS, {}).get(HAPLOCOV, {})
    dist = haplocov_params.get(DIST, 0)
    size = haplocov_params.get(SIZE, 0)
    param_string = f"dist{dist}size{size}"
    
    # Determine Input File (Conditional Logic)
    results_dir = Path(paths_config.get(RESULTS))
    if virus_name == 'sars-cov-2':
        logging.info("Virus is 'sars-cov-2'. Using Nextstrain reformatted output as input.")
        # This path might need adjustment based on where the sars-cov-2 processing script saves its output
        input_dir = results_dir / "nextstrain_output" / virus_name
        input_file = input_dir / "nextstrain_reformatted.tsv" 
    else:
        logging.info(f"Virus is '{virus_name}'. Using HaploCoV reformatted output as input.")
        input_dir = results_dir / "haplocov_output" / virus_name / param_string
        input_file = input_dir / "haplocov_reformatted.tsv"
    
    # Determine Output Directory
    samples_dir = Path(paths_config.get(SAMPLES)) / virus_name / param_string
    samples_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Sample files will be saved to: {samples_dir}")

    # --- 2. Load Data ---
    logging.info(f"Loading data from: {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        initial_length = len(df)
    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file}. Please run previous pipeline steps.")
        sys.exit(1)

    # --- 3. Data Cleaning and Reporting ---
    logging.info("Cleaning data by dropping rows with missing pangoLin or mutations...")
    
    # Report and drop rows with missing pangoLin
    missing_pango_mask = df['pangoLin'].isna()
    if missing_pango_mask.any():
        missing_pango_ids = df.loc[missing_pango_mask, 'genomeID'].tolist()
        logging.warning(f"Found {len(missing_pango_ids)} sequences with missing 'pangoLin'. They will be dropped.")
        logging.warning(f"  Dropped IDs (missing pangoLin): {missing_pango_ids}")
        df.dropna(subset=['pangoLin'], inplace=True)

    # Report and drop rows with missing mutations
    missing_mutations_mask = df['mutations'].isna()
    if missing_mutations_mask.any():
        missing_mutations_ids = df.loc[missing_mutations_mask, 'genomeID'].tolist()
        logging.warning(f"Found {len(missing_mutations_ids)} sequences with missing 'mutations'. They will be dropped.")
        logging.warning(f"  Dropped IDs (missing mutations): {missing_mutations_ids}")
        df.dropna(subset=['mutations'], inplace=True)
        
    logging.info(f"Initially, there were {initial_length} rows.")
    logging.info(f"Data cleaning complete. {len(df)} rows remaining for sample creation.")

    # --- MODIFICATION: Add special filtering for SARS-CoV-2 ---
    if virus_name == 'sars-cov-2':
        logging.info("Applying special filter for SARS-CoV-2: Keeping only recombinant lineages.")
        recombinant_mask = df['pangoLin'].apply(is_recombinant)
        df = df[recombinant_mask]
        logging.info(f"{len(df)} rows remaining after keeping only recombinant lineages.")

    if df.empty:
        logging.warning("No data remains after cleaning. No sample files will be created.")
        return

    # --- 4. Create Sample Files and JSON Summary ---
    logging.info("Grouping by lineage and creating sample files...")
    lineage_counts_for_json = {}
    
    # Group by the pangoLin column to process each lineage
    for lineage_name, group_df in df.groupby('pangoLin'):
        
        # Prepare the DataFrame for output
        output_df = group_df[['genomeID', 'pangoLin', 'mutations']].copy()
        output_df.rename(columns={
            'pangoLin': 'true_lineage',
            'mutations': 'nuc_changes'
        }, inplace=True)
        
        # Store count for the JSON summary
        num_sequences = len(output_df)
        lineage_counts_for_json[lineage_name] = num_sequences
        
        # Construct output filename and save as tab-separated
        output_filename = f"samples_{lineage_name}.csv"
        output_filepath = samples_dir / output_filename
        
        try:
            output_df.to_csv(output_filepath, sep='\t', index=False)
            logging.info(f"  Saved: {output_filepath} ({num_sequences} rows)")
        except Exception as e:
            logging.error(f"Failed to save sample file for lineage '{lineage_name}': {e}")
            
    # --- 5. Save JSON Summary File ---
    json_filepath = samples_dir / "samples_total.json"
    logging.info(f"Saving lineage count summary to: {json_filepath}")
    try:
        with open(json_filepath, 'w') as f_json:
            json.dump(lineage_counts_for_json, f_json, indent=4, sort_keys=True)
    except Exception as e:
        logging.error(f"Failed to save JSON summary file: {e}")


def main():
    """Main function to run the sample creation step."""
    parser = argparse.ArgumentParser(description="Create per-lineage sample files for RecombinHunt.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    try:
        config = yaml.safe_load(Path(args.config).read_text())
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_05.2_create_samples")

    try:
        create_samples_step(args.virus, config)
        logging.info("Sample creation step finished successfully.")
    except Exception as e:
        logging.critical(f"Sample creation step failed with a critical error: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
