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

        # Get analysis window from config for dynamic filename
        try:
            virus_config = config.get(VIRUSES).get(virus_name)
            analysis_window_months = virus_config.get("analysis_window_months", 6)
        except Exception:
            analysis_window_months = 6
        
        logging.info(f"Loading data from the last {analysis_window_months} months...")
        input_file_l6m = input_dir / f"nextstrain_reformatted_last_{analysis_window_months}_months.tsv"
    else:
        logging.info(f"Virus is '{virus_name}'. Using HaploCoV reformatted output as input.")
        input_dir = results_dir / "haplocov_output" / virus_name / param_string
        input_file = input_dir / "haplocov_reformatted.tsv"
    
    # Determine Output Directory
    samples_dir = Path(paths_config.get(SAMPLES)) / virus_name / param_string
    samples_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Sample files will be saved to: {samples_dir}")

    if virus_name == 'sars-cov-2':
        # Get analysis window from config for dynamic directory naming
        try:
            virus_config = config.get(VIRUSES).get(virus_name)
            analysis_window_months = virus_config.get("analysis_window_months", 6)
        except Exception:
            analysis_window_months = 6
        
        samples_dir_l6m = Path(paths_config.get(SAMPLES)) / virus_name / param_string / f"last_{analysis_window_months}_months"
        samples_dir_l6m.mkdir(parents=True, exist_ok=True)
        logging.info(f"FOR last {analysis_window_months} months, Sample files will be saved to: {samples_dir_l6m}")

    # --- 2. Load Data ---
    logging.info(f"Loading data from: {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
        initial_length = len(df)

        if virus_name == 'sars-cov-2':
            df_l6m = pd.read_csv(input_file_l6m, sep='\t', low_memory=False)
            initial_length_l6m = len(df_l6m)
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

    if virus_name == 'sars-cov-2':
        missing_pango_mask_l6m = df_l6m['pangoLin'].isna()
        if missing_pango_mask_l6m.any():
            missing_pango_ids_l6m = df_l6m.loc[missing_pango_mask_l6m, 'genomeID'].tolist()
            logging.warning(f"Found {len(missing_pango_ids_l6m)} sequences with missing 'pangoLin'. They will be dropped.")
            logging.warning(f"  Dropped IDs (missing pangoLin): {missing_pango_ids_l6m}")
            df_l6m.dropna(subset=['pangoLin'], inplace=True)

    # Report and drop rows with missing mutations
    missing_mutations_mask = df['mutations'].isna()
    if missing_mutations_mask.any():
        missing_mutations_ids = df.loc[missing_mutations_mask, 'genomeID'].tolist()
        logging.warning(f"Found {len(missing_mutations_ids)} sequences with missing 'mutations'. They will be dropped.")
        logging.warning(f"  Dropped IDs (missing mutations): {missing_mutations_ids}")
        df.dropna(subset=['mutations'], inplace=True)

    if virus_name == 'sars-cov-2':
        missing_mutations_mask_l6m = df_l6m['mutations'].isna()
        if missing_mutations_mask_l6m.any():
            missing_mutations_ids_l6m = df_l6m.loc[missing_mutations_mask_l6m, 'genomeID'].tolist()
            logging.warning(f"Found {len(missing_mutations_ids_l6m)} sequences with missing 'mutations'. They will be dropped.")
            logging.warning(f"  Dropped IDs (missing mutations): {missing_mutations_ids_l6m}")
            df_l6m.dropna(subset=['mutations'], inplace=True)
        
    logging.info(f"Initially, there were {initial_length} rows.")
    logging.info(f"Data cleaning complete. {len(df)} rows remaining for sample creation.")

    if virus_name == 'sars-cov-2':
        logging.info(f"Initially, there were {initial_length_l6m} rows.")
        logging.info(f"Data cleaning complete. {len(df_l6m)} rows remaining for sample creation.")

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
        if virus_name == 'sars-cov-2':
            logging.info(f"Processing lineage: {lineage_name} (SARS-CoV-2 special handling)")
            # + Collection date
            output_df = group_df[['genomeID', 'pangoLin', 'mutations', 'Collection date']].copy()
            output_df.rename(columns={
                'pangoLin': 'true_lineage',
                'mutations': 'nuc_changes',
                'Collection date': 'collection_date'
            }, inplace=True)

        else:
            logging.info(f"Processing lineage: {lineage_name}")
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

    lineage_counts_for_json_l6m = {}
    if virus_name == 'sars-cov-2':
        for lineage_name_l6m, group_df_l6m in df_l6m.groupby('pangoLin'):
            # Prepare the DataFrame for output
            output_df_l6m = group_df_l6m[['genomeID', 'pangoLin', 'mutations', 'Collection date']].copy()
            output_df_l6m.rename(columns={
                'pangoLin': 'true_lineage',
                'mutations': 'nuc_changes',
                'Collection date': 'collection_date'
            }, inplace=True)
            num_sequences_l6m = len(output_df_l6m)
            lineage_counts_for_json_l6m[lineage_name_l6m] = num_sequences_l6m
            output_filename_l6m = f"samples_{lineage_name_l6m}.csv"
            output_filepath_l6m = samples_dir_l6m / output_filename_l6m
            try:
                output_df_l6m.to_csv(output_filepath_l6m, sep='\t', index=False)
                logging.info(f"  Saved: {output_filepath_l6m} ({num_sequences_l6m} rows)")
            except Exception as e:
                logging.error(f"Failed to save sample file for lineage '{lineage_name_l6m}': {e}")
            
    # --- 5. Save JSON Summary File ---
    json_filepath = samples_dir / "samples_total.json"
    logging.info(f"Saving lineage count summary to: {json_filepath}")
    try:
        with open(json_filepath, 'w') as f_json:
            json.dump(lineage_counts_for_json, f_json, indent=4, sort_keys=True)
    except Exception as e:
        logging.error(f"Failed to save JSON summary file: {e}")

    if virus_name == 'sars-cov-2':
        json_filepath_l6m = samples_dir_l6m / "samples_total.json"
        logging.info(f"Saving lineage count summary to: {json_filepath_l6m}")
        try:
            with open(json_filepath_l6m, 'w') as f_json:
                json.dump(lineage_counts_for_json_l6m, f_json, indent=4, sort_keys=True)
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
