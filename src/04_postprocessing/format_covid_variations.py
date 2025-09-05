# src/04_postprocessing/format_covid_variations.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import yaml
import re

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    sys.exit(1)

# --- Transformation Helper Functions ---

def _transform_deletion_part(part: str) -> str:
    """Transforms a single deletion entry (e.g., '11288-11297' -> '11288_11297')."""
    return part.replace('-', '_')

def _transform_insertion_part(part: str) -> str:
    """Transforms a single insertion entry (e.g., '28263:A' -> '28263_.|A')."""
    return part.replace(':', '_.|')

def _transform_substitution_part(part: str) -> str:
    """Transforms a single substitution entry using regex (e.g., 'C241T' -> '241_C|T')."""
    try:
        match = re.fullmatch(r'([A-Z-]+)(\d+)([A-Z-]+)', part, re.IGNORECASE)
        if match:
            base1, number, base2 = match.groups()
            return f"{number}_{base1}|{base2}"
    except Exception as e:
        logging.warning(f"Error processing substitution part '{part}': {e}")    
    return part

# --- Main Processing Functions ---

def process_mutations(df: pd.DataFrame) -> pd.DataFrame:
    """
    Processes the 'substitutions', 'deletions', and 'insertions' columns,
    transforms them, and combines them into a single 'mutations' column.
    """
    logging.info("Starting transformation of mutation columns...")
    df_processed = df.copy()

    def transform_column(column_series: pd.Series, part_transformer_func) -> pd.Series:
        series_str = column_series.fillna('').astype(str)
        def process_cell(cell_value: str):
            if not cell_value: return ""
            parts = cell_value.split(',')
            transformed_parts = [part_transformer_func(part.strip()) for part in parts if part.strip()]
            return ",".join(transformed_parts)
        return series_str.apply(process_cell)

    mutation_source_cols = ['substitutions', 'deletions', 'insertions']
    transformed_cols = []

    for col in mutation_source_cols:
        transformed_col_name = f"{col}_transformed"
        transformed_cols.append(transformed_col_name)
        if col in df_processed.columns:
            logging.info(f"Processing '{col}'...")
            if col == 'substitutions': transformer = _transform_substitution_part
            elif col == 'deletions': transformer = _transform_deletion_part
            else: transformer = _transform_insertion_part
            df_processed[transformed_col_name] = transform_column(df_processed[col], transformer)
        else:
            logging.warning(f"Column '{col}' not found. An empty transformed column will be created.")
            df_processed[transformed_col_name] = ''

    logging.info("Combining transformed columns into 'mutations' column...")
    df_processed['mutations'] = df_processed[transformed_cols].apply(
        lambda row: ",".join(filter(None, row)),
        axis=1
    )
    
    logging.info("Mutation processing complete.")
    return df_processed

# --- Validation Function (from your notebook) ---

def validate_formats(df: pd.DataFrame):
    """Runs assertions to validate the format of intermediate transformed columns."""
    logging.info("Validating intermediate column formats...")
    
    SUBSTITUTION_PATTERN = re.compile(r"^\d+_[A-Z-]+\|[A-Z-]+$", re.IGNORECASE)
    DELETION_PATTERN = re.compile(r"^\d+(_\d+)?$")
    INSERTION_PATTERN = re.compile(r"^\d+_\.\|[A-Z-]+$", re.IGNORECASE)

    patterns = {
        'substitutions_transformed': SUBSTITUTION_PATTERN,
        'deletions_transformed': DELETION_PATTERN,
        'insertions_transformed': INSERTION_PATTERN
    }

    def validate_column_format(series: pd.Series, pattern: re.Pattern) -> bool:
        series_to_check = series[series != '']
        if series_to_check.empty: return True
        
        for index, mutation_list_str in series_to_check.items():
            for mutation in mutation_list_str.split(','):
                if not pattern.fullmatch(mutation.strip()):
                    logging.error(f"Invalid format in row index {index}, column '{series.name}'.")
                    logging.error(f"  -> Problematic value: '{mutation}'")
                    return False
        return True

    try:
        assert validate_column_format(df['substitutions_transformed'], patterns['substitutions_transformed']), "Substitution format validation failed!"
        assert validate_column_format(df['deletions_transformed'], patterns['deletions_transformed']), "Deletion format validation failed!"
        assert validate_column_format(df['insertions_transformed'], patterns['insertions_transformed']), "Insertion format validation failed!"
        logging.info("All format assertions passed successfully!")
    except AssertionError as e:
        logging.critical(f"VALIDATION FAILED: {e}")
        sys.exit(1)

# --- Main Execution Block ---

def main():
    parser = argparse.ArgumentParser(description="Post-process Nextstrain SARS-CoV-2 metadata to create a standardized 'mutations' column.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process (should be 'sars-cov-2').")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    if args.virus != 'sars-cov-2':
        print(f"Error: This script is specifically designed for 'sars-cov-2'. You provided '{args.virus}'.")
        sys.exit(1)

    config = yaml.safe_load(Path(args.config).read_text())
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_04_format_covid_variations")

    # --- Get Paths ---
    paths_config = config.get(PATHS, {})
    processed_data_dir = Path(paths_config.get(PROCESSED_DATA, 'data/processed')) / args.virus
    results_dir = Path(paths_config.get(RESULTS, 'results'))
    
    input_file = processed_data_dir / "metadata.tsv"
    
    # Define the specific output directory for this run
    output_dir = results_dir / "nextstrain_output" / args.virus
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "nextstrain_reformatted.tsv"
    output_file_last_6_months = output_dir / "nextstrain_reformatted_last_6_months.tsv"

    # --- Load and Process Data ---
    logging.info(f"Loading processed metadata from: {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file}. Please run the preprocessing scripts first.")
        sys.exit(1)
        
    # Run the main processing function from your notebook
    df_with_mutations = process_mutations(df)
    
    # Run validation on the intermediate columns
    validate_formats(df_with_mutations)
    
    # Drop intermediate and original mutation columns to create the final DataFrame
    logging.info("Dropping intermediate and original mutation columns...")
    cols_to_drop = [
        'substitutions_transformed', 'deletions_transformed', 'insertions_transformed',
        'substitutions', 'deletions', 'insertions'
    ]
    # Drop only the columns that actually exist
    existing_cols_to_drop = [col for col in cols_to_drop if col in df_with_mutations.columns]
    df_final = df_with_mutations.drop(columns=existing_cols_to_drop)
    
    # --- MODIFICATION: Rename 'pango_lineage' to 'pangoLin' ---
    if 'pango_lineage' in df_final.columns:
        logging.info("Renaming 'pango_lineage' column to 'pangoLin' for consistency...")
        df_final.rename(columns={'pango_lineage': 'pangoLin'}, inplace=True)

    if 'Pango lineage' in df_final.columns:
        logging.info("Renaming 'Pango lineage' column to 'pangoLin' for consistency...")
        df_final.rename(columns={'Pango lineage': 'pangoLin'}, inplace=True)

    # --- MODIFICATION: Rename 'Virus name' to genomeID ---
    if 'Virus name' in df_final.columns:
        logging.info("Renaming 'Virus name' column to 'genomeID' for consistency...")
        df_final.rename(columns={'Virus name': 'genomeID'}, inplace=True)
    
    logging.info(f"Saving reformatted data with {len(df_final)} rows to: {output_file}")
    df_final.to_csv(output_file, sep='\t', index=False)

    current_date = pd.Timestamp.now()
    six_months_ago = current_date - pd.DateOffset(months=6)
    df_last_6_months = df_final[df_final['Collection date'] >= six_months_ago.strftime('%Y-%m-%d')]
    logging.info(f"Saving reformatted data from the last 6 months with {len(df_last_6_months)} rows to: {output_file_last_6_months}")
    df_last_6_months.to_csv(output_file_last_6_months, sep='\t', index=False)

    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()
