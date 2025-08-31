# src/02_preprocessing/nextstrain/prep_metadata_nextstrain.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import re
import yaml
from typing import Tuple, List, Set

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    sys.exit(1)

# --- Constants ---
# The primary ID column in the raw Nextstrain metadata that maps to our standard "Virus name"
SPECIAL_CHARS_PATTERN = re.compile(r'[&/\\:]')

# Hardcoded column mapping for Nextstrain data
# Defines the final output columns and their source.
# Location and Pango lineage are handled dynamically in the processing function.
COLUMN_MAPPING = {
    "accession": "Virus name",
    "date": "Collection date",
    "date_submitted": "Submission date"
}

# --- Formatting Functions ---

def _format_single_date(date_str: str) -> str:
    """
    Helper function to format a single date string according to specific rules.
    - YYYY-MM-DD -> YYYY-MM-DD
    - YYYY-MM    -> YYYY-MM-01
    - YYYY       -> YYYY-01-01
    - 'XX' in YYYY-XX-XX or YYYY-MM-XX -> replaced with '01'
    Returns pd.NA for invalid or unparseable formats.
    """
    if pd.isna(date_str):
        return pd.NA
    
    date_str = str(date_str).strip()
    date_str = date_str.split('T')[0] # Remove time part if present

    date_str = date_str.upper().replace('XX', '01')
    
    # Check for YYYY-MM-DD format
    if re.fullmatch(r'\d{4}-\d{2}-\d{2}', date_str):
        return date_str
    # Check for YYYY-MM format
    elif re.fullmatch(r'\d{4}-\d{2}', date_str):
        return f"{date_str}-01"
    # Check for YYYY format
    elif re.fullmatch(r'\d{4}', date_str):
        return f"{date_str}-01-01"
    else:
        # As a fallback, try pandas to_datetime for other potential formats
        # If it fails, return NA
        try:
            return pd.to_datetime(date_str).strftime('%Y-%m-%d')
        except (ValueError, TypeError):
            return pd.NA

def format_date_columns(df: pd.DataFrame) -> Tuple[pd.DataFrame, Set[str]]:
    """
    Applies standardized date formatting to the relevant columns and logs format statistics.
    """
    logging.info("Formatting date columns...")
    
    date_cols_to_process = ['Collection date', 'Submission date']
    dropped_ids = set()
    df_formatted = df.copy()
    
    for col_name in date_cols_to_process:
        if col_name in df_formatted.columns:
            logging.info(f"  - Analyzing and formatting '{col_name}'")

            # First, identify and drop rows where the date starts with 'XXXX'
            # Ensure the column is treated as string for this check
            is_xxxx_mask = df_formatted[col_name].astype(str).str.startswith('XXXX', na=False)
            if is_xxxx_mask.any():
                xxxx_dropped_ids = set(df_formatted.loc[is_xxxx_mask, 'Virus name'].dropna())
                dropped_ids.update(xxxx_dropped_ids)
                df_formatted = df_formatted[~is_xxxx_mask].copy()
                logging.info(f"    Removed {len(xxxx_dropped_ids)} rows with incomplete 'XXXX' dates.")
            
            # --- New Logging Logic ---
            # Work with a string representation of the series for pattern matching
            series_str = df_formatted[col_name].astype(str)
            
            # Count initial NAs by checking for pandas NA and common null-like strings
            na_mask = df_formatted[col_name].isna() | series_str.str.lower().isin(['nan', 'na', 'nat', '<na>', ''])
            na_count = na_mask.sum()
            
            # Analyze only the non-NA part of the series
            non_na_series = series_str[~na_mask]
            
            ymd_count = non_na_series.str.fullmatch(r'\d{4}-\d{2}-\d{2}', na=False).sum()
            ym_count = non_na_series.str.fullmatch(r'\d{4}-\d{2}', na=False).sum()
            y_count = non_na_series.str.fullmatch(r'\d{4}', na=False).sum()
            
            # Count how many will fall into the 'else' block (unparseable or other formats)
            else_count = len(non_na_series) - (ymd_count + ym_count + y_count)

            logging.info(f"    Date format statistics for '{col_name}':")
            logging.info(f"      YYYY-MM-DD format: {ymd_count}")
            logging.info(f"      YYYY-MM format   : {ym_count}")
            logging.info(f"      YYYY format      : {y_count}")
            logging.info(f"      Other/unparseable: {else_count}")
            logging.info(f"      Missing (NA)     : {na_count}")
            # --- End of New Logging Logic ---

            # Apply the transformation to each cell in the column
            df_formatted[col_name] = df_formatted[col_name].apply(_format_single_date)
            
    return df_formatted, dropped_ids

def format_location_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Combines 'region' and 'country' into a single 'Location' column.
    The final format is 'region / country'. Handles missing values.
    """
    logging.info("Formatting 'Location' column by combining 'region' and 'country'...")
    
    region_col = 'region'
    country_col = 'country'

    # Handle potentially missing columns by creating empty Series if needed
    region_series = df[region_col].astype(str).fillna('') if region_col in df.columns else pd.Series([''] * len(df), index=df.index)
    country_series = df[country_col].astype(str).fillna('') if country_col in df.columns else pd.Series([''] * len(df), index=df.index)
    country_series = country_series.str.split(':', n=1).str[0].str.strip()

    # Combine using apply
    def combine_location(row):
        region = row['region']
        country = row['country']
        if region and country:
            return f"{region} / {country}"
        elif region:
            return region
        elif country:
            return country
        else:
            return pd.NA # Use pandas NA for missing location

    # Create the new 'Location' series using your logic
    location_series = pd.DataFrame({'region': region_series, 'country': country_series}).apply(combine_location, axis=1)
    
    # Assign the new series to the 'Location' column of the DataFrame
    df['Location'] = location_series

    COLUMN_MAPPING['Location'] = 'Location' 
    
    return df

# --- Main Processing Function ---

def process_nextstrain_metadata(df: pd.DataFrame, filter_rules: list, source_lineage_col: str, column_mapping: dict, isCovid: bool) -> Tuple[pd.DataFrame, Set[str], Set[str]]:
    """Applies a full preprocessing workflow to the raw Nextstrain metadata."""
    df_processed = df.copy()

    logging.info(f"Starting preprocessing of Nextstrain metadata with {len(df_processed)} rows.")
    
    # Determine the source ID column from the mapping, with a fallback
    source_id_col = next((k for k, v in column_mapping.items() if v == 'Virus name'))
    if source_id_col not in df_processed.columns:
        logging.error(f"Source ID column '{source_id_col}' not found in the data. Exiting.")
        sys.exit(1)

    filtered_ids = set()
    non_classified_ids = set()

    # --- 1. Apply Quality Filters from Config ---
    logging.info("Step 1: Applying quality filters from config...")
    for rule in filter_rules:
        col, op, val = rule.get('column'), rule.get('operator'), rule.get('value')
        if col not in df_processed.columns:
            logging.warning(f"Filter rule column '{col}' not found. Skipping filter.")
            continue
        
        initial_rows = len(df_processed)
        if op == "notna": keep_mask = df_processed[col].notna()
        elif op == ">=": numeric_col = pd.to_numeric(df_processed[col], errors='coerce'); keep_mask = numeric_col >= val
        elif op == ">": numeric_col = pd.to_numeric(df_processed[col], errors='coerce'); keep_mask = numeric_col > val
        elif op == "<=": numeric_col = pd.to_numeric(df_processed[col], errors='coerce'); keep_mask = numeric_col <= val
        elif op == "<": numeric_col = pd.to_numeric(df_processed[col], errors='coerce'); keep_mask = numeric_col < val
        elif op == "==": keep_mask = df_processed[col] == val
        elif op == "!=": keep_mask = df_processed[col] != val
        else: logging.warning(f"Filter operator '{op}' not implemented. Skipping."); continue
            
        dropped_mask = ~keep_mask
        ids_to_drop_this_step = set(df_processed.loc[dropped_mask, source_id_col].dropna())
        filtered_ids.update(ids_to_drop_this_step)
        df_processed = df_processed[keep_mask]
        logging.info(f"  Filter '{col} {op} {val}': Removed {initial_rows - len(df_processed)} rows.")
        
    if not isCovid:
        # --- 2. Special Character Filter on ID ---
        logging.info("Step 2: Checking for special characters in ID...")
        initial_rows = len(df_processed)
        special_char_mask = df_processed[source_id_col].astype(str).str.contains(SPECIAL_CHARS_PATTERN, regex=True, na=True)
        ids_to_drop_special_chars = set(df_processed.loc[special_char_mask, source_id_col].dropna())
        filtered_ids.update(ids_to_drop_special_chars)
        df_processed = df_processed[~special_char_mask]
        logging.info(f"  Filter 'Special Chars in ID': Removed {initial_rows - len(df_processed)} rows.")
    
    # --- 3. Handle Pango Lineage ---
    logging.info("Step 3: Processing lineage column...")
    if source_lineage_col and source_lineage_col.upper() != 'NONE':
        logging.info(f"Using source column '{source_lineage_col}' for lineage information.")
        if source_lineage_col in df_processed.columns:
            lineage_na_mask = df_processed[source_lineage_col].isna()
            ids_to_drop_no_lineage = set(df_processed.loc[lineage_na_mask, source_id_col].dropna())
            non_classified_ids.update(ids_to_drop_no_lineage)
            df_processed = df_processed[~lineage_na_mask]
            column_mapping[source_lineage_col] = 'Pango lineage'
            logging.info(f"  Identified and removed {len(ids_to_drop_no_lineage)} non-classified rows (missing '{source_lineage_col}').")
        else:
            logging.error(f"Source lineage column '{source_lineage_col}' not found. Exiting.")
            sys.exit(1)
    else:
        logging.info("No source lineage column specified ('NONE'). Assigning default lineage 'A.1'.")
        df_processed['Pango lineage'] = 'A.1'

    # --- 4. Final Formatting ---
    logging.info("Step 4: Performing final formatting...")
    
    # Rename columns first to get standard names like 'Collection date'
    df_processed = df_processed.rename(columns=column_mapping)
    
    # Format dates (includes dropping 'XXXX' dates and adding those IDs to filtered_ids)
    df_processed, xxxx_dropped_ids = format_date_columns(df_processed)
    filtered_ids.update(xxxx_dropped_ids)
    
    # Format location
    df_processed = format_location_column(df_processed)
    
    # Select and order final columns
    if isCovid: expected_final_columns = ["Virus name", "Collection date", "Submission date", "Location", "Pango lineage", "substitutions", "deletions", "insertions"]
    else:       expected_final_columns = ["Virus name", "Collection date", "Submission date", "Location", "Pango lineage"]
    
    final_columns_to_keep = [col for col in expected_final_columns if col in df_processed.columns]
    df_final = df_processed[final_columns_to_keep]
    
    logging.info(f"Preprocessing complete. Final dataset has {len(df_final)} rows.")
    
    return df_final, sorted(list(filtered_ids)), sorted(list(non_classified_ids))


def main():
    """Main function to run the Nextstrain metadata preprocessing."""
    parser = argparse.ArgumentParser(description="Preprocess raw Nextstrain metadata based on rules in config.yaml.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    config = yaml.safe_load(Path(args.config).read_text())
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_02.1_prep_metadata_nextstrain")

    virus_config = config.get(VIRUSES, {}).get(args.virus, {})
    paths_config = config.get(PATHS, {})
    
    raw_dir = Path(paths_config.get(RAW_DATA, 'data/raw')) / args.virus
    processed_dir = Path(paths_config.get(PROCESSED_DATA, 'data/processed')) / args.virus
    processed_dir.mkdir(parents=True, exist_ok=True)
    
    input_file = raw_dir / "raw_metadata.tsv"
    output_file = processed_dir / "metadata.tsv"
    filtered_ids_file = processed_dir / "filtered_sequence_ids.txt"
    non_classified_ids_file = processed_dir / "non_classified_sequence_ids.txt"

    logging.info(f"Loading raw metadata from {input_file}...")
    try:
        df_raw = pd.read_csv(input_file, sep='\t', low_memory=False)
    except FileNotFoundError:
        logging.error(f"Raw metadata file not found. Please run '01_fetch_data' first.")
        sys.exit(1)

    processing_params = virus_config.get(PARAMETERS, {}).get(METADATA_PROCESSING, {})
    filter_rules = processing_params.get(FILTERS, [])
    source_lineage_col = processing_params.get(SOURCE_LINEAGE_COLUMN, 'NONE')

    isCovid = False
    if virus_config.get(NAME) == 'sars-cov-2': 
        # rename accession to strain in COLUMN_MAPPING
        COLUMN_MAPPING['strain'] = COLUMN_MAPPING.pop('accession', 'Virus name')
        isCovid = True
    
    # Pass a copy of the mapping constant to the processing function
    df_processed, filtered_ids, non_classified_ids = process_nextstrain_metadata(df_raw, filter_rules, source_lineage_col, COLUMN_MAPPING.copy(), isCovid=isCovid)
    
    logging.info(f"Saving processed metadata ({len(df_processed)} rows) to {output_file}...")
    df_processed.to_csv(output_file, sep='\t', index=False, na_rep='NA')

    logging.info(f"Saving {len(filtered_ids)} quality-filtered sequence IDs to {filtered_ids_file}...")
    with open(filtered_ids_file, 'w') as f:
        f.write('\n'.join(filtered_ids))
        
    logging.info(f"Saving {len(non_classified_ids)} non-classified sequence IDs to {non_classified_ids_file}...")
    with open(non_classified_ids_file, 'w') as f:
        f.write('\n'.join(non_classified_ids))

    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()
