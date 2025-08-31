# src/02_preprocessing/ncbi/prep_metadata_ncbi.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import re
from typing import Tuple, List, Set
import yaml

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    sys.exit(1)

# --- Hardcoded Column Mappings & Constants ---
# This map defines the final output columns and their source.
COLUMN_MAPPING = {
    "Accession": "Virus name",
    "Isolate Collection date": "Collection date",
    "Release date": "Submission date"
}
# The primary ID column in the raw NCBI metadata
SOURCE_ID_COL = "Accession"
# Special characters to look for in the ID to drop the row
SPECIAL_CHARS_PATTERN = re.compile(r'[&/\\:]')

# --- Placeholder Formatting Functions ---

def _format_single_date(date_str: str) -> str:
    """
    Helper function to format a single date string according to specific rules.
    - YYYY-MM-DD -> YYYY-MM-DD
    - YYYY-MM    -> YYYY-MM-01
    - YYYY       -> YYYY-01-01
    Returns pd.NA for invalid or unparseable formats.
    """
    if pd.isna(date_str):
        return pd.NA
    
    date_str = str(date_str).strip()
    date_str = date_str.split('T')[0] # Remove time part if present
    
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

def format_date_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Applies standardized date formatting to the relevant columns and logs format statistics.
    """
    logging.info("Formatting date columns...")
    
    date_cols_to_process = ['Collection date', 'Submission date']
    
    for col_name in date_cols_to_process:
        if col_name in df.columns:
            logging.info(f"  - Analyzing and formatting '{col_name}'")
            
            # --- New Logging Logic ---
            # Work with a string representation of the series for pattern matching
            series_str = df[col_name].astype(str)
            
            # Count initial NAs by checking for pandas NA and common null-like strings
            na_mask = df[col_name].isna() | series_str.str.lower().isin(['nan', 'na', 'nat', '<na>', ''])
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
            df[col_name] = df[col_name].apply(_format_single_date)
            
    return df

def format_location_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Combines 'Geographic Region' and 'Geographic Location' into a single 'Location' column.
    The final format is 'region / country'. Handles missing values.
    """
    logging.info("Formatting 'Location' column by combining 'Geographic Region' and 'Geographic Location'...")
    
    region_col = 'Geographic Region'
    country_col = 'Geographic Location'

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

def process_ncbi_metadata(df: pd.DataFrame, filter_rules: list, source_lineage_col: str) -> Tuple[pd.DataFrame, Set[str], Set[str]]:
    """
    Applies a full preprocessing workflow to the raw NCBI metadata.
    """
    df_processed = df.copy()
    
    # --- 1. Apply Quality & Special Character Filters ---
    logging.info("Step 1: Applying quality filters from config and checking for special characters in ID...")
    filtered_ids = set()

    # Apply filters from config file
    for rule in filter_rules:
        col, op, val = rule.get(COLUMN), rule.get(OPERATOR), rule.get(VALUE)
        if col not in df_processed.columns:
            logging.warning(f"Filter rule column '{col}' not found. Skipping filter.")
            continue
        
        initial_rows = len(df_processed)
        ids_to_drop_this_step = set()
        
        if op == "notna":
            keep_mask = df_processed[col].notna()
        elif op == ">=":
            numeric_col = pd.to_numeric(df_processed[col], errors='coerce')
            keep_mask = numeric_col >= val
        else: # Add other operators as needed
            logging.warning(f"Filter operator '{op}' not implemented. Skipping.")
            continue
            
        dropped_mask = ~keep_mask
        ids_to_drop_this_step = set(df_processed.loc[dropped_mask, SOURCE_ID_COL].dropna())
        filtered_ids.update(ids_to_drop_this_step)
        df_processed = df_processed[keep_mask]
        logging.info(f"  Filter '{col} {op} {val}': Removed {initial_rows - len(df_processed)} rows.")
        
    # Special character filter on the ID column
    initial_rows = len(df_processed)
    special_char_mask = df_processed[SOURCE_ID_COL].astype(str).str.contains(SPECIAL_CHARS_PATTERN, regex=True, na=True)
    ids_to_drop_special_chars = set(df_processed.loc[special_char_mask, SOURCE_ID_COL].dropna())
    filtered_ids.update(ids_to_drop_special_chars)
    df_processed = df_processed[~special_char_mask]
    logging.info(f"  Filter 'Special Chars in ID': Removed {initial_rows - len(df_processed)} rows.")

    # --- 2. Handle Pango Lineage based on config ---
    logging.info("Step 2: Processing lineage column...")
    non_classified_ids = set()
    
    if source_lineage_col and source_lineage_col.upper() != 'NONE':
        logging.info(f"Using source column '{source_lineage_col}' for lineage information.")
        if source_lineage_col in df_processed.columns:
            is_na_mask = df_processed[source_lineage_col].isna()
            is_unclassifiable_mask = df_processed[source_lineage_col].str.lower().isin(['unclassified', 'unknown', 'na', 'n/a', 'not applicable'])
            combined_drop_mask = is_na_mask | is_unclassifiable_mask

            ids_to_drop_no_lineage = set(df_processed.loc[combined_drop_mask, SOURCE_ID_COL].dropna())
            non_classified_ids.update(ids_to_drop_no_lineage)
            df_processed = df_processed[~combined_drop_mask]
            logging.info(f"  Identified and removed {len(ids_to_drop_no_lineage)} non-classified rows (missing '{source_lineage_col}').")
            # Add the mapping for the Pango lineage column for the final step
            COLUMN_MAPPING[source_lineage_col] = 'Pango lineage'
        else:
            logging.error(f"Source lineage column '{source_lineage_col}' not found in data. Exiting.")
            sys.exit(1)
    else:
        logging.info("No source lineage column specified ('NONE'). Assigning default lineage 'A.1' to all rows.")
        df_processed['Pango lineage'] = 'A.1'

    # --- 3. Final Formatting on the Cleaned DataFrame ---
    logging.info("Step 3: Performing final formatting on date and location columns...")
    
    # Rename columns to standard format
    df_processed = df_processed.rename(columns=COLUMN_MAPPING)
    
    # Format specific columns
    df_processed = format_date_columns(df_processed)
    df_processed = format_location_column(df_processed)
    
    # Keep only the columns defined in the mapping's values
    expected_final_columns = ["Virus name", "Collection date", "Submission date", "Location", "Pango lineage"]
    final_columns_to_keep = [col for col in expected_final_columns if col in df_processed.columns]
    df_final = df_processed[final_columns_to_keep]
    
    logging.info(f"Preprocessing complete. Final dataset has {len(df_final)} rows.")
    
    return df_final, sorted(list(filtered_ids)), sorted(list(non_classified_ids))


def main():
    """Main function to run the NCBI metadata preprocessing."""
    parser = argparse.ArgumentParser(description="Preprocess raw NCBI metadata based on rules in config.yaml.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    config = yaml.safe_load(Path(args.config).read_text())
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_02.1_prep_metadata_ncbi")

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
        logging.error(f"Raw metadata file not found. Please run the '01_fetch_data' script first.")
        sys.exit(1)

    filter_rules = virus_config.get(PARAMETERS, {}).get(METADATA_PROCESSING, {}).get(FILTERS, [])
    source_lineage_col = virus_config.get(PARAMETERS, {}).get(METADATA_PROCESSING, {}).get(SOURCE_LINEAGE_COLUMN, "NONE")
    
    df_processed, filtered_ids, non_classified_ids = process_ncbi_metadata(df_raw, filter_rules, source_lineage_col)
    
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
