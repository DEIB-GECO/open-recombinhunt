# src/02_preprocessing/ftp/prep_metadata_ftp.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import yaml
import re
from typing import Tuple, List

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    sys.exit(1)

# --- Formatting Helper Functions ---

def format_location_column(series: pd.Series) -> pd.Series:
    """
    Formats the location column by keeping everything before the first ':'
    if a ':' is present.
    """
    logging.info("  - Formatting 'Location' column...")
    # Ensure the series is treated as string, handle NAs by filling them
    # so they are not processed by .str methods.
    return series.astype(str).str.split(':', n=1).str[0].str.strip()

def _format_single_date(date_str: str) -> str:
    """
    Helper to format a single date string by padding it.
    - YYYY-MM-DD -> YYYY-MM-DD
    - YYYY-MM    -> YYYY-MM-01
    - YYYY       -> YYYY-01-01
    """
    if pd.isna(date_str):
        return pd.NA
    
    date_str = str(date_str).strip().split('T')[0]
    
    if re.fullmatch(r'\d{4}-\d{2}-\d{2}', date_str): return date_str
    elif re.fullmatch(r'\d{4}-\d{2}', date_str): return f"{date_str}-01"
    elif re.fullmatch(r'\d{4}', date_str): return f"{date_str}-01-01"
    else: return pd.NA

def format_date_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Applies standardized date formatting to the relevant columns."""
    logging.info("  - Formatting date columns...")
    if 'Collection date' in df.columns:
        df['Collection date'] = df['Collection date'].apply(_format_single_date)
    if 'Submission date' in df.columns:
        df['Submission date'] = df['Submission date'].apply(_format_single_date)
    return df

# --- Main Processing Function ---

def process_ftp_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """
    Applies the full filtering and formatting workflow for the FTP metadata.
    """
    logging.info("Starting metadata processing...")
    initial_rows = len(df)
    
    # --- 1. Apply Filters ---
    logging.info("Step 1: Applying filtering rules...")
    
    # To handle potential errors with .str accessor on mixed-type columns,
    # we convert columns to string type for filtering where needed.
    df['Genotype'] = df['Genotype'].astype(str)
    df['Segment'] = df['Segment'].astype(str)
    df['Geo_Location'] = df['Geo_Location'].astype(str)

    filtered_df = df[
        (df["Species"] == "Alphainfluenzavirus influenzae") &
        (df["Genotype"].str.contains("H5N1|h5n1", na=False)) &
        (df["Segment"].str.contains("HA|segment 4|Segment 4|4", na=False)) &
        (
            (df["Geo_Location"].str.contains("USA|Canada|Mexico", na=False) & ~df["Geo_Location"].str.contains(":", na=False)) |
            (df["Geo_Location"].str.contains("USA:|Canada:|Mexico:", na=False) & df["Geo_Location"].str.contains(":", na=False))
        ) &
        (df["Length"] > 1700) &
        (df["Collection_Date"].notna()) &
        (df["Release_Date"].notna())
    ].copy()
    
    logging.info(f"  Filtering complete. Kept {len(filtered_df)} rows out of {initial_rows}.")

    if filtered_df.empty:
        logging.warning("No rows remained after filtering. Output will be empty.")
        return filtered_df

    # --- 2. Rename Columns ---
    logging.info("Step 2: Renaming columns to standard format...")
    column_mapping = {
        "#Accession": "Virus name",
        "Collection_Date": "Collection date",
        "Release_Date": "Submission date",
        "Geo_Location": "Location"
    }
    df_renamed = filtered_df.rename(columns=column_mapping)
    
    # --- 3. Create Pango Lineage Column ---
    logging.info("Step 3: Creating 'Pango lineage' column and setting to 'A.1'...")
    df_renamed['Pango lineage'] = 'A.1'

    # --- 4. Format Columns ---
    logging.info("Step 4: Formatting Location and Date columns...")
    df_renamed['Location'] = format_location_column(df_renamed['Location'])
    df_formatted = format_date_columns(df_renamed)

    # --- 5. Select and Order Final Columns ---
    final_columns = [
        "Virus name", "Collection date", "Submission date",
        "Location", "Pango lineage"
    ]
    # Ensure we only select columns that exist to prevent KeyErrors
    existing_final_columns = [col for col in final_columns if col in df_formatted.columns]
    df_final = df_formatted[existing_final_columns]

    # Include Continent in the Location column.
    # Add "North America / " to the beginning of the Location values.
    df_final['Location'] = "North America / " + df_final['Location']

    logging.info(f"Processing complete. Final DataFrame has {len(df_final)} rows.")
    return df_final


def main():
    """Main function to run the FTP metadata preprocessing."""
    parser = argparse.ArgumentParser(description="Preprocess raw FTP metadata for Influenza.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process (e.g., 'influenza').")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    config = yaml.safe_load(Path(args.config).read_text())
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_02.1_prep_metadata_ftp")

    paths_config = config.get(PATHS, {})
    
    raw_dir = Path(paths_config.get(RAW_DATA, 'data/raw')) / args.virus
    processed_dir = Path(paths_config.get(PROCESSED_DATA, 'data/processed')) / args.virus
    processed_dir.mkdir(parents=True, exist_ok=True)
    
    input_file = raw_dir / "raw_metadata.csv"
    output_file = processed_dir / "metadata.tsv"
    accession_ids_file = processed_dir / f"{args.virus}-accession-ids.txt"

    logging.info(f"Loading raw metadata from {input_file}...")
    try:
        # Use comment='#' to handle the '#Accession' column name if it causes issues
        df_raw = pd.read_csv(input_file, low_memory=False)
        # If the first column is now 'Accession' instead of '#Accession', fix it
        if 'Accession' in df_raw.columns and '#Accession' not in df_raw.columns:
            df_raw.rename(columns={'Accession': '#Accession'}, inplace=True)

    except FileNotFoundError:
        logging.error(f"Raw metadata file not found. Please run the '01_fetch_data' script first.")
        sys.exit(1)

    df_processed = process_ftp_metadata(df_raw)
    
    logging.info(f"Saving processed metadata ({len(df_processed)} rows) to {output_file}...")
    df_processed.to_csv(output_file, sep='\t', index=False, na_rep='NA')

    # Extract and save the accession IDs for the next step
    if 'Virus name' in df_processed.columns:
        accession_ids = df_processed['Virus name'].dropna().unique().tolist()
        logging.info(f"Saving {len(accession_ids)} unique accession IDs to {accession_ids_file}...")
        with open(accession_ids_file, 'w') as f:
            f.write('\n'.join(accession_ids))
    else:
        logging.warning("Could not find 'Virus name' column in final DataFrame. Accession ID file will be empty.")
        accession_ids_file.touch() # Create empty file

    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()
