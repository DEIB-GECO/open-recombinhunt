# src/05_prepare_recombinhunt/create_environment.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import yaml
from collections import Counter
from typing import Dict, Tuple

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    sys.exit(1)

# --- Helper Functions ---

def is_recombinant(lineage_name: str) -> bool:
    """Identifies top-level recombinant lineages (e.g., 'XA', 'XB')."""
    if not isinstance(lineage_name, str):
        return False
    lineage_name = lineage_name.upper()
    return lineage_name.startswith('X') and '.' not in lineage_name

def extract_pos(variation_string: str) -> int:
    """Extracts the genomic position from a mutation string."""
    try:
        return int(variation_string.split('_', 1)[0])
    except (ValueError, IndexError):
        return pd.NA

# --- Core Calculation Functions ---

def calculate_change_probability_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the baseline probability of each mutation across the entire dataset.
    
    Args:
        df (pd.DataFrame): The cleaned DataFrame containing a 'mutations' column.

    Returns:
        pd.DataFrame: A DataFrame with 'variation' as the index and 
                        columns for 'probability' and 'pos'.
    """
    logging.info("Calculating change probabilities...")
    total_sequences = len(df)
    if total_sequences == 0:
        logging.warning("Input DataFrame is empty. Cannot calculate change probabilities.")
        return pd.DataFrame(columns=['probability', 'pos'])
    logging.info(f"  Total sequences for calculation: {total_sequences:,}")

    all_mutations_list = df['mutations'].dropna().str.cat(sep=',').split(',')
    mutation_counts = Counter(all_mutations_list)
    assert '' not in mutation_counts, "Data Quality Error: Empty string found as a mutation. Check for trailing commas or empty parts in the 'mutations' column of the input file."
    logging.info(f"  Found {len(mutation_counts):,} unique mutations across all sequences.")

    if not mutation_counts:
        logging.warning("No mutations found in the input data. Returning empty DataFrame.")
        return pd.DataFrame(columns=['probability', 'pos'])
    
    prob_df = pd.DataFrame.from_dict(mutation_counts, orient='index', columns=['count'])
    prob_df['probability'] = prob_df['count'] / total_sequences
    prob_df['pos'] = prob_df.index.to_series().apply(extract_pos)
    
    prob_df = prob_df[['probability', 'pos']].sort_values(by='pos')
    prob_df.index.name = None
    logging.info(f"  Finished calculating change probabilities. Shape of result: {prob_df.shape}")
    return prob_df

def get_lineage_and_variation_counts(df: pd.DataFrame) -> Tuple[Dict, Dict]:
    """
    Calculates genome and variation counts per lineage.
    
    Calculates the necessary count dictionaries from the main DataFrame to be used
    as input for the lineage characterization function.

    Args:
        df (pd.DataFrame): The cleaned DataFrame containing 'pango_lineage' and 'mutations'.

    Returns:
        A tuple containing two dictionaries:
        - genome_counter_by_lineage (dict): {lineage: count_of_genomes}
        - variation_counter_by_lineage (dict): {lineage: Counter({mutation: count})}
    """
    logging.info("Calculating genome and variation counts per lineage...")
    genome_counter_by_lineage = df['pangoLin'].value_counts().to_dict()
    logging.info(f"  Found {len(genome_counter_by_lineage):,} unique lineages in the dataset.")

    # Calculate the counts of each mutation within each lineage
    variation_counter_by_lineage: Dict[str, Counter] = {
        lineage: Counter() for lineage in genome_counter_by_lineage.keys()
    }

    for index, row in df.iterrows():
        current_lineage = row['pangoLin']
        mutations = str(row['mutations'])
        
        # Split by comma, strip whitespace, and filter out empty strings
        variations_in_row = [var.strip() for var in mutations.split(',') if var.strip()]
        
        if variations_in_row:
            variation_counter_by_lineage[current_lineage].update(variations_in_row)

    logging.info("  Finished calculating genome and variation counts.")
    return genome_counter_by_lineage, variation_counter_by_lineage

def _calculate_lineage_frequencies_df(
    genome_counter_by_lineage: dict,
    variation_counter_by_lineage: dict,
    min_genome_count: int = 10,
    threshold: float = 0.75
) -> pd.DataFrame:
    """
    Internal helper to calculate a DataFrame of mutation frequencies per lineage.
    """
    logging.info("\n--- Running common logic: calculating lineage mutation frequencies ---")
    
    selected_lineages = {
        lineage: count for lineage, count in genome_counter_by_lineage.items() 
        if count >= min_genome_count
    }
    logging.info(f"Selected {len(selected_lineages)} lineages with at least {min_genome_count} genomes.")

    if not selected_lineages:
        logging.info("Warning: No lineages met the min_genome_count requirement.")
        return pd.DataFrame()

    lc = {lineage:  
        {variation: count/genome_counter_by_lineage[lineage] 
         for (variation, count) in variation_counts.items()}
         for (lineage, variation_counts) in variation_counter_by_lineage.items() 
         if lineage in selected_lineages}

    lc = {lineage: 
        {variation: freq 
         for (variation,freq) in variation_frequencies.items() if freq >= threshold}
         for (lineage, variation_frequencies) in lc.items()}
    
    distinct_variations = sorted({c for (_,variations) in lc.items() for c in variations.keys()})
    sorted_lineages = sorted(lc.keys())

    if not lc: # Handle case where lc is empty after filtering
        logging.info("Warning: No lineages/variations met the criteria. Returning an empty DataFrame.")
        return pd.DataFrame()

    table = []
    for l in sorted_lineages:
        all_variations_frequency_in_l = list(map(lambda c: lc[l].get(c, 0.0), distinct_variations))
        table.append(all_variations_frequency_in_l)

    lc_df_ = pd.DataFrame(table, columns=distinct_variations, index=sorted_lineages).T
    del table

    return lc_df_

def calculate_lc_df(
    genome_counter_by_lineage: dict,
    variation_counter_by_lineage: dict,
    min_genome_count: int = 10,
    char_threshold: float = 0.75
) -> pd.DataFrame:
    """
    Creates the lc_df (lineage characterization) DataFrame. The cells are
    boolean, indicating if a mutation is characteristic of a lineage.
    """
    logging.info("\nCalculating 'lc_df.parquet'...")
    
    # Get the base frequency data
    frequency_df = _calculate_lineage_frequencies_df(
        genome_counter_by_lineage, variation_counter_by_lineage, min_genome_count, threshold=char_threshold
    )
    if frequency_df.empty:
        return pd.DataFrame()
    
    # Convert frequencies to boolean: True if frequency >= char_threshold, else False
    lc_df = frequency_df >= char_threshold
    
    # Remove rows that are not characteristic of any lineage
    initial_mutation_count = len(lc_df)
    lc_df = lc_df[lc_df.any(axis=1)]
    logging.info(f"Removed {initial_mutation_count - len(lc_df)} mutations that were not characteristic of any lineage.")

    lc_df['lc_pos'] = lc_df.index.to_series().apply(extract_pos)
    
    # Final formatting
    cols = sorted([col for col in lc_df.columns if col != 'lc_pos']) + ['lc_pos']
    lc_df = lc_df[cols]
    lc_df.sort_values(by='lc_pos', inplace=True)
    lc_df.index.name = None

    logging.info("lc_df calculation complete.")
    return lc_df

def calculate_change2lineage_probability_df(
    genome_counter_by_lineage: dict,
    variation_counter_by_lineage: dict,
    min_genome_count: int = 10
) -> pd.DataFrame:
    """
    Creates the change2lineage_probability DataFrame. The cells are numeric,
    representing the frequency of a mutation within a lineage.
    """
    logging.info("\nCalculating 'change2lineage_probability.parquet'...")
    
    # Get the base frequency data
    prob_df = _calculate_lineage_frequencies_df(
        genome_counter_by_lineage, variation_counter_by_lineage, min_genome_count, threshold=0.0
    )
    if prob_df.empty:
        return pd.DataFrame()

    prob_df['pos'] = prob_df.index.to_series().apply(extract_pos)
    
    # Final formatting
    cols = sorted([col for col in prob_df.columns if col != 'pos']) + ['pos']
    prob_df = prob_df[cols]
    prob_df.sort_values(by='pos', inplace=True)
    prob_df.index.name = None

    logging.info("change2lineage_probability calculation complete.")
    return prob_df

def calculate_lc_quality_df(genome_counter_by_lineage: dict) -> pd.DataFrame:
    """
    Creates the lc_quality_df, which contains the number of sequences for each lineage.

    Args:
        genome_counter_by_lineage (dict): A dictionary mapping lineage to its genome count.

    Returns:
        pd.DataFrame: A DataFrame with lineages as the index and 'n_sequences' as the column.
    """
    logging.info("\nCalculating 'lc_quality_df.parquet'...")
    if not genome_counter_by_lineage:
        print("Warning: genome_counter_by_lineage is empty. Returning an empty DataFrame.")
        return pd.DataFrame()

    # Convert the dictionary to a DataFrame
    # orient='index' makes the dictionary keys the DataFrame index
    lc_quality_df = pd.DataFrame.from_dict(
        {k: genome_counter_by_lineage[k] for k in sorted(genome_counter_by_lineage) if not is_recombinant(k)}, orient='index', columns=["num"]
    ).T
    
    logging.info("lc_quality_df calculation complete.")
    return lc_quality_df

def main():
    """Main function to orchestrate environment creation."""
    parser = argparse.ArgumentParser(description="Create environment files for RecombinHunt.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    config = yaml.safe_load(Path(args.config).read_text())
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_05.1_create_environment")

    # --- 1. Determine Paths and Parameters ---
    virus_name = args.virus
    virus_config = config.get(VIRUSES, {}).get(virus_name, {})
    paths_config = config.get(PATHS, {})
    
    haplocov_params = virus_config.get(PARAMETERS, {}).get(HAPLOCOV, {})
    recombinhunt_params = virus_config.get(PARAMETERS, {}).get(RECOMBINHUNT, {})
    
    param_string = f"dist{haplocov_params.get(DIST, 0)}size{haplocov_params.get(SIZE, 0)}"
    
    # --- 2. Determine Input File (Conditional Logic) ---
    results_dir = Path(paths_config.get(RESULTS))
    if virus_name == 'sars-cov-2':
        logging.info("Virus is 'sars-cov-2'. Using Nextstrain reformatted as input.")
        input_file = results_dir / "nextstrain_output" / virus_name / "nextstrain_reformatted.tsv"
    else:
        logging.info(f"Virus is '{virus_name}'. Using HaploCoV reformatted output as input.")
        input_file = results_dir / "haplocov_output" / virus_name / param_string / "haplocov_reformatted.tsv"

    # Determine Output Directory
    env_dir = Path(paths_config.get(ENVIRONMENTS)) / virus_name / param_string
    env_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Environment files will be saved to: {env_dir}")

    # --- 3. Load and Prepare Data ---
    logging.info(f"Loading data from: {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file}. Please run previous pipeline steps.")
        sys.exit(1)

    # --- NEW: Clean the DataFrame before any calculations ---
    logging.info(f"Initial row count: {len(df):,}")
    initial_rows = len(df)
    # Drop rows where 'mutations' column is NaN or an empty string
    df.dropna(subset=['mutations'], inplace=True)
    df = df[df['mutations'] != '']
    rows_dropped = initial_rows - len(df)
    if rows_dropped > 0:
        logging.info(f"Dropped {rows_dropped:,} rows with missing or empty 'mutations' column.")
    logging.info(f"Row count after cleaning: {len(df):,}")
    # --- End of New Section ---

    # --- 4. Run Calculations ---
    min_genomes = recombinhunt_params.get(MIN_GENOME_COUNT, 10)
    lc_thresh = recombinhunt_params.get(LC_THRESHOLD, 0.75)

    # Calculate change probabilities
    change_prob_df = calculate_change_probability_df(df)

    # Calculate base statistics
    genome_counts, variation_counts = get_lineage_and_variation_counts(df)
    
    # Calculate final environment DataFrames
    lc_df = calculate_lc_df(
        genome_counts, 
        variation_counts,
        min_genomes,
        lc_thresh
    )

    c2lp_df = calculate_change2lineage_probability_df(
        genome_counts,
        variation_counts,
        min_genomes
    )

    lc_quality_df = calculate_lc_quality_df(genome_counts)

    # --- 5. Save Outputs (Conditional Logic for SARS-CoV-2) ---
    logging.info("Saving environment files...")
    change_prob_df.to_parquet(env_dir / "change_probability.parquet")
    lc_quality_df.to_parquet(env_dir / "lc_quality_df.parquet")

    if virus_name == 'sars-cov-2':
        logging.info("Virus is 'sars-cov-2'. Saving two versions of lineage characterisations.")
        # Save the version with recombinants
        c2lp_df.to_parquet(env_dir / 'change2lineage_probability_with_X.parquet')
        lc_df.to_parquet(env_dir / "lc_df_with_X.parquet")

        # Calculate the versions without recombinants
        genome_counts_non_recombinants, variation_counts_non_recombinants = get_lineage_and_variation_counts(
            df[~df["pangoLin"].apply(is_recombinant)]
        )

        c2lp_df_non_recombinants = calculate_change2lineage_probability_df(
            genome_counts_non_recombinants,
            variation_counts_non_recombinants,
            min_genomes
        )

        lc_df_non_recombinants = calculate_lc_df(
            genome_counts_non_recombinants,
            variation_counts_non_recombinants,
            min_genomes,
            lc_thresh
        )

        c2lp_df_non_recombinants.to_parquet(env_dir / 'change2lineage_probability.parquet')
        lc_df_non_recombinants.to_parquet(env_dir / "lc_df.parquet")


    else:
        # For all other viruses, save the single file
        c2lp_df.to_parquet(env_dir / 'change2lineage_probability.parquet')
        lc_df.to_parquet(env_dir / "lc_df.parquet")
        
    logging.info("Environment creation complete.")

if __name__ == "__main__":
    main()
