# src/06_run_recombinhunt/run_recombinhunt.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import yaml
import json
import os
import random
import glob
from collections import defaultdict
from tqdm import tqdm

# --- Path Setup ---
# Add the project's 'src' directory to the path to find local modules like utils and constants.
SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    print("Please ensure RecombinHunt is installed and that 'utils.py' and 'constants.py' are in the 'src/utils' directory.")
    sys.exit(1)

# Add the path to the recombinhunt library source code.
# This assumes a specific directory structure: openrecombinhunt/libs/recombinhunt-cov-7.0.0/src/
PROJECT_ROOT = SRC_PATH.parent
RECOMBINHUNT_LIB_PATH = PROJECT_ROOT / LIBS / "recombinhunt-cov-7.0.0" / SRC
sys.path.append(str(RECOMBINHUNT_LIB_PATH))

try:
    from recombinhunt.core.method import Experiment
    from recombinhunt.core.environment import Environment
    from recombinhunt.validation.utils import compute_X_perc_characterization
    from recombinhunt.validation.case_analysis import CaseAnalysis
except ImportError as e:
    print(f"Error: Could not import RecombinHunt core modules. {e}")
    print("Please ensure RecombinHunt is installed and that the 'recombinhunt' package is accessible.")
    sys.exit(1)

# --- Helper Functions from Notebook ---

def load_lineage_genome_counts(json_filepath: Path) -> dict:
    """Loads the JSON file containing total genome counts for each lineage."""
    logging.info(f"Loading lineage totals from '{json_filepath}'")
    try:
        with open(json_filepath, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        logging.error(f"JSON totals file not found at '{json_filepath}'")
        return {}
    except Exception as e:
        logging.error(f"An unexpected error occurred loading '{json_filepath}': {e}")
        return {}

def get_sample_files(samples_dir: Path) -> list:
    """Finds all 'samples_*.csv' files in the specified directory."""
    if not samples_dir.is_dir():
        logging.error(f"Samples directory '{samples_dir}' not found.")
        return []
    
    pattern = str(samples_dir / "samples_*.csv")
    sample_files = []
    for filepath in glob.glob(pattern):
        filename = os.path.basename(filepath)
        lineage_name = filename[len("samples_"):-len(".csv")]
        sample_files.append((lineage_name, filepath))
    
    logging.info(f"Found {len(sample_files)} sample files in '{samples_dir}'.")
    return sample_files

def calculate_consensus_sequence(all_nuc_changes: list, threshold: float = 0.5) -> list:
    """Calculates the consensus sequence based on nucleotide changes."""
    return compute_X_perc_characterization(strings=all_nuc_changes, threshold=threshold)

# --- Main Workhorse Function from Notebook ---

def run_recombinhunt_for_lineage(
    true_lineage_name: str, 
    sample_filepath: str,
    total_genomes_for_lineage: int,
    case_number_start: int,
    experiment_environment,
    run_mode: str,
    num_samples_to_run: int,
    random_seed: int,
    threshold: float,
    report_file_handle,
    html_report_handle,
    virus_name_for_report: str,
    structured_reports_dir: Path,
    virus_name: str
) -> dict:
    """
    Reads a per-lineage sample file, runs RecombinHunt, and returns detailed statistics.
    """
    results_summary = {
        'true_lineage': true_lineage_name, 'run_mode': run_mode,
        'samples_available': 0, 'samples_run': 0, 'recombinhunt_runs_successful': 0,
        'recombinhunt_runs_failed': 0, 'count_0BP_single_candidate': 0,
        'count_1BP_recombinant': 0, 'count_2BP_recombinant': 0,
        'count_other_recombinant': 0, 'recombinant_genome_ids': [],
        'stat1_exact_match_success': 0, 'stat2_in_candidates_success': 0
    }

    summary_for_filtering = []
    summary_for_extensive_table = []

    try:
        df_lineage = pd.read_csv(sample_filepath, sep='\t', dtype=str)
    except Exception as e:
        logging.error(f"Error reading sample file {sample_filepath}: {e}")
        return results_summary

    """
    # ALREADY HANDLED IN CREATE_SAMPLES.PY WITH NEW UPDATE
    if virus_name == "sars-cov-2" and run_mode.lower() == ALL:
        logging.warning(f"Warning: Running 'all' mode for SARS-CoV-2 may be resource-intensive due to large sample sizes. Reducing to last 6 months.")
        if 'collection_date' in df_lineage.columns:
            logging.info("Filtering samples to those collected in the last 6 months based on 'Collection date'.")
            logging.info(f"Initial sample count: {len(df_lineage)}")
            df_lineage['collection_date'] = pd.to_datetime(df_lineage['collection_date'], errors='coerce')
            six_months_ago = pd.Timestamp.now() - pd.DateOffset(months=6)
            df_lineage = df_lineage[df_lineage['collection_date'] >= six_months_ago]
            logging.info(f"Filtered samples to last 6 months. {len(df_lineage)} samples remain.")
    """

    if virus_name == "sars-cov-2":
        logging.info("get rid of the collection_date column for all rows - SARS-CoV-2 specific handling.")
        if 'collection_date' in df_lineage.columns:
            df_lineage.drop(columns=['collection_date'], inplace=True)

    nuc_to_ids_map = defaultdict(list)
    for _, row in df_lineage.iterrows():
        if pd.notna(row['nuc_changes']) and pd.notna(row['genomeID']):
            nuc_to_ids_map[row['nuc_changes']].append(row['genomeID'])
            
    all_nuc_changes = df_lineage['nuc_changes'].dropna().tolist()
    results_summary['samples_available'] = len(all_nuc_changes)

    if not all_nuc_changes:
        return results_summary
    
    nuc_changes_to_process = []
    log_prefix = f"  Lineage '{true_lineage_name}', Mode '{run_mode}':"

    if run_mode.lower() == ALL:
        nuc_changes_to_process = all_nuc_changes
        logging.info(f"{log_prefix} Processing all {len(all_nuc_changes)} samples.")
    elif run_mode.lower() == CONSENSUS:
        consensus_list = calculate_consensus_sequence(all_nuc_changes, threshold=threshold)
        if consensus_list:
            nuc_changes_to_process = [','.join(consensus_list)]
            logging.info(f"{log_prefix} Processing 1 consensus sequence.")
    elif run_mode.lower() == RANDOM:
        if num_samples_to_run and len(all_nuc_changes) > num_samples_to_run > 0:
            random.seed(random_seed)
            nuc_changes_to_process = random.sample(all_nuc_changes, num_samples_to_run)
            logging.info(f"{log_prefix} Sampled {len(nuc_changes_to_process)} of {len(all_nuc_changes)} available.")
        else:
            nuc_changes_to_process = all_nuc_changes
            logging.info(f"{log_prefix} Requested {num_samples_to_run}, available {len(all_nuc_changes)}. Processing all available.")

    results_summary['samples_run'] = len(nuc_changes_to_process)
    
    if not nuc_changes_to_process:
        logging.warning(f"No valid nucleotide changes to process for lineage '{true_lineage_name}' to run experiments with mode '{run_mode}'.")

    current_case_number = case_number_start

    for nuc_change_data in nuc_changes_to_process:
        try:
            nuc_change_data_list = [s.strip() for s in nuc_change_data.split(',') if s.strip()]
            assert set(nuc_change_data_list).issubset(set(experiment_environment.c2lp_df.index))
            
            experiment = Experiment(environment=experiment_environment) 
            experiment.set_target(nuc_change_data_list)
            result = experiment.run()

            recombinhunt_candidates = []
            if result and result.genome_view and result.genome_view.describe:
                description = result.genome_view.describe()
                recombinhunt_candidates = description.get('designated candidates', [])
                if isinstance(recombinhunt_candidates, str):
                    recombinhunt_candidates = list(map(str.strip, recombinhunt_candidates.split('+')))
            
            num_candidates = len(recombinhunt_candidates)
            if num_candidates == 1:
                results_summary['count_0BP_single_candidate'] += 1
            elif num_candidates >= 2:
                if num_candidates == 2: results_summary['count_1BP_recombinant'] += 1
                elif num_candidates == 3: results_summary['count_2BP_recombinant'] += 1
                else: results_summary['count_other_recombinant'] += 1
                
                ids_for_this_pattern = nuc_to_ids_map.get(nuc_change_data, [])
                results_summary['recombinant_genome_ids'].extend(ids_for_this_pattern)

                if html_report_handle:
                    representative_id = ids_for_this_pattern[0] if ids_for_this_pattern else "Unknown_ID"
                    case_report_folder = structured_reports_dir / f"case_{current_case_number}_{true_lineage_name}_{representative_id}"

                    ca = CaseAnalysis(
                        experiment=experiment,
                        case_name=representative_id,
                        number_of_sequences=total_genomes_for_lineage,
                        case_number=current_case_number,
                        case_group_name=true_lineage_name,
                        lineage_hierarchy=experiment.lh
                    )
                    ca.print_case_details(html_report_handle)
                    ca.save_structured_report(case_report_folder)

                    summary_for_filtering.append({
                        'genomeIDs': ",".join(ids_for_this_pattern),
                        'original_lineage': true_lineage_name,
                        'recombinant_parents': ca.best_candidates_str,
                        'breakpoint_count': f"{num_candidates - 1}BP",
                        'case_report_folder': str(case_report_folder)
                    })

                    summary_for_extensive_table.append(ca.analysis_table_row())

            if len(recombinhunt_candidates) == 1 and recombinhunt_candidates[0] == true_lineage_name:
                results_summary['stat1_exact_match_success'] += 1
            if true_lineage_name in recombinhunt_candidates:
                results_summary['stat2_in_candidates_success'] += 1
            results_summary['recombinhunt_runs_successful'] += 1
        except Exception as e:
            results_summary['recombinhunt_runs_failed'] += 1
            logging.error(f"  Error on one sample in lineage {true_lineage_name}: {e}", exc_info=True)

    if results_summary['recombinhunt_runs_successful'] > 0:
        results_summary['stat1_perc'] = (results_summary['stat1_exact_match_success'] / results_summary['recombinhunt_runs_successful']) * 100
        results_summary['stat2_perc'] = (results_summary['stat2_in_candidates_success'] / results_summary['recombinhunt_runs_successful']) * 100
    else:
        results_summary['stat1_perc'], results_summary['stat2_perc'] = 0.0, 0.0

    results_summary['summary_for_filtering'] = summary_for_filtering
    results_summary['summary_for_extensive_table'] = summary_for_extensive_table
        
    return results_summary

# --- Orchestrator ---

def run_experiments(virus_name: str, config: dict):
    """Main orchestration logic for running RecombinHunt experiments."""
    logging.info(f"--- Starting Step 6: Run RecombinHunt for '{virus_name}' ---")

    # --- 1. Get Paths and Parameters ---
    paths_config = config.get(PATHS, {})
    virus_config = config.get(VIRUSES, {}).get(virus_name, {})
    haplocov_params = virus_config.get(PARAMETERS, {}).get(HAPLOCOV, {})
    recombinhunt_params = virus_config.get(PARAMETERS, {}).get(RECOMBINHUNT, {})

    run_mode = recombinhunt_params.get(RUN_MODE)
    num_samples = recombinhunt_params.get(NUM_SAMPLES_TO_RUN, 50)
    threshold = recombinhunt_params.get(CONSENSUS_THRESHOLD, 0.75)
    
    param_string = f"dist{haplocov_params.get(DIST, 0)}size{haplocov_params.get(SIZE, 0)}"
    
    env_dir = Path(paths_config.get(ENVIRONMENTS)) / virus_name / param_string
    samples_dir = Path(paths_config.get(SAMPLES)) / virus_name / param_string

    if virus_name == "sars-cov-2" and run_mode.lower() == ALL:
        samples_dir = samples_dir / "last_6_months"
    
    results_dir = Path(paths_config.get(RESULTS)) / "recombinhunt_output" / virus_name / param_string
    if virus_name == "sars-cov-2":
        results_dir = results_dir / run_mode
    results_dir.mkdir(parents=True, exist_ok=True)

    structured_reports_dir = results_dir / "structured_reports"
    structured_reports_dir.mkdir(parents=True, exist_ok=True)

    # --- 2. Load RecombinHunt Environment ---
    logging.info(f"Loading RecombinHunt environment from: {env_dir}")
    try:
        env = Environment(str(env_dir))
        env.assert_valid()
    except Exception as e:
        logging.critical(f"Failed to load a valid RecombinHunt environment from {env_dir}: {e}")
        sys.exit(1)

    # --- 3. Setup Report Files ---
    if run_mode == RANDOM:      file_suffix = f"{run_mode}-{num_samples}"
    elif run_mode == CONSENSUS: file_suffix = f"{run_mode}-{threshold}"
    elif run_mode == ALL:       file_suffix = run_mode

    report_filepath = results_dir / f"experiment-output.{file_suffix}.txt"
    html_report_filepath = results_dir / f"experiment-detail.{file_suffix}.html"
    
    # --- 4. Run Experiments ---
    with open(report_filepath, 'w') as report_file, open(html_report_filepath, 'w') as html_file:
        #html_file.write(f"<html><head><title>RecombinHunt Report: {virus_name}</title></head><body>\n")
        #html_file.write(f"<h1>RecombinHunt Detailed Report: {virus_name}</h1>\n<h2>Mode: {run_mode}</h2><hr>\n")
        
        def write_to_report(msg):
            report_file.write(str(msg) + '\n')

        number_of_genomes_per_lineage = load_lineage_genome_counts(samples_dir / "samples_total.json")
        sample_files = get_sample_files(samples_dir)
        all_results_data, all_summary_for_filtering, all_summary_for_extensive_table = [], [], []
        recombinants_by_lineage, case_counter = {}, 1

        for lineage_name, filepath in tqdm(sample_files, desc="Running Experiments"):
            write_to_report(f"\n{'='*10} Processing Lineage: {lineage_name} {'='*10}")
            total_genomes = number_of_genomes_per_lineage.get(lineage_name, 0)
            write_to_report(f"Total genomes available for this lineage: {total_genomes}")

            results = run_recombinhunt_for_lineage(
                true_lineage_name=lineage_name, sample_filepath=filepath,
                total_genomes_for_lineage=total_genomes, case_number_start=case_counter,
                experiment_environment=env, run_mode=run_mode,
                num_samples_to_run=num_samples, random_seed=42, threshold=threshold,
                report_file_handle=report_file, html_report_handle=html_file,
                virus_name_for_report=virus_name,
                structured_reports_dir=structured_reports_dir,
                virus_name=virus_name
            )
            case_counter += 1
            
            if results:
                recombinant_ids = results.get('recombinant_genome_ids', [])
                if recombinant_ids: recombinants_by_lineage[lineage_name] = recombinant_ids
                
                write_to_report(f"  Results for {lineage_name}:")
                write_to_report(f"    Samples Run: {results['samples_run']}")
                write_to_report(f"    Successful Runs: {results['recombinhunt_runs_successful']}")
                write_to_report(f"    Result - No Match: {results['recombinhunt_runs_failed']}")
                write_to_report(f"    Result - 0BP (Single Candidate): {results['count_0BP_single_candidate']}")
                write_to_report(f"    Result - 1BP (Recombinant): {results['count_1BP_recombinant']}")
                write_to_report(f"    Result - 2BP (Recombinant): {results['count_2BP_recombinant']}")
                all_results_data.append(results)
                all_summary_for_filtering.extend(results.get('summary_for_filtering', []))
                all_summary_for_extensive_table.extend(results.get('summary_for_extensive_table', []))

        # --- 5. Write Final Summaries ---
        if all_results_data:
            summary_df = pd.DataFrame(all_results_data)
            summary_df['No Match'] = summary_df['recombinhunt_runs_failed']
            summary_df['True Lineage in Candidates %'] = summary_df.apply(
                lambda row: (row['stat2_in_candidates_success'] / row['recombinhunt_runs_successful'] * 100) if row['recombinhunt_runs_successful'] > 0 else 0,
                axis=1
            )
            summary_df.rename(columns={
                'count_0BP_single_candidate': '0BP (Single Candidate)',
                'count_1BP_recombinant': '1BP', 
                'count_2BP_recombinant': '2BP'
            }, inplace=True)
            
            write_to_report("\n\n" + "="*20 + " Overall Summary Table " + "="*20)
            display_cols = ['true_lineage', 'samples_run', 'True Lineage in Candidates %', 
                            '0BP (Single Candidate)', '1BP', '2BP', 'No Match']
            existing_cols = [col for col in display_cols if col in summary_df.columns]
            write_to_report(summary_df[existing_cols].to_string(index=False))

            write_to_report("\n\n" + "="*20 + " Recombinant Genomes Detected by Lineage " + "="*20)
            if recombinants_by_lineage:
                for lineage, ids in recombinants_by_lineage.items():
                    write_to_report(f"\nLineage: {lineage}")
                    write_to_report("  " + ", ".join(ids))
            else:
                write_to_report("\nNo recombinant genomes were detected in this run.")

        if all_summary_for_filtering:
            filtering_df = pd.DataFrame(all_summary_for_filtering)
            filtering_df.to_csv(results_dir / "recombinant_summary.tsv", index=False, sep='\t')
            write_to_report("\n\n" + "="*20 + " recombinant_summary " + "="*20)
            write_to_report(filtering_df.to_string(index=False))
            logging.info(f"Recombinant summary saved to: {results_dir / 'recombinant_summary.tsv'}")
        else:
            logging.warning("No recombinant summaries were generated. Check if any recombinants were detected.")

        if all_summary_for_extensive_table:
            # Define headers for the extensive table based on the CaseAnalysis output
            extensive_headers = [
                "case_number", "case_name", "n_sequences", "n_changes", "group_name",
                "gt_candidates", "test_status", "initial_region_span", "gap_history",
                "best_candidates", "direction_L1", "rank_L1_L2",
                "bc_breakpoints_target", "gt_breakpoints_target",
                "p_value_vs_L1", "p_value_vs_L2", "alternative_candidates"
            ]
            extensive_df = pd.DataFrame(all_summary_for_extensive_table, columns=extensive_headers)
            extensive_df.to_csv(results_dir / "recombinant_summary_extensive.tsv", sep='\t', index=False)
            write_to_report("\n\n" + "="*20 + " Extensive Summary Table " + "="*20)
            write_to_report(extensive_df.to_string(index=False))
            logging.info(f"Saved extensive summary table to: {results_dir / 'recombinant_summary_extensive.tsv'}")
        else:
            logging.warning("No extensive summary table was generated. Check if any cases were analyzed.")
            
        html_file.write("</body></html>\n")
        logging.info(f"Text report saved to: {report_filepath}")
        logging.info(f"HTML report saved to: {html_report_filepath}")

def main():
    """Main function to run the RecombinHunt experiment step."""
    parser = argparse.ArgumentParser(description="Run RecombinHunt experiments on prepared samples.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    try:
        config = yaml.safe_load(Path(args.config).read_text())
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_06_run_recombinhunt")

    try:
        run_experiments(args.virus, config)
        logging.info("RecombinHunt experiment step finished successfully with given config and run mode.")

        if args.virus.lower() == "sars-cov-2" and config.get(VIRUSES).get(args.virus).get(PARAMETERS).get(RECOMBINHUNT).get(RUN_MODE) != CONSENSUS:
            logging.info("For sars-cov-2, it's recommended to run the CONSENSUS mode next for a comprehensive analysis.")
            logging.info("Now running the CONSENSUS mode as a follow-up step.")
            config.get(VIRUSES).get(args.virus).get(PARAMETERS).get(RECOMBINHUNT)[RUN_MODE] = CONSENSUS
            run_experiments(args.virus, config)
            logging.info("Follow-up CONSENSUS mode run finished successfully.")

    except Exception as e:
        logging.critical(f"RecombinHunt experiment step failed: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
