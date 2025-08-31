# src/04_postprocessing/format_haplocov_variations.py

import argparse
import pandas as pd
from pathlib import Path
import sys
import logging
import yaml
from tqdm import tqdm
from typing import List

# Add the project's 'src' directory to the Python path.
SRC_PATH = Path(__file__).resolve().parent.parent
sys.path.append(str(SRC_PATH))

try:
    from utils.utils import setup_logging
    from utils.constants import *
except ImportError as e:
    print(f"Error: Could not import a required module. {e}")
    sys.exit(1)

# --- Caching Dictionary ---
# This global dictionary will store already computed transformations to speed up processing.
FIXES = {}

# --- Formatting Functions for Each Mutation Type ---

def format_deletion(mutation: str) -> str:
    """
    Transforms a HaploCoV deletion format to the RecombinHunt format.
    - Single deletion: '3622_G|.' -> '3622'
    - Multiple deletion: '3622_GATACC|......' -> '3622_3627'
    """
    try:
        parts = mutation.split('_')
        start_pos = int(parts[0])
        
        deleted_bases = parts[1].split('|')[0]
        num_deleted = len(deleted_bases)
        
        if num_deleted == 1:
            return str(start_pos)
        else:
            end_pos = start_pos + num_deleted - 1
            return f"{start_pos}_{end_pos}"
    except (ValueError, IndexError):
        logging.warning(f"Could not parse supposed deletion '{mutation}'. Returning original.")
        return mutation


def format_insertion(mutation: str) -> str:
    """
    Transforms a HaploCoV insertion format to the RecombinHunt format.
    Example: '7124_....|TTTA' -> '7124_.|TTTA'
    """
    try:
        left_full, right = mutation.split('|', 1)
        pos, left_dots = left_full.split('_', 1)
        return f"{pos}_.|{right}"
    except (ValueError, IndexError):
        logging.warning(f"Could not parse supposed insertion '{mutation}' with standard format. Returning original.")
        return mutation


def format_compound_mutation(mutation: str) -> List[str]:
    """
    Deconstructs a complex compound mutation string into a list of simple mutations,
    now correctly handling multi-base substitutions.
    """
    try:
        left_full, right_bases = mutation.split('|', 1)
        start_pos_str, left_bases = left_full.split('_', 1)
        start_pos = int(start_pos_str)
        
        if len(left_bases) != len(right_bases):
            logging.warning(f"Compound mutation '{mutation}' has mismatched lengths. Returning original.")
            return [mutation]

        results = []
        i = 0
        while i < len(left_bases):
            current_pos = start_pos + i
            left_char = left_bases[i]
            right_char = right_bases[i]

            # Case 1: Substitution block (base -> base)
            if left_char != '.' and right_char != '.':
                sub_start_pos = current_pos
                sub_end_index = i
                sub_left = ""
                sub_right = ""
                # Find the end of this substitution block
                while sub_end_index < len(left_bases) and left_bases[sub_end_index] != '.' and right_bases[sub_end_index] != '.':
                    sub_left += left_bases[sub_end_index]
                    sub_right += right_bases[sub_end_index]
                    sub_end_index += 1
                
                results.append(f"{sub_start_pos}_{sub_left}|{sub_right}")
                i = sub_end_index
            
            # Case 2: Deletion block (base -> dot)
            elif left_char != '.' and right_char == '.':
                del_start_pos = current_pos
                del_end_index = i
                # Find the end of this deletion block
                while del_end_index < len(left_bases) and left_bases[del_end_index] != '.' and right_bases[del_end_index] == '.':
                    del_end_index += 1
                
                num_deleted = del_end_index - i
                if num_deleted == 1:
                    results.append(str(del_start_pos))
                else:
                    del_end_pos = del_start_pos + num_deleted - 1
                    results.append(f"{del_start_pos}_{del_end_pos}")
                i = del_end_index

            # Case 3: Insertion block (dot -> base)
            elif left_char == '.' and right_char != '.':
                ins_pos = current_pos
                ins_end_index = i
                inserted_sequence = ""
                # Find the end of this insertion block
                while ins_end_index < len(left_bases) and left_bases[ins_end_index] == '.' and right_bases[ins_end_index] != '.':
                    inserted_sequence += right_bases[ins_end_index]
                    ins_end_index += 1
                
                results.append(f"{ins_pos}_.|{inserted_sequence}")
                i = ins_end_index
            
            else: # Should not happen (dot -> dot)
                i += 1
        
        return results

    except Exception as e:
        logging.warning(f"Error parsing compound mutation '{mutation}': {e}. Returning original.")
        return [mutation]


def format_listv_cell(cell_value: str) -> str:
    """
    Processes a single cell from the 'listV' column. It splits the comma-separated
    mutations, identifies their type, formats them, and joins them back.
    """
    if pd.isna(cell_value):
        return ""

    final_reformatted_mutations = []
    mutations = [mut.strip() for mut in str(cell_value).split(',') if mut.strip()]

    for mut in mutations:
        if mut in FIXES:
            final_reformatted_mutations.extend(FIXES[mut])
            continue

        parsed_mutations = []
        
        if '.' not in mut:
            parsed_mutations = [mut]
        else:
            try:
                left_full, right = mut.split('|', 1)
                
                is_pure_deletion = (right.strip() == '.' * len(right.strip()))
                
                if '_' in left_full:
                    pos, left_bases = left_full.split('_', 1)
                    is_pure_insertion = (left_bases.strip() == '.' * len(left_bases.strip()))
                else:
                    is_pure_insertion = False

                if is_pure_deletion and not is_pure_insertion:
                    parsed_mutations = [format_deletion(mut)]
                elif is_pure_insertion and not is_pure_deletion:
                    parsed_mutations = [format_insertion(mut)]
                else:
                    # It's a compound mutation if it contains dots but isn't a pure deletion/insertion
                    parsed_mutations = format_compound_mutation(mut)

            except ValueError:
                logging.warning(f"Could not parse mutation '{mut}' with '.'. Returning original.")
                parsed_mutations = [mut]
        
        FIXES[mut] = parsed_mutations
        final_reformatted_mutations.extend(parsed_mutations)

    return ",".join(final_reformatted_mutations)


def main():
    """Main function to run the HaploCoV output post-processing."""
    parser = argparse.ArgumentParser(description="Post-process HaploCoV output to reformat mutation strings.")
    parser.add_argument("--virus", required=True, help="The name of the virus to process.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to the main YAML configuration file.")
    args = parser.parse_args()

    try:
        config = yaml.safe_load(Path(args.config).read_text())
    except FileNotFoundError:
        print(f"CRITICAL ERROR: Config file not found at '{args.config}'", file=sys.stderr)
        sys.exit(1)
    
    log_dir = Path(config.get(PATHS, {}).get(LOGS, 'logs'))
    setup_logging(log_dir=log_dir, log_name_prefix=f"{args.virus}_04_postprocess_haplocov")

    # --- Get Paths ---
    paths_config = config.get(PATHS, {})
    results_base_dir = Path(paths_config.get(RESULTS, 'results'))
    
    virus_config = config.get(VIRUSES).get(args.virus)
    haplocov_params = virus_config.get(PARAMETERS).get(HAPLOCOV)
    dist = haplocov_params.get(DIST)
    size = haplocov_params.get(SIZE)
    param_string = f"dist{dist}size{size}"
    
    input_dir = results_base_dir / "haplocov_output" / args.virus / param_string
    input_file = input_dir / "haplocov_assigned.tsv"
    output_file = input_dir / "haplocov_reformatted.tsv"

    # --- Load and Process Data ---
    logging.info(f"Loading HaploCoV output from: {input_file}")
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
    except FileNotFoundError:
        logging.error(f"Input file not found: {input_file}. Please run the '03_run_haplocov' script first.")
        sys.exit(1)

    if 'listV' not in df.columns:
        logging.error(f"Required column 'listV' not found in {input_file}. Cannot proceed.")
        sys.exit(1)

    logging.info("Reformatting mutations in the 'listV' column...")
    
    tqdm.pandas(desc="Formatting mutations")
    df['mutations'] = df['listV'].progress_apply(format_listv_cell)

    df.drop(columns=['listV'], inplace=True)
    
    logging.info(f"Saving reformatted data to: {output_file}")
    df.to_csv(output_file, sep='\t', index=False, na_rep='NA')

    logging.info("Script finished successfully.")


if __name__ == "__main__":
    try:
        from tqdm import tqdm
    except ImportError:
        print("Error: tqdm library not found. Please install it: pip install tqdm")
        sys.exit(1)
    main()
    
    #mock = "6959_CTT...|GGAGTT, 5122_C|A,8701_G|T, 3622_G|., 10415_AACCTGAAACCGGGA|...............,5776_A|T, 10368_AGGAAT|.....C,354_C|T"
    #print(format_listv_cell(mock))
