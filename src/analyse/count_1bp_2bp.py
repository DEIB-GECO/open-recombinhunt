import pandas as pd
from pathlib import Path

# Configuration - change these parameters for different viruses
virus = 'sars-cov-2'
dist = 0
size = 0

# Construct the file path based on virus, dist, and size parameters
paramset = f"dist{dist}size{size}"
file_path = Path(f"results/recombinhunt_output/{virus}/{paramset}/all/recombinant_summary.tsv")

print(f"Reading file: {file_path}")

try:
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Check if breakpoint_count column exists
    if 'breakpoint_count' not in df.columns:
        print("Error: 'breakpoint_count' column not found in the file")
        print(f"Available columns: {list(df.columns)}")
        exit(1)
    
    # Count 1BP and 2BP
    num_rows = len(df)
    count_1bp = (df['breakpoint_count'] == '1BP').sum()
    count_2bp = (df['breakpoint_count'] == '2BP').sum()
    total_recombinants = count_1bp + count_2bp
    
    print(f"\nBreakpoint Count Results for {virus} (dist={dist}, size={size}):")
    print(f"Total count: {num_rows}")
    print(f"1BP count: {count_1bp}")
    print(f"2BP count: {count_2bp}")
    print(f"Total recombinants: {total_recombinants}")
    print(f"Total sequences: {len(df)}")
    
    if total_recombinants > 0:
        print(f"1BP percentage: {(count_1bp/total_recombinants)*100:.2f}%")
        print(f"2BP percentage: {(count_2bp/total_recombinants)*100:.2f}%")
    
except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
    print("Please check the virus name, dist, and size parameters")
except Exception as e:
    print(f"Error reading file: {e}")

