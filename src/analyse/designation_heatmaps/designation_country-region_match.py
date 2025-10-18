import argparse
import pandas as pd
import sys
import os 
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_grouped_percentages(df: pd.DataFrame, group_column_name: str, pango_column_name: str = 'pangoLin') -> pd.DataFrame:
    """
    Calculates the percentage of designations within each group (e.g., continent or country).
    Returns a DataFrame with group, pangoLin, count, total_in_group, and percentage.
    """
    analysis_df = df[[group_column_name, pango_column_name]].dropna().copy()

    if analysis_df.empty:
        print(f"Warning: No valid data for {group_column_name}/{pango_column_name} analysis after dropping rows with missing values.", file=sys.stderr)
        return pd.DataFrame(columns=[group_column_name, pango_column_name, 'count', 'total_in_group', 'percentage'])

    designation_counts_df = analysis_df.groupby([group_column_name, pango_column_name]).size().reset_index(name='count')
    if designation_counts_df.empty:
        return pd.DataFrame(columns=[group_column_name, pango_column_name, 'count', 'total_in_group', 'percentage'])

    total_entries_per_group_df = designation_counts_df.groupby(group_column_name)['count'].sum().reset_index(name='total_in_group')
    merged_df = pd.merge(designation_counts_df, total_entries_per_group_df, on=group_column_name)

    if merged_df.empty:
        return pd.DataFrame(columns=[group_column_name, pango_column_name, 'count', 'total_in_group', 'percentage'])

    merged_df['percentage'] = (merged_df['count'] / merged_df['total_in_group']) * 100
    
    merged_df_sorted = merged_df.sort_values(
        by=[group_column_name, 'percentage', pango_column_name],
        ascending=[True, False, True]
    )
    return merged_df_sorted

def save_text_summary_to_file(continent_data: pd.DataFrame, country_data: pd.DataFrame, output_filepath: str):
    """
    Formats and saves the textual summary of designation percentages to a file.
    """
    print(f"Saving text summary to: {output_filepath}")
    try:
        with open(output_filepath, 'w') as f_out:
            f_out.write("--- Analyzing Designations by Continent (Percentages) ---\n")
            if continent_data.empty:
                f_out.write("No data for continent summary.\n")
            else:
                current_group_val = None
                for _, row in continent_data.iterrows():
                    if row['continent'] != current_group_val:
                        current_group_val = row['continent']
                        f_out.write(f"\nContinent: {current_group_val} (Total observed designations in this group: {row['total_in_group']})\n")
                    f_out.write(f"  Designation: {row['pangoLin']}, Percentage: {row['percentage']:.2f}%\n")

            f_out.write("\n\n--- Analyzing Designations by Country (Percentages) ---\n")
            if country_data.empty:
                f_out.write("No data for country summary.\n")
            else:
                current_group_val = None
                for _, row in country_data.iterrows():
                    if row['country'] != current_group_val:
                        current_group_val = row['country']
                        f_out.write(f"\nCountry: {current_group_val} (Total observed designations in this group: {row['total_in_group']})\n")
                    f_out.write(f"  Designation: {row['pangoLin']}, Percentage: {row['percentage']:.2f}%\n")
        print(f"Text summary successfully saved.")
    except Exception as e:
        print(f"Error saving text summary to {output_filepath}: {e}", file=sys.stderr)


def visualize_continent_summary_and_save(continent_data_df: pd.DataFrame, style: str, output_dir: str):
    """
    Visualizes the continent summary data and saves the plot to a file.
    """
    if continent_data_df.empty:
        print("No continent data to visualize.")
        return

    if not {'continent', 'pangoLin', 'percentage'}.issubset(continent_data_df.columns):
        print("Error: DataFrame for visualization is missing required columns ('continent', 'pangoLin', 'percentage').", file=sys.stderr)
        return

    print(f"\n--- Generating and saving '{style}' plot for Continent Summary ---")
    plt.style.use('seaborn-v0_8-whitegrid')
    fig = None # Initialize fig variable

    if style == "stacked":
        output_filename = "stackbar.png"
        try:
            pivot_df_stacked = continent_data_df.pivot_table(
                index='continent', columns='pangoLin', values='percentage', fill_value=0
            )
            if pivot_df_stacked.empty:
                print("Pivoted data for stacked bar chart is empty.")
                return

            fig, ax = plt.subplots(figsize=(15, 8)) # Create figure and axes
            pivot_df_stacked.plot(kind='bar', stacked=True, ax=ax, colormap='Spectral')
            ax.set_title('Pango Lineage Percentage Distribution within Continents')
            ax.set_ylabel('Percentage (%)')
            ax.set_xlabel('Continent')
            plt.xticks(rotation=45, ha='right')
            ax.legend(title='Pango Lineage', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize='small')
            fig.tight_layout(rect=[0, 0, 0.85, 1])
            
            filepath = os.path.join(output_dir, output_filename)
            plt.savefig(filepath)
            print(f"Saved stacked bar plot to: {filepath}")
        except Exception as e:
            print(f"Error generating or saving stacked bar plot: {e}", file=sys.stderr)
        finally:
            if fig:
                plt.close(fig) # Close the specific figure

    elif style == "heatmap":
        output_filename = "heatmap.png"
        try:
            pivot_df_heatmap = continent_data_df.pivot_table(
                index='pangoLin', columns='continent', values='percentage', fill_value=0
            )
            if pivot_df_heatmap.empty:
                print("Pivoted data for heatmap is empty.")
                return

            fig_height = max(8, len(pivot_df_heatmap.index) * 0.4)
            fig_width = max(10, len(pivot_df_heatmap.columns) * 1.2) # Increased width factor a bit
            
            fig = plt.figure(figsize=(fig_width, fig_height)) # Assign to fig
            sns.heatmap(pivot_df_heatmap, annot=True, fmt=".1f", cmap="Blues", linewidths=.5, cbar_kws={'label': 'Percentage (%)'})
            plt.title('Pango Lineage Percentage Heatmap across Continents')
            plt.ylabel('Pango Lineage')
            plt.xlabel('Continent')
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout()

            filepath = os.path.join(output_dir, output_filename)
            plt.savefig(filepath)
            print(f"Saved heatmap to: {filepath}")
        except Exception as e:
            print(f"Error generating or saving heatmap: {e}", file=sys.stderr)
        finally:
            if fig:
                plt.close(fig) # Close the specific figure
            
    else:
        print(f"Unknown visualization style: '{style}'. Please choose 'stacked' or 'heatmap'.")


def main():
    parser = argparse.ArgumentParser(
        description="Analyzes designation occurrences (as percentages) by geographic location and saves outputs."
    )
    parser.add_argument("input_file", metavar="INPUT_HAPLOCOV_FILE", type=str, help="Path to the HaploCoV output file.")
    parser.add_argument("--sep", type=str, default='\t', help="Delimiter used in the input file. Default is tab.")
    parser.add_argument(
        "--visualize_style_continent", type=str, default="heatmap",
        choices=["none", "stacked", "heatmap", "both"],
        help="Style of visualization for continent data: 'none', 'stacked', 'heatmap', or 'both'. Default is 'both'."
    )
    args = parser.parse_args()

    try:
        # Determine output directory from input file path
        abs_input_filepath = os.path.abspath(args.input_file)
        output_directory = os.path.dirname(abs_input_filepath)

        df = pd.read_csv(abs_input_filepath, sep=args.sep, low_memory=False)
    except FileNotFoundError:
        print(f"Error: Input file not found at '{args.input_file}'", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: Input file '{args.input_file}' is empty.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading or parsing file '{args.input_file}': {e}", file=sys.stderr)
        sys.exit(1)

    required_columns = ['continent', 'country', 'pangoLin']
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        print(f"Error: Input file '{args.input_file}' is missing required columns: {', '.join(missing_cols)}", file=sys.stderr)
        sys.exit(1)

    # --- Calculate Percentages ---
    continent_summary_df = calculate_grouped_percentages(df, 'continent')
    country_summary_df = calculate_grouped_percentages(df, 'country')

    # --- Save Text Summary ---
    text_output_filepath = os.path.join(output_directory, "distribution-of-designations-over-regions.txt")
    save_text_summary_to_file(continent_summary_df, country_summary_df, text_output_filepath)

    # --- Visualize and Save Continent Data ---
    if args.visualize_style_continent != "none":
        if args.visualize_style_continent in ["stacked", "both"]:
            visualize_continent_summary_and_save(continent_summary_df, "stacked", output_directory)
        if args.visualize_style_continent in ["heatmap", "both"]:
            visualize_continent_summary_and_save(continent_summary_df, "heatmap", output_directory)

    print("\n\nAnalysis and output saving complete.")

if __name__ == "__main__":
    main()