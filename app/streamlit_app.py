import datetime as dt
import json
import sys
from pathlib import Path
import base64
import yaml
import os
import re
import pandas as pd
import streamlit as st
from streamlit_option_menu import option_menu
import streamlit.components.v1 as components
import plotly.express as px
import plotly.graph_objects as go
from geopy.geocoders import Nominatim
from agstyler import draw_grid, PINLEFT, PRECISION_TWO
from about_virus import describe

st.set_page_config(
    page_title="OpenRecombinHunt",
    page_icon="🧬",
    layout="wide"
)

# --- Path Setup for Imports ---
# Add the project root to the Python path to allow imports from 'src'
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(PROJECT_ROOT))

mapping = {
    "sars-cov-2": "SARS-CoV-2",
    "yellow-fever": "Yellow Fever",
    "zika": "Zika",
    "rsv-a": "RSV-A",
    "rsv-b": "RSV-B",
    "monkeypox": "Monkeypox",
    "influenza": "Influenza (H5N1)"
}

try:
    from src.utils.constants import *
except ImportError as e:
    print(f"CRITICAL ERROR: Could not import from 'src/utils/constants.py'. {e}")
    sys.exit(1)

# --- Configuration & Path Setup ---
CONFIG_PATH = Path("config/config.yaml")

try:
    with open(CONFIG_PATH, 'r') as f:
        config = yaml.safe_load(f)
    RESULTS_DIR_BASE = Path(config.get(PATHS).get(RESULTS))
    REFERENCE_LENGTHS = {virus: config.get(VIRUSES).get(virus).get(REFERENCE).get(LENGTH) for virus in mapping.keys()}
except FileNotFoundError:
    st.error(f"Configuration file not found at {CONFIG_PATH}. Please ensure the file exists.")
    st.stop()
except Exception as e:
    st.error(f"Error loading or parsing configuration file: {e}")
    st.stop()

def visualize(virus):
    return mapping.get(virus.lower())

def loading_animation(virus_name):
    # write text on top of loading animation
    loader_html = f"""
        <style>
        .loader-overlay {{
            position: fixed;
            top: 0; left: 0;
            width: 100%; height: 100%;
            background: rgba(255, 255, 255, 0.95);
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            z-index: 9999;
            font-family: 'Segoe UI', sans-serif;
        }}
        .loader {{
            border: 8px solid #f3f3f3;
            border-top: 8px solid #3498db;
            border-radius: 50%;
            width: 70px;
            height: 70px;
            animation: spin 1s linear infinite;
            margin-bottom: 20px;
        }}
        @keyframes spin {{
            0% {{ transform: rotate(0deg); }}
            100% {{ transform: rotate(360deg); }}
        }}
        .loading-text {{
            font-size: 20px;
            color: #333;
            text-align: center;
            animation: fadein 2s infinite alternate;
        }}
        @keyframes fadein {{
            from {{ opacity: 0.4; }}
            to {{ opacity: 1; }}
        }}
        </style>
        <div class="loader-overlay">
            <div class="loader"></div>
            <div class="loading-text">
                The {virus_name} data is heavy.<br>
                Pulling the sequences from the archives...
            </div>
        </div>
    """
    return loader_html

@st.cache_data
def discover_viruses():
    """discover available viruses by scanning recombinhunt_output/ directory"""
    recombinhunt_output_path = RESULTS_DIR_BASE / RECOMBINHUNT_OUTPUT

    if not recombinhunt_output_path.exists():
        st.warning(f"No recombinhunt output found at {recombinhunt_output_path}.")
        return []

    virus_dirs = [d for d in recombinhunt_output_path.iterdir() if d.is_dir()]
    if not virus_dirs:
        st.warning(f"No virus directories found in {recombinhunt_output_path}.")
        return []

    # return only the names, not paths
    return [d.name for d in virus_dirs]

@st.cache_data
def load_master_data(virus):
    """load and merge master data for the specified virus"""

    # load param set
    try:
        if virus not in ["sars-cov-2"]:
            params = config.get(VIRUSES).get(virus).get(PARAMETERS).get(HAPLOCOV)
            dist = params.get(DIST)
            size = params.get(SIZE)
        elif virus in ["sars-cov-2"]:
            dist, size = 0, 0
        paramset = f"dist{dist}size{size}"
    except Exception as e:
        st.error(f"Error loading parameters for {virus}: {e}")
        return None

    # load source data
    try:
        if virus.lower() == "sars-cov-2":
            # Get analysis window from config to construct the correct filename
            try:
                virus_config = config.get(VIRUSES).get(virus)
                analysis_window_months = virus_config.get("analysis_window_months", 6)
            except Exception:
                analysis_window_months = 6
            
            source_file = RESULTS_DIR_BASE / NEXTSTRAIN_OUTPUT / virus / f"nextstrain_reformatted_last_{analysis_window_months}_months.tsv"
            columns = ["genomeID", "Collection date", "Submission date", "Location", "pangoLin"]
        else:
            source_file = RESULTS_DIR_BASE / HAPLOCOV_OUTPUT / virus / paramset / "haplocov_reformatted.tsv"
            columns = ["genomeID", "collectionD", "continent", "country", "pangoLin"]

        if source_file.exists():
            source_df = pd.read_csv(source_file, sep="\t", usecols=columns)
        else:
            st.warning(f"Source file not found: {source_file}")
    except Exception as e:
        st.error(f"Error loading master data for {virus}: {e}")
        return None
    
    # load recombinant summary
    try:
        recombinant_summary_file_base = RESULTS_DIR_BASE / RECOMBINHUNT_OUTPUT / virus / paramset
    
        if virus == "sars-cov-2": recombinant_summary_file = recombinant_summary_file_base / ALL / "recombinant_summary.tsv"
        else:                     recombinant_summary_file = recombinant_summary_file_base / "recombinant_summary.tsv" 

        if recombinant_summary_file.exists():
            recombinant_summary_df = pd.read_csv(recombinant_summary_file, sep="\t")
            recombinant_summary_df.rename(columns={"genomeIDs": "genomeID"}, inplace=True)
        else:
            st.warning(f"Recombinant summary file not found: {recombinant_summary_file}")
    except Exception as e:
        st.error(f"Error loading recombinant summary for {virus}: {e}")
        return None

    # merge dataframes
    try:
        merged_df = pd.merge(source_df, recombinant_summary_df, on="genomeID", how="left")
    except Exception as e:
        st.error(f"Error merging dataframes for {virus}: {e}")
        return None

    # add is_recombinant column to merged_df
    # if breakpoint_count.notnull()
    merged_df["is_recombinant"] = merged_df["breakpoint_count"].notnull()

    # rename some columns of merged_df
    mapping = {
        "collectionD": "collection_date",
        "Collection date": "collection_date",
        "Submission date": "submission_date",
    }
    merged_df.rename(columns=mapping, inplace=True)

    return merged_df

@st.cache_data
def load_consensus_data(virus):
    if virus == "sars-cov-2":
        dist, size = 0, 0
        paramset = f"dist{dist}size{size}"
        recombinant_summary_file_base = RESULTS_DIR_BASE / RECOMBINHUNT_OUTPUT / virus / paramset / CONSENSUS / "recombinant_summary.tsv"

    if recombinant_summary_file_base.exists():
        df = pd.read_csv(recombinant_summary_file_base, sep="\t")
    else:
        st.warning(f"Recombinant summary file not found: {recombinant_summary_file_base}")
        return None

    return df

@st.cache_data
def load_complete_data_stats(virus):
    stats = {}
    if virus == "sars-cov-2":
        source_file = RESULTS_DIR_BASE / NEXTSTRAIN_OUTPUT / virus / "nextstrain_reformatted.tsv"
        columns = ["genomeID", "Collection date", "Submission date", "Location", "pangoLin"]
        df = pd.read_csv(source_file, sep="\t", usecols=columns)
        
        stats["total_records"] = len(df)

        df["Collection date"] = pd.to_datetime(df["Collection date"], errors="coerce")
        min_date = df["Collection date"].min()
        max_date = df["Collection date"].max()
        stats["min_collection_date"] = min_date.strftime("%Y-%m-%d") if pd.notnull(min_date) else "N/A"
        stats["max_collection_date"] = max_date.strftime("%Y-%m-%d") if pd.notnull(max_date) else "N/A"

        stats["unique_countries"] = df["Location"].apply(lambda x: x.split("/")[1].strip() if isinstance(x, str) else x).nunique()
        stats["country_distribution"] = df["Location"].apply(lambda x: x.split("/")[1].strip() if isinstance(x, str) else x).value_counts().reset_index()
        stats["unique_lineages"] = df["pangoLin"].nunique()
        stats["lineage_distribution"] = df["pangoLin"].value_counts().reset_index()

        return stats
    else:
        pass

def apply_time_filter(df, virus):
    """
    applies time based filtering to the dataframe
    user can select between three options:
        -no time filtering (returns all data)
        -filter by specific date selections (choose start and end dates: end date defaults to today)
        -filter by the number of latest X sequences (user inputs an integer number, we sort the dataframe by collection date and return the latest X rows)
    """
    # select filter type
    filter_type = st.pills(
        "Choose filter type:",
        options = ["No filtering", "Filter by Date Range", "Filter by Latest Sequences"],
        selection_mode = "single"
    )

    # Get the download date and analysis window from config for the virus
    try:
        virus_config = config.get(VIRUSES).get(virus)
        download_date_str = virus_config.get("download_date")
        analysis_window_months = virus_config.get("analysis_window_months", 6)  # Default to 6 months if not specified
        
        if download_date_str:
            download_date = pd.to_datetime(download_date_str).date()
        else:
            st.error(f"No download date found for {virus}")
            download_date = dt.date.today()  # fallback to today if no download_date found
    except Exception as e:
        st.error(f"Error getting download date for {virus}: {e}")
        download_date = dt.date.today()  # fallback to today if any error occurs
        analysis_window_months = 6
    
    # Calculate the actual date range in the data for reference
    try:
        collection_dates = pd.to_datetime(df["collection_date"], errors='coerce').dropna()
        if not collection_dates.empty:
            min_data_date = collection_dates.min().date()
            max_data_date = collection_dates.max().date()
            st.info(f"Data contains collection dates from {min_data_date} to {max_data_date}")
        else:
            st.warning("No valid collection dates found in the data")
    except Exception as e:
        st.warning(f"Could not analyze collection date range: {e}")

    # select filter value
    filter_value = None
    if filter_type == "Filter by Date Range":
        col1, col2 = st.columns(2)
        with col1:
            if virus == "sars-cov-2":
                # Use pandas DateOffset for accurate time window calculation (same as format_covid_variations.py)
                download_date_pd = pd.to_datetime(download_date)
                window_start = download_date_pd - pd.DateOffset(months=analysis_window_months)
                min_allowed_date = window_start.date()
                start_date = st.date_input("Start Date", value=min_allowed_date, min_value=min_allowed_date, max_value=download_date)
                st.info(f"For SARS-CoV-2, the start date cannot be earlier than {min_allowed_date}. Only the last {analysis_window_months} months of data is available.")
            else:
                start_date = st.date_input("Start Date", value=pd.to_datetime("2020-01-01"), max_value=download_date)
        with col2:
            end_date = st.date_input("End Date", value=download_date, min_value=start_date, max_value=download_date)
            st.info(f"Data was downloaded on: {download_date}")
        filter_value = (start_date, end_date)
    elif filter_type == "Filter by Latest Sequences":
        filter_value = st.number_input("Number of Latest Sequences", min_value=1, value=100)

    # apply filtering to the df
    filtered_df = df
    if filter_type == "Filter by Date Range" and filter_value:
        start_date, end_date = filter_value
        # convert string df['collection_date'] to datetime object
        collection_date = pd.to_datetime(df["collection_date"]).dt.date

        filtered_df = df[(collection_date >= start_date) & (collection_date <= end_date)]
    elif filter_type == "Filter by Latest Sequences" and filter_value:
        filtered_df = df.sort_values("collection_date", ascending=False).head(filter_value)

    return filtered_df

def create_key_metrics(df):
    """
    creates key metrics as cards for the summary dashboard
        -# Total Sequences
        -# Recombinant Events
        -Recombination Rate
        -Top Recombinant Lineage
        -Most Common Parents
    """
    st.title("Key Metrics")

    total_sequences = len(df)
    total_recombinants = df["is_recombinant"].sum()
    num_1BP = df[df["breakpoint_count"] == "1BP"].shape[0]
    num_2BP = df[df["breakpoint_count"] == "2BP"].shape[0]
    recombination_rate = (total_recombinants / total_sequences * 100) if total_sequences > 0 else 0

    # get top recombinant lineage
    # and top recombinant lineage count
    top_recombinant_value_counts = df[df["is_recombinant"]]["pangoLin"].value_counts()
    top_recombinant_lineage = top_recombinant_value_counts.idxmax() if not top_recombinant_value_counts.empty else "N/A"
    top_recombinant_count = top_recombinant_value_counts.max() if not top_recombinant_value_counts.empty else 0

    # get most common parents
    # and most common parents count
    most_common_parents_value_counts = df[df["is_recombinant"]]["recombinant_parents"].value_counts()
    most_common_parents = most_common_parents_value_counts.idxmax() if not most_common_parents_value_counts.empty else "N/A"
    most_common_parents_count = most_common_parents_value_counts.max() if not most_common_parents_value_counts.empty else 0

    # Calculate unique patterns (unique recombinant parents across all lineages)
    unique_patterns = df[df["is_recombinant"]]["recombinant_parents"].nunique() if not df[df["is_recombinant"]].empty else 0

    # Create compact metric cards using custom layout
    col1, col2 = st.columns(2)
    col3, col4, col5 = st.columns(3)
    col6, col7, col8 = st.columns(3)
    
    with col1:
        st.markdown(f"""
        <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
            <div style="font-size: 0.8rem; color: #666;">Total Sequences</div>
            <div style="font-size: 1.2rem; font-weight: bold;">{total_sequences:,}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown(f"""
        <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
            <div style="font-size: 0.8rem; color: #666;">Recombinant Sequences</div>
            <div style="font-size: 1.2rem; font-weight: bold;">{total_recombinants:,}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown(f"""
        <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
            <div style="font-size: 0.8rem; color: #666;">1BP</div>
            <div style="font-size: 1.2rem; font-weight: bold;">{num_1BP:,}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col4:
        st.markdown(f"""
        <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
            <div style="font-size: 0.8rem; color: #666;">2BP</div>
            <div style="font-size: 1.2rem; font-weight: bold;">{num_2BP:,}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col5:
        st.markdown(f"""
        <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
            <div style="font-size: 0.8rem; color: #666;">Unique Patterns</div>
            <div style="font-size: 1.2rem; font-weight: bold;">{unique_patterns:,}</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col6:
        st.markdown(f"""
        <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
            <div style="font-size: 0.8rem; color: #666;">Top Recombinant Lineage</div>
            <div style="font-size: 1.0rem; font-weight: bold;">{top_recombinant_lineage} ({top_recombinant_count})</div>
        </div>
        """, unsafe_allow_html=True)
    
    with col7:
        st.markdown(f"""
        <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
            <div style="font-size: 0.8rem; color: #666;">Most Common Patterns</div>
            <div style="font-size: 1.0rem; font-weight: bold;">{most_common_parents} ({most_common_parents_count})</div>
        </div>
        """, unsafe_allow_html=True)

def create_summary_tables(df):
    "create summary and hotspots tables"
    st.title("Summary Tables")

    with st.spinner("Creating lineage breakdown table..."):
        st.subheader("Lineage Breakdown")
        with st.expander("", expanded=True):
            # create a table
            # group by each lineage (pangoLin)
            # count 1BP in breakpoint_count as 1BP Count
            # rate of 1BP Count / total sequences
            # count 2BP in breakpoint_count as 2BP Count
            # rate of 2BP Count / total sequences
            # rest is No Recombination
            # total sequences in each lineage
            lineage_breakdown = df.groupby("pangoLin").agg(
                BP1_Count=("breakpoint_count", lambda x: (x == "1BP").sum()),
                BP1_Rate=("breakpoint_count", lambda x: (x == "1BP").mean() * 100 if len(x) > 0 else 0.00),
                BP2_Count=("breakpoint_count", lambda x: (x == "2BP").sum()),
                BP2_Rate=("breakpoint_count", lambda x: (x == "2BP").mean() * 100 if len(x) > 0 else 0.00),
                No_Recombination=("breakpoint_count", lambda x: ((x != "1BP") & (x != "2BP")).sum()),
                Total_Sequences=("breakpoint_count", "size")
            ).reset_index()

            lineage_breakdown.rename(columns={
                "pangoLin": "Lineage",
                "BP1_Count": "1BP Count",
                "BP1_Rate": "1BP Rate",
                "BP2_Count": "2BP Count",
                "BP2_Rate": "2BP Rate",
                "No_Recombination": "No Recombination",
                "Total_Sequences": "Total Sequences"
            }, inplace=True)

            lineage_breakdown.set_index("Lineage", inplace=True)

            lineage_breakdown.sort_values(by="1BP Rate", ascending=False, inplace=True)

            st.dataframe(
                lineage_breakdown.style.format({
                    "1BP Rate": "{:.2f}%",
                    "2BP Rate": "{:.2f}%"
                }),
                use_container_width=True
            )

    with st.spinner("Creating recombination hotspots table..."):
        st.subheader("Breakdown of Detected Recombinations")
        with st.expander("", expanded=True):
            # group by recombinant_parents
            # display frequency
            recombination_hotspots = df.groupby("recombinant_parents").size().reset_index(name="Frequency")
            recombination_hotspots.rename(columns={"recombinant_parents": "Recombinant Parents"}, inplace=True)
            recombination_hotspots.set_index("Recombinant Parents", inplace=True)
            recombination_hotspots.sort_values(by="Frequency", ascending=False, inplace=True)
            st.write(recombination_hotspots)

def create_temporal_plot(df, virus):
    if "collection_date" not in df.columns:
        st.warning("Collection date information is not available.")
        return

    with st.spinner("Generating temporal distribution plot..."):
        df = df.dropna(subset=["collection_date"])
        df["collection_date"] = pd.to_datetime(df["collection_date"])

        date_range = df["collection_date"].max() - df["collection_date"].min()

        if date_range > pd.Timedelta(days=180):
            freq = "M" 
            freq_label = "Month"  
        else:
            freq = "W" 
            freq_label = "Week" 

        df["year-month"] = df["collection_date"].dt.to_period(freq)

        monthly_data = df.groupby("year-month").agg(
            {
                "is_recombinant": ["count", "sum"]
            }
        )

        monthly_data.columns = [
            "total_sequences", "recombinations"
        ]

        monthly_data = monthly_data.reset_index()
        monthly_data["year-month"] = monthly_data["year-month"].astype(str)

        fig = go.Figure()

        fig.add_trace(
            go.Scatter(
                x=monthly_data["year-month"],
                y=monthly_data["total_sequences"],
                fill="tonexty",
                mode="none",
                name=f"Total Sequences per {freq_label}",
                fillcolor="rgba(74, 144, 226, 0.15)",
                line=dict(width=0)
            )
        )

        fig.add_trace(
            go.Scatter(
                x=monthly_data["year-month"],
                y=monthly_data["recombinations"],
                mode="lines+markers",
                name="Recombination Events",
                line=dict(color="#4A90E2",
                                width=3,
                                shape="spline",
                                smoothing=1.3),
                marker=dict(size=6, color="#4A90E2")
            )
        )

        fig.update_layout(
            xaxis_title="Time",
            yaxis_title="Number of Sequences",
            hovermode="x unified",
            height=500,
            showlegend=True
        )

        if virus.lower() == "sars-cov-2":
            fig.update_layout(yaxis_type="log")
            fig.update_yaxes(title="Number of Sequences (log scale)")

        st.plotly_chart(fig, width="stretch")

def create_geographic_map(df, virus):
    # only keep is_recombinant
    df = df[df["is_recombinant"] == True]
    if df.empty:
        st.warning("No recombinant sequences found.")
        return

    with st.spinner("Generating geographic distribution map..."):
        if virus == "sars-cov-2":
            df['country'] = df["Location"].apply(lambda x: x.split("/")[1].strip() if isinstance(x, str) else x)

        df['country'] = df['country'].str.strip()
        geo_data = df["country"].value_counts().reset_index()
        geo_data.columns = ["country", "count"]

        lat_lon_df = pd.read_csv("app/country.csv")

        geo_data = geo_data.merge(lat_lon_df, on="country", how="left")
        
        missing_countries = geo_data[geo_data["latitude"].isna()]["country"].tolist()

        if missing_countries:
            print("Missing countries (will geocode):", missing_countries)

            for country in missing_countries:
                try:
                    loc = Nominatim(user_agent="GetLoc")
                    getLoc = loc.geocode(country)

                    geo_data.loc[geo_data["country"] == country, "latitude"] = getLoc.latitude
                    geo_data.loc[geo_data["country"] == country, "longitude"] = getLoc.longitude
                except Exception as e:
                    print(f"Could not geocode {country}: {e}")

        geo_data = geo_data.dropna(subset=["latitude", "longitude"])

        fig = px.scatter_mapbox(
            geo_data,
            lat="latitude",
            lon="longitude",
            size="count",
            hover_name="country",
            color="count",
            color_continuous_scale=px.colors.sequential.Blues,
            size_max=40,
            zoom=1
        )

        fig.update_layout(mapbox_style="carto-positron")
        fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})

        a, _ = st.columns([3,2])
        with a:
            st.plotly_chart(fig, use_container_width=True)

def create_distribution_plots(df, virus):
    "creates temporal and locational distributions"
    st.title("Distribution Plots")

    st.subheader("Temporal Distribution")
    create_temporal_plot(df, virus)

    st.subheader("Locational Distribution")
    create_geographic_map(df, virus)

def apply_user_filter(df, virus):
    """Apply user-defined filters to the DataFrame."""
    # filter types
    a, b, c = st.columns(3)

    with a:
        lineage_filter = st.selectbox(
            'Select Lineage:',
            ["All"] + sorted(df["pangoLin"].dropna().unique().tolist()) if "pangoLin" in df.columns else ["NA"]
        )

    with b:
        breakpoint_filter = st.selectbox(
            "Breakpoint Count:",
            ["All"] + sorted(df["breakpoint_count"].dropna().unique().tolist()) if "breakpoint_count" in df.columns else ["NA"]
        )

    with c:
        if virus == "sars-cov-2": 
            df["continent"] = df["Location"].apply(lambda x: x.split("/")[0].strip() if isinstance(x, str) else x)
            df["country"] = df["Location"].apply(lambda x: x.split("/")[1].strip() if isinstance(x, str) else x)
        else: df["continent"] = df["continent"].apply(lambda x: x.strip() if isinstance(x, str) else x)
        location_filter = st.selectbox(
            "Select Continent:",
            ["All"] + sorted(df["continent"].dropna().unique().tolist()) if "continent" in df.columns else ["NA"]
        )

    # apply filters
    if lineage_filter not in ["All", "NA"]:
        df = df[df["pangoLin"] == lineage_filter]

    if breakpoint_filter not in ["All", "NA"]:
        df = df[df["breakpoint_count"] == breakpoint_filter]

    if location_filter not in ["All", "NA"]:
        df = df[df["continent"] == location_filter]

    return df

@st.cache_data
def load_report_data(path):
    """Load report data from a specified path."""
    
    report = {}

    with st.spinner("Loading detailed report for the genome..."):
        # load summary.json
        summary_path = os.path.join(path, "summary.json")
        if os.path.exists(summary_path):
            with open(summary_path, "r") as f:
                report["summary"] = json.load(f)

        # load all region tables
        # files named region_*_table.csv: 
        # * in [1, 2] if 1BP
        # * in [1, 2, 3] if 2BP
        region_files = [f for f in os.listdir(path) if f.startswith("region_") and f.endswith("_table.csv")]
        for region_file in region_files:
            region_path = os.path.join(path, region_file)
            if os.path.exists(region_path):
                with open(region_path, "r") as f:
                    report[region_file] = pd.read_csv(f)

        # load plots (in json format)
        # plot_per_region.json
        # plot_whole_genome.json
        plot_files = [f for f in os.listdir(path) if f.startswith("plot_") and f.endswith(".json")]
        for plot_file in plot_files:
            plot_path = os.path.join(path, plot_file)
            if os.path.exists(plot_path):
                with open(plot_path, "r") as f:
                    report[plot_file] = json.load(f)

        # load target_mutations.txt
        target_mutations_path = os.path.join(path, "target_mutations.txt")
        if os.path.exists(target_mutations_path):
            with open(target_mutations_path, "r") as f:
                report["target_mutations"] = f.read().splitlines()

    return report

def format_region_table(df: pd.DataFrame, region: str) -> pd.DataFrame:
    # Drop index column if it exists (Streamlit shows it by default otherwise)
    df = df.reset_index(drop=True)

    # Rename columns according to mapping
    rename_map = {
        "Unnamed: 0": "Alternative Lineage",
        "num_seq": "#Seq",
        "t_ch_MAX": "max-L1" if region == "1" else "max-L2" if region == "2" else "max-L3",
        "max_CL": "LR",
        "aic": "AIC",
        "PV": "p-value",
        "PV_OK": "C1",
        "t_ch_MAX_OK": "C2",
        "phyl_OK": "C3",
    }
    df = df.rename(columns=rename_map)

    # Drop unnecessary column(s)
    cols_to_drop = [c for c in df.columns if "CL@BC" in c]  # matches CL@BC_t_ch_MAX
    df = df.drop(columns=cols_to_drop, errors="ignore")

    # Format p-value column
    if "p-value" in df.columns:
        def format_pval(x):
            if x is None or pd.isna(x):
                return ""
            try:
                return f"{float(x):.0e}"  # scientific notation
            except:
                return str(x)

        df["p-value"] = df["p-value"].apply(format_pval)

    # Replace *, None in C1, C2, C3 with tick mark or empty
    tick = "✔️"
    for c in ["C1", "C2", "C3"]:
        if c in df.columns:
            df[c] = df[c].apply(lambda x: tick if str(x).strip() == "*" else "")

    return df

def display_detailed_report_summary(report, virus, analysis_mode):
    if not report:
        st.warning("No report data available.")
        return
    
    # summary
    if "summary" in report:
        with st.expander("Case Summary", expanded=True):
            summary = report["summary"]
            
            # Original st.metric code (commented for future reference)
            # if analysis_mode == "Consensus Sequence Analysis":
            #     [a] = st.columns(1)
            #     a.metric("Lineage Name", summary["group_name"], border=True)
            # else:
            #     a, b = st.columns(2)
            #     a.metric("Genome ID", summary["case_name"], border=True)
            #     b.metric("Lineage Name", summary["group_name"], border=True)
            # [b] = st.columns(1)
            # b.metric(f"Number of Mutations (Reference Sequence Length: {REFERENCE_LENGTHS.get(virus)})", summary["number_of_changes"], border=True)
            # [a] = st.columns(1)
            # a.metric("Recombinant Parents", summary["best_candidates"], border=True)
            # [a] = st.columns(1)
            # a.metric("Breakpoints Location in mutations-space", summary["best_candidates_breakpoints_target"], border=True)
            # [a] = st.columns(1)
            # a.metric("Breakpoints Location in Genomic Coordinates", summary["best_candidates_breakpoints_genomic"], border=True)
            # BC = summary["best_candidates"]
            # BP = f"{BC.count('+')}BP"
            # C1 = BC.split("+")[0]
            # C2 = BC.split("+")[1]
            # a, b = st.columns(2)
            # a.metric(f"Recombinant Confidence:\n{BP} Rec. vs {C1}", summary["p_value_vs_L1"], border=True)
            # b.metric(f"Recombinant Confidence:\n{BP} Rec. vs {C2}", summary["p_value_vs_L2"], border=True)
            
            # Compact metric cards using custom layout
            if analysis_mode == "Consensus Sequence Analysis":
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"""
                    <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                        <div style="font-size: 0.8rem; color: #666;">Lineage Name</div>
                        <div style="font-size: 1.0rem; font-weight: bold;">{summary["group_name"]}</div>
                    </div>
                    """, unsafe_allow_html=True)
            else:
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"""
                    <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                        <div style="font-size: 0.8rem; color: #666;">Genome ID</div>
                        <div style="font-size: 1.0rem; font-weight: bold;">{summary["case_name"]}</div>
                    </div>
                    """, unsafe_allow_html=True)
                with col2:
                    st.markdown(f"""
                    <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                        <div style="font-size: 0.8rem; color: #666;">Lineage Name</div>
                        <div style="font-size: 1.0rem; font-weight: bold;">{summary["group_name"]}</div>
                    </div>
                    """, unsafe_allow_html=True)
            
            # Number of mutations
            col3, col4 = st.columns(2)
            with col3:
                st.markdown(f"""
                <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                    <div style="font-size: 0.8rem; color: #666;">Number of Mutations</div>
                    <div style="font-size: 1.0rem; font-weight: bold;">{summary["number_of_changes"]}</div>
                    <div style="font-size: 0.7rem; color: #888;">(Ref Length: {REFERENCE_LENGTHS.get(virus)})</div>
                </div>
                """, unsafe_allow_html=True)
            
            # Recombinant Parents
            with col4:
                st.markdown(f"""
                <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                    <div style="font-size: 0.8rem; color: #666;">Recombination Pattern</div>
                    <div style="font-size: 1.0rem; font-weight: bold;">{summary["best_candidates"]}</div>
                </div>
                """, unsafe_allow_html=True)
            
            # Breakpoints
            col5, col6 = st.columns(2)
            with col5:
                st.markdown(f"""
                <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                    <div style="font-size: 0.8rem; color: #666;">Breakpoints Location in Mutations-space</div>
                    <div style="font-size: 1.0rem; font-weight: bold;">{summary["best_candidates_breakpoints_target"]}</div>
                </div>
                """, unsafe_allow_html=True)
            
            with col6:
                st.markdown(f"""
                <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                    <div style="font-size: 0.8rem; color: #666;">Breakpoints Location in Genomic Coordinates</div>
                    <div style="font-size: 1.0rem; font-weight: bold;">{summary["best_candidates_breakpoints_genomic"]}</div>
                </div>
                """, unsafe_allow_html=True)
            
            # Confidence metrics
            BC = summary["best_candidates"]
            BP = f"{BC.count('+')}BP"
            C1 = BC.split("+")[0]
            C2 = BC.split("+")[1]
            
            col7, col8 = st.columns(2)
            with col7:
                st.markdown(f"""
                <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                    <div style="font-size: 0.8rem; color: #666;">Confidence: {BP} vs {C1}</div>
                    <div style="font-size: 1.0rem; font-weight: bold;">{summary["p_value_vs_L1"]}</div>
                </div>
                """, unsafe_allow_html=True)
            
            with col8:
                st.markdown(f"""
                <div style="background-color: #f0f2f6; padding: 0.5rem; border-radius: 0.5rem; margin: 0.2rem;">
                    <div style="font-size: 0.8rem; color: #666;">Confidence: {BP} vs {C2}</div>
                    <div style="font-size: 1.0rem; font-weight: bold;">{summary["p_value_vs_L2"]}</div>
                </div>
                """, unsafe_allow_html=True)    

def  display_detailed_report(report):
    if not report:
        st.warning("No report data available.")
        return

    # region tables
    region_tables = {
        k: v
        for k, v in report.items()
        if k.startswith("region_") and isinstance(v, pd.DataFrame)
    }

    information = """
        Tables report the number of sequences,
        breakpoint, maximum likelihood ratio. The next columns illustrate the comparison
        of the candidate with the first of the table: value of one-sided AIC comparison
        between recombination model and non-recombination model (lower values are
        when row candidate is similar to the first candidate); p-value of AIC -- without
        multiple comparison corrections; and three conditions: (C1) marked if p-value is
        ≥10−5; (C2) marked if row breakpoint is at most one mutation apart from the one of
        the first candidate; (C3) marked if candidate belongs to the same phylogenetic
        branch as the first one.

        The most plausible candidates are the ones that have all three conditions (C1, C2, C3) marked.
    """
    if region_tables:
        with st.expander("Region Analysis Tables", expanded=True):
            st.write(information)
            for region_name, df in region_tables.items():
                title_mapping = {
                    "1": "First",
                    "2": "Second",
                    "3": "Third"
                }
                region = region_name.split("_")[1]
                region_title = title_mapping.get(region)
                st.markdown(f"#### {region_title} Region Candidates")

                # Format table
                formatted_df = format_region_table(df, region).reset_index(drop=True)

                # --- Filtering UI ---
                _, col0 = st.columns([11, 5])
                with col0:
                    filter_all = st.checkbox("Display only most plausible candidates", value=True, key=f"{region_name}_all")
                # with col1:
                #     filter_c1 = st.checkbox("C1", key=f"{region_name}_c1")
                # with col2:
                #     filter_c2 = st.checkbox("C2", key=f"{region_name}_c2")
                # with col3:
                #     filter_c3 = st.checkbox("C3", key=f"{region_name}_c3")


                # --- Apply filters ---
                filtered_df = formatted_df.copy()
                if filter_all:
                    filtered_df = filtered_df[
                        (filtered_df["C1"] == "✔️") &
                        (filtered_df["C2"] == "✔️") &
                        (filtered_df["C3"] == "✔️")
                    ]
                # else:
                #     if filter_c1:
                #         filtered_df = filtered_df[filtered_df["C1"] == "✔️"]
                #     if filter_c2:
                #         filtered_df = filtered_df[filtered_df["C2"] == "✔️"]
                #     if filter_c3:
                #         filtered_df = filtered_df[filtered_df["C3"] == "✔️"]

                st.dataframe(filtered_df, use_container_width=True, hide_index=True)

    # plot visualization graphs from JSON
    plot_files = [f for f in report.keys() if f.startswith("plot_") and f.endswith(".json")]
    plot_per_region = plot_files[0]
    if plot_per_region:
        with st.expander("Visualization", expanded=True):
            st.markdown(f"#### {plot_per_region.replace('_', ' ').replace('.json', '').title()}")
            fig = go.Figure(report[plot_per_region])
            st.plotly_chart(fig, use_container_width=True)

    # target mutations list
    if "target_mutations" in report:
        with st.expander("Target Mutations", expanded=True):

            mutations = [m.strip() for m in report["target_mutations"][0].split(",")]

            # --- Download button ---
            mutations_text = "\n".join(mutations)
            st.download_button(
                label="Download Mutations List",
                data=mutations_text,
                file_name="target_mutations.txt",
                mime="text/plain"
            )

            def classify_mutation(mutation: str):
                # Deletion: 12345_12349 or 12345
                if re.fullmatch(r"\d+_\d+", mutation) or re.fullmatch(r"\d+", mutation):
                    return "deletion"
                # Insertion: 12345_.|ATG
                elif re.fullmatch(r"\d+_\.\|[A-Za-z]+", mutation):
                    return "insertion"
                # Substitution: 12345_T|A
                elif re.fullmatch(r"\d+_[A-Za-z]+\|[A-Za-z]+", mutation):
                    return "substitution"
                else:
                    return "other"

            colors = {
                "deletion":   {"bg": "#ffebee", "fg": "#c62828"},   # red
                "insertion":  {"bg": "#e8f5e9", "fg": "#2e7d32"},   # green
                "substitution": {"bg": "#e3f2fd", "fg": "#1565c0"}, # blue
                "other":      {"bg": "#eeeeee", "fg": "#424242"},   # grey
            }

            def make_chip(text, kind):
                style = colors[kind]
                return (
                    f"<span style='background:{style['bg']}; color:{style['fg']}; "
                    f"padding:3px 8px; border-radius:12px; margin:2px; "
                    f"display:inline-block; font-size:90%; font-weight:500'>{text}</span>"
                )

            # --- Legend chips ---
            legend = " ".join([
                make_chip("Deletion", "deletion"),
                make_chip("Insertion", "insertion"),
                make_chip("Substitution", "substitution"),
            ])
            st.markdown(legend, unsafe_allow_html=True)

            # --- Mutation chips ---
            chips = []
            for m in mutations:
                mtype = classify_mutation(m)
                chips.append(make_chip(m, mtype))

            st.markdown(" ".join(chips), unsafe_allow_html=True)

@st.dialog("Genome Details", width="large")
def genome_details_dialog(selected_id, report, virus, analysis_mode):
    """Dialog to display detailed genome analysis."""
    st.subheader(f"Details of: {selected_id}")
    
    # Case Summary
    display_detailed_report_summary(report, virus, analysis_mode)
    
    # Region Analysis Tables
    display_region_analysis_tables(report)
    
    # Visualization
    display_visualization_section(report)
    
    # Target Mutations
    display_target_mutations_section(report)

def display_region_analysis_tables(report):
    # region tables
    region_tables = {
        k: v
        for k, v in report.items()
        if k.startswith("region_") and isinstance(v, pd.DataFrame)
    }

    information = """
        Tables report the number of sequences,
        breakpoint, maximum likelihood ratio. The next columns illustrate the comparison
        of the candidate with the first of the table: value of one-sided AIC comparison
        between recombination model and non-recombination model (lower values are
        when row candidate is similar to the first candidate); p-value of AIC -- without
        multiple comparison corrections; and three conditions: (C1) marked if p-value is
        ≥10−5; (C2) marked if row breakpoint is at most one mutation apart from the one of
        the first candidate; (C3) marked if candidate belongs to the same phylogenetic
        branch as the first one.

        The most plausible candidates are the ones that have all three conditions (C1, C2, C3) marked.
    """
    if region_tables:
        with st.expander("Region Analysis Tables", expanded=True):
            st.write(information)
            for region_name, df in region_tables.items():
                title_mapping = {
                    "1": "First",
                    "2": "Second",
                    "3": "Third"
                }
                region = region_name.split("_")[1]
                region_title = title_mapping.get(region)
                st.markdown(f"#### {region_title} Region Candidates")

                # Format table
                formatted_df = format_region_table(df, region).reset_index(drop=True)

                # --- Filtering UI ---
                _, col0 = st.columns([11, 5])
                with col0:
                    filter_all = st.checkbox("Display only most plausible candidates", value=True, key=f"{region_name}_all")
                # with col1:
                #     filter_c1 = st.checkbox("C1", key=f"{region_name}_c1")
                # with col2:
                #     filter_c2 = st.checkbox("C2", key=f"{region_name}_c2")
                # with col3:
                #     filter_c3 = st.checkbox("C3", key=f"{region_name}_c3")


                # --- Apply filters ---
                filtered_df = formatted_df.copy()
                if filter_all:
                    filtered_df = filtered_df[
                        (filtered_df["C1"] == "✔️") &
                        (filtered_df["C2"] == "✔️") &
                        (filtered_df["C3"] == "✔️")
                    ]
                # else:
                #     if filter_c1:
                #         filtered_df = filtered_df[filtered_df["C1"] == "✔️"]
                #     if filter_c2:
                #         filtered_df = filtered_df[filtered_df["C2"] == "✔️"]
                #     if filter_c3:
                #         filtered_df = filtered_df[filtered_df["C3"] == "✔️"]

                st.dataframe(filtered_df, use_container_width=True, hide_index=True)

    else:
        st.info("No region analysis tables available.")

def display_visualization_section(report):
    """Display visualization plots from the report."""
    with st.expander("Visualization", expanded=True):
        plot_files = [f for f in report.keys() if f.startswith("plot_") and f.endswith(".json")]
        if plot_files:
            plot_per_region = plot_files[0]
            st.markdown(f"#### {plot_per_region.replace('_', ' ').replace('.json', '').title()}")
            fig = go.Figure(report[plot_per_region])
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No visualization data available.")

def display_target_mutations_section(report):
    """Display target mutations from the report."""
    with st.expander("Target Mutations", expanded=True):
        if "target_mutations" in report:
            mutations = [m.strip() for m in report["target_mutations"][0].split(",")]
            
            # Download button
            mutations_text = "\n".join(mutations)
            st.download_button(
                label="Download Mutations List",
                data=mutations_text,
                file_name="target_mutations.txt",
                mime="text/plain"
            )
            
            # Mutation classification and display
            def classify_mutation(mutation: str):
                if re.fullmatch(r"\d+_\d+", mutation) or re.fullmatch(r"\d+", mutation):
                    return "deletion"
                elif re.fullmatch(r"\d+_\.\|[A-Za-z]+", mutation):
                    return "insertion"
                elif re.fullmatch(r"\d+_[A-Za-z]+\|[A-Za-z]+", mutation):
                    return "substitution"
                else:
                    return "other"

            colors = {
                "deletion":   {"bg": "#ffebee", "fg": "#c62828"},
                "insertion":  {"bg": "#e8f5e9", "fg": "#2e7d32"},
                "substitution": {"bg": "#e3f2fd", "fg": "#1565c0"},
                "other":      {"bg": "#eeeeee", "fg": "#424242"},
            }

            def make_chip(text, kind):
                style = colors[kind]
                return (
                    f"<span style='background:{style['bg']}; color:{style['fg']}; "
                    f"padding:3px 8px; border-radius:12px; margin:2px; "
                    f"display:inline-block; font-size:90%; font-weight:500'>{text}</span>"
                )

            # Legend chips
            legend = " ".join([
                make_chip("Deletion", "deletion"),
                make_chip("Insertion", "insertion"),
                make_chip("Substitution", "substitution"),
            ])
            st.markdown(legend, unsafe_allow_html=True)

            # Mutation chips
            chips = []
            for m in mutations:
                mtype = classify_mutation(m)
                chips.append(make_chip(m, mtype))

            st.markdown(" ".join(chips), unsafe_allow_html=True)
        else:
            st.info("No target mutations data available.")

def create_recombinant_cases_table(df, virus, analysis_mode):
    """Create a table to display recombinant cases."""
    st.subheader("Recombinant Cases")

    if df.empty:
        st.warning("No recombinant cases found.")
        return
    
    # Full width table
    with st.spinner("Loading recombinant cases..."):
        if analysis_mode == "Consensus Sequence Analysis":
            formatter = {
                "original_lineage": ("Lineage", {"width": 150}),
                "breakpoint_count": ("BP Count", {"width": 80}),
                "recombinant_parents": ("Recombinant Parents", {"width": 250}),
            }
        else:
            formatter = {
                "genomeID": ("Genome ID", PINLEFT),
                "breakpoint_count": ("BP Count", {"width": 80}),
                "original_lineage": ("Assigned Lineage", {"width": 150}),
                "recombinant_parents": ("Recombinant Parents", {"width": 250}),
                "country": ("Country", {"width": 100}),
                "collection_date": ("Collection Date", {"width": 100}),
            }

        custom_css = {
            ".ag-root": {"font-family": "inherit"}, 
            ".ag-cell": {"font-family": "inherit"},
            ".ag-header-cell": {"font-family": "inherit"}
        }

        response = draw_grid(
            df,
            formatter=formatter,
            fit_columns=True,
            selection="single",     
            use_checkbox=True,     
            max_height=600,
            css=custom_css,
        )

    # Handle table selection and show dialog
    if response:
        selected = response["selected_rows"]
        if selected is not None:
            if analysis_mode == "Consensus Sequence Analysis": 
                selected_id = selected["original_lineage"].iloc[0]
            else: 
                selected_id = selected["genomeID"].iloc[0]
            
            # Load report data
            path_to_the_case_report_folder = selected["case_report_folder"].iloc[0]
            report = load_report_data(path_to_the_case_report_folder)
            
            # Show dialog with details - handle potential conflicts gracefully
            try:
                genome_details_dialog(selected_id, report, virus, analysis_mode)
            except Exception as e:
                if "Only one dialog is allowed" not in str(e):
                    st.warning(f"Error opening dialog: {e}")
        else:
            word = "lineage" if analysis_mode == "Consensus Sequence Analysis" else "genome"
            st.info(f"Select a recombinant {word} from the table above to see its details.")

def sidebar(virus_list):
    """sidebar navigation for the streamlit"""
    with st.sidebar:

        menu_options = ["Home"] + sorted([visualize(v) for v in virus_list if v in mapping.keys()])
        menu_icons = ["house"] + ["virus2"] * (len(menu_options) - 1)

        selected = option_menu(
            menu_title="OpenRecombinHunt",
            options=menu_options,
            icons=menu_icons,
            menu_icon="",
            default_index=0,
            orientation="vertical",
            styles={
                "container": {"padding": "0!important"},
                "icon": {"color": "#000000", "font-size": "18px"},
                "nav-link": {"font-size": "16px", "text-align": "left", "margin": "0px", "--hover-color": "#818181"},
                "nav-link-selected": {"background-color": "#4A90E2"},
                "menu-title": {"display": "none"},  # Hide the menu title
            }
        )

        return selected


def show_home_page():
    """Displays the new, comprehensive welcome page."""
    
    st.title("OpenRecombinHunt: Automatic Detection of Recombination From Publicly Available Sequences")
    st.markdown("---")

    # 1. A More Compelling Introduction
    st.header("Welcome to the OpenRecombinHunt Analysis Dashboard")
    st.markdown("""
    This application provides an interactive interface to explore the results of the OpenRecombinHunt pipeline, an automated tool for detecting recombination in viral genomes. The goal of this project is to empower researchers and public health officials by making complex genomic analysis accessible and easy to interpret.

    Use the menu on the left to select a virus and begin your exploration.
    """)

    # 2. The Pipeline Diagram (Now with PNG)
    st.header("The OpenRecombinHunt Pipeline")
    
    png_path = "app/entire-pipeline.png"
    try:
        # Display the PNG image directly
        st.image(png_path, use_container_width=False)
        
        # Keep the download button for user convenience
        with open(png_path, "rb") as png_file:
            st.download_button(
                label="Download Pipeline Diagram (PNG)",
                data=png_file,
                file_name="OpenRecombinHunt_Pipeline.png",
                mime='image/png'
            )
    except FileNotFoundError:
        st.warning("Pipeline diagram PNG not found. Please ensure 'app/entire-pipeline.png' exists.")

    st.markdown("---")

    # 3. Interactive Pipeline Explanation
    with st.expander("Module 1: Data Acquisition"):
        st.markdown("""
        This module is responsible for the automated download of raw data. It uses a central configuration file to fetch metadata and sequences from public databases like NCBI and Nextstrain, handling various download methods (CLI, URL, FTP) to produce a standardized set of raw files.
        """)
    with st.expander("Module 2: Preprocessing"):
        st.markdown("""
        Raw data is subjected to a rigorous, source-specific preprocessing workflow. This module cleans, filters, and standardizes the data, applying quality control rules defined in the configuration file to ensure only high-quality, complete records are used for analysis.
        """)
    with st.expander("Module 3: HaploCoV"):
        st.markdown("""
        For viruses without an existing nomenclature, this module uses HaploCoV to perform *de novo* lineage classification. It clusters sequences into "haplogroups" based on shared mutation profiles, providing the essential lineage assignments needed for the core analysis. For viruses with existing classifications, it can also augment them by identifying novel sub-clusters.
        """)
    with st.expander("Module 4: Postprocessing"):
        st.markdown("""
        This module standardizes the mutation notation from different sources (Nextstrain and HaploCoV) into a single, consistent format required by the RecombinHunt tool. It correctly parses substitutions, insertions, deletions, and complex "compound" mutations.
        """)
    with st.expander("Module 5: Prepare for RecombinHunt"):
        st.markdown("""
        Before the final analysis, this module prepares two key sets of inputs: the environment, which characterizes the genetic landscape of the virus, and the per-lineage samples that will be tested for recombination.
        """)
    with st.expander("Module 6: RecombinHunt"):
        st.markdown("""
        This is the core analysis module. It uses the prepared environment and samples to run the RecombinHunt tool, a data-driven method that uses a statistical framework to detect genomes with one or two recombination breakpoints.
        """)
    with st.expander("Module 7: Streamlit"):
         st.markdown("""
        This is the final visualization layer of the pipeline. The Streamlit application you are currently using reads the outputs from the pipeline and presents them in a structured, interactive dashboard format for exploration.
        """)

    st.markdown("---")
    
    # 4. Core Technologies Section
    st.header("Core Technologies")
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("What is HaploCoV?")
        st.markdown("""
        HaploCoV is a software framework for the unsupervised classification of viral variants. It clusters viral genomes into "haplogroups" based on shared, high-frequency mutations, making it ideal for assigning lineages to viruses that lack an established nomenclature.
        
        *Reference: Chiara, M., et al. (2023). Commun Biol.*
        """)
    with col2:
        st.subheader("What is RecombinHunt?")
        st.markdown("""
        RecombinHunt is a data-driven method for identifying recombinant viral genomes. It uses a likelihood-based approach to compare recombinant and non-recombinant models, allowing it to detect mosaic genomes with high accuracy and within reduced turn-around times.

        *Reference: Alfonsi, T., et al. (2024). Nat Commun.*
        """)

    st.markdown("---")

    # 5. "How to Use" and "About" Sections
    st.header("About This Project")
    st.markdown("""
    This dashboard is the final output of a Master's thesis project in Computer Science and Engineering at Politecnico di Milano. The entire pipeline is designed to be fully automated and is run on a monthly schedule to provide up-to-date analysis. The source code is publicly available on [GitHub](https://github.com/yavuuzsameet/openrecombinhunt/tree/master).
    
    **Author:** Yavuz Samet Topcuoglu
    
    **Advisors:** Prof. Anna Bernasconi, Dr. Tommaso Alfonsi
    """)

def show_virus_page(virus):
    """display virus-specific analysis and visualizations"""
    virus_name = visualize(virus)
    st.title(f"{virus_name} Dashboard")

    loader_placeholder = st.empty()
    loader_placeholder.markdown(loading_animation(virus_name), unsafe_allow_html=True)

    # load master data
    with st.spinner(f"Loading data for {virus}..."):
        master_df = load_master_data(virus)
        loader_placeholder.markdown("<style>.loader-overlay{display:none;}</style>", unsafe_allow_html=True)

    if master_df is None or master_df.empty:
        st.error(f"No data found for {virus}.")

        # display a card with information
        st.markdown(f"""
        <div style="border: 1px solid #eee; border-radius: 5px; padding: 10px;">
            <h4>{virus_name} is in the oven...</h4>
            <p>Please wait for some time.</p>
        </div>
        """, unsafe_allow_html=True)    

        return

    # tabs - different structure for SARS-CoV-2
    if virus == "sars-cov-2":
        # Get analysis window months from config for dynamic tab naming
        try:
            virus_config = config.get(VIRUSES).get(virus)
            analysis_window_months = virus_config.get("analysis_window_months", 6)
        except Exception:
            analysis_window_months = 6
        
        tab0, tab1, tab2, tab3 = st.tabs([
            "About the Virus", 
            f"Summary Dashboard - Last {analysis_window_months} Months", 
            f"Recombinant Explorer - Last {analysis_window_months} Months", 
            "Recombinant Explorer - Consensus Sequences"
        ])
    else:
        tab0, tab1, tab2 = st.tabs(["About the Virus", "Summary Dashboard", "Recombinant Explorer"])

    with tab0:
        st.header(f"About {virus_name}")

        if "stats" not in st.session_state:
            st.session_state.stats = {}

        if virus == "sars-cov-2":
            if "sars-cov-2" not in st.session_state.stats:
                loader_placeholder = st.empty()
                loader_placeholder.markdown(loading_animation(virus_name), unsafe_allow_html=True)
                stats = load_complete_data_stats(virus)
                st.session_state.stats["sars-cov-2"] = stats
                loader_placeholder.markdown("<style>.loader-overlay{display:none;}</style>", unsafe_allow_html=True)
            else:
                stats = st.session_state.stats["sars-cov-2"]

            describe(virus, config, stats)
        else: describe(virus, config, master_df)

    with tab1:
        if virus == "sars-cov-2":
            st.info(f"Due to the vast amount of SARS-CoV-2 data, the Summary Dashboard is limited to the most recent {analysis_window_months} months of sequences. For a comprehensive analysis of available SARS-CoV-2 sequences, please utilize the Recombinant Explorer tabs.")

        # time-based filtering
        with st.spinner("Applying time filter..."):
            summary_df = apply_time_filter(master_df, virus)

        st.markdown("---")

        # create key metrics and display
        with st.spinner("Creating key metrics..."):
            create_key_metrics(summary_df)

        st.markdown("---")

        # create summary tables and display
        create_summary_tables(summary_df)

        st.markdown("---")

        create_distribution_plots(summary_df, virus)

    # Handle different tab structures based on virus
    if virus == "sars-cov-2":
        # SARS-CoV-2 has 4 tabs - separate recombinant explorer tabs
        with tab2:
            # Last X months analysis (genome-level)
            st.info(f"This tab shows recombinant cases from the last {analysis_window_months} months analysis, allowing for genome-level analysis of the most recent data.")
            
            # filtering
            with st.spinner("Applying filters..."):
                recombinant_df = master_df[master_df["is_recombinant"]]
                explorer_df = apply_user_filter(recombinant_df, virus,)

            st.markdown("---")

            # create interactive table with radio buttons as the index column
            create_recombinant_cases_table(explorer_df, virus, None)

        with tab3:
            # Consensus sequence analysis (lineage-level)
            st.info("This tab shows recombinant cases from consensus sequence analysis. The consensus sequence analysis encompasses all available sequences that belong to the same lineage into a 'consensus sequence' and allows for lineage-level analysis rather than genome-level analysis.")
            
            df = load_consensus_data(virus)
            create_recombinant_cases_table(df, virus, "Consensus Sequence Analysis")
    else:
        # Other viruses have 3 tabs - original structure
        with tab2:
            # filtering
            with st.spinner("Applying filters..."):
                recombinant_df = master_df[master_df["is_recombinant"]]
                explorer_df = apply_user_filter(recombinant_df, virus,)

            st.markdown("---")

            # create interactive table with radio buttons as the index column
            create_recombinant_cases_table(explorer_df, virus, None)

    if "initial_rerun_done" not in st.session_state:
        st.session_state.initial_rerun_done = True
        st.rerun()

def main():
    # discover available viruses
    viruses = discover_viruses()

    viruses_visualized = [visualize(v) for v in viruses]

    # sidebar navigation
    selected = sidebar(viruses)

    if selected == "Home":
        show_home_page()
    elif selected in viruses_visualized:
        show_virus_page(viruses[viruses_visualized.index(selected)])

if __name__ == "__main__":
    main()