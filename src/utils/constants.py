# src/utils/constants.py

# This file contains constant variables for keys used throughout the pipeline,
# primarily for accessing values from the config.yaml file.
# Using these constants prevents typos and makes the code more maintainable.

# --- Top-level keys in config.yaml ---
PATHS = "paths"
VIRUSES = "viruses"
NCBI_CLI = "ncbi-cli"
NEXTSTRAIN_URL = "nextstrain-url"
FTP_URL = "ftp-url"
FTP_CLI = "ftp-cli"

# --- Keys for path configuration ---
RAW_DATA = "raw_data"
PROCESSED_DATA = "processed_data"
RESULTS = "results"
RECOMBINHUNT_OUTPUT = "recombinhunt_output"
NEXTSTRAIN_OUTPUT = "nextstrain_output"
HAPLOCOV_OUTPUT = "haplocov_output"
ENVIRONMENTS = "environments"
SAMPLES = "samples"
LOGS = "logs"
LIBS = "libs"
SRC = "src"

# --- Keys within each virus's configuration ---
NAME = "name"
SOURCE = "source"
NCBI = "ncbi"
NEXTSTRAIN = "nextstrain"
FTP = "ftp"
METHOD = "method"
CLI = "cli"
URL = "url"
TAXON_ID = "taxon_id"
REFERENCE = "reference"
PARAMETERS = "parameters"

# --- Keys within the 'reference' block ---
ACCESSION_ID = "accession_id"

# --- Keys for tool parameters ---
PARAMETERS = "parameters"
METADATA_PROCESSING = "metadata_processing"
HAPLOCOV = "haplocov"
RECOMBINHUNT = "recombinhunt"

# --- Keys for Metadata Processing parameters ---
SOURCE_LINEAGE_COLUMN = "source_lineage_column"
FILTERS = "filters"
COLUMN = "column"
VALUE = "value"
OPERATOR = "operator"

# --- Keys for HaploCov parameters ---
DIST = "dist"
SIZE = "size"
DESIGNATION_MODE = "designation_mode"
HAPLOCOV = "haplocov"
NOMENCLATURE = "nomenclature"
BOTH = "both"

# --- Keys for RecombinHunt parameters ---
LC_THRESHOLD = "lc_threshold"
MIN_GENOME_COUNT = "min_genome_count"
RUN_MODE = "run_mode"
NUM_SAMPLES_TO_RUN = "num_samples_to_run"
CONSENSUS_THRESHOLD = "consensus_threshold"
RANDOM = "random"
CONSENSUS = "consensus"
ALL = "all"

# --- Keys for CLI command blocks ---
SEQUENCES = "sequences"
# The 'REFERENCE' key is already defined above and can be reused here.

# --- Keys for Nextstrain URL parameters ---
METADATA = "metadata"
SEQUENCES = "sequences"

