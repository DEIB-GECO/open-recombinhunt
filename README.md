# OpenRecombinHunt: An Automated Recombination Analysis Pipeline

## Overview
OpenRecombinHunt is an end-to-end bioinformatics pipeline designed to automate the detection and analysis of viral recombination events. The pipeline handles everything from data acquisition from public repositories like NCBI and Nextstrain, to data preprocessing, environment creation, running recombination analysis using HaploCoV and RecombinHunt, and generating final reports and visualizations.

The entire workflow is configurable via a central config/config.yaml file and is designed to be run on a monthly schedule to process the latest available data.

## Directory Structure
The project is organized into a modular structure to separate code, data, configuration, and results.


```
openrecombinhunt/
|
├── app/
|   └── streamlit_app.py      # (Future) The Streamlit web server application.
|
├── config/
|   └── config.yaml           # Main configuration for paths, viruses, and tool parameters.
|
├── data/
|   ├── raw/{virus_name}/       # Original, immutable data from NCBI, Nextstrain, etc.
|   └── processed/{virus_name}/ # Cleaned & prepared data, ready for haplocov input.
|
├── environments/
|   └── {virus_name}/         # Recombinhunt environments, organized by virus.
|
├── samples/
|   └── {virus_name}/         # Per-lineage sample files for RecombinHunt.
|
├── libs/
|   ├── haplocov/                # Source code for the HaploCoV tool.
|   └── recombinhunt-cov-7.0.0/  # Source code for the RecombinHunt tool.
|
├── results/
|   |── analysis
|   |── nextstrain_output
|   ├── haplocov_output/
|   └── recombinhunt_output/
|
├── src/
|   ├── 00_master/pipeline.py                   # The main orchestrator script for the entire pipeline.
|   ├── 01_data_acquisition/fetch_data.py       # Script for downloading data.
|   ├── 02_preprocessing/                       # Scripts for cleaning and preparing data.
|   ├── 03_haplocov/run_haplocov.py             # Script for running the HaploCoV tool.
|   ├── 04_postprocessing/                      # Scripts for formatting the tsv files before recombinhunt.
|   ├── 05_prepapre_recombinhunt/               # Scripts for creating the environment and the samples for recombinhunt.
|   ├── 06_recombinhunt/run_recombinhunt.py     # Script for running the Recombinhunt tool.
|   └── utils/                                  # Shared helper functions and constants.
|
├── logs/
|   └── ... (Log files for each pipeline run are stored here).
|
├── environment.yml         # Conda environment definition.
└── README.md               
```

## Setup and Installation
All the commands are intended to be used in a terminal window inside the directory where this file resides. 

Software dependencies are listed in the environment.yml and installed through Conda and PIP within a Docker container (openrecombinhunt-base) that prepares the vitual environment for running the software. To prepare the virtual environment, Start Docker and simply run:
```bash
docker compose build
```

## Configuration
The entire pipeline is controlled by the config/config.yaml file. Before running, you must configure:

Paths: Define the base paths for data, results, logs, etc.

Viruses: Add or modify entries for each virus, specifying the data source, method, taxon ID, reference sequence, and parameters for HaploCoV and RecombinHunt.

## How to Run the Pipeline (updating the data)
First, start the dedicated virtual environment by running the command:
```bash
docker compose run --rm pipeline
```
The above command will open a terminal window in the virtualized environment, where you can execute the pipeline commands.

The pipeline is executed via the main orchestrator script src/00_master/pipeline.py. It should be run from the root directory of the openrecombinhunt project. 

To run the entire pipeline for a single virus:
```bash
python src/00_master/pipeline.py --virus yellow-fever
```
To run the entire pipeline for all viruses defined in the config file:
```bash
python src/00_master/pipeline.py --virus all
```

Some files and folders in the following directories will be overwritten as a result of the pipeline updates:
- environments/
- results/
- samples/

## Start the Openrecombinhunt web app (inspect the data)
Run
```bash
docker compose up frontend
```

The above command will copy the content of the results folder within the virtual environment and start a web-server accessible through a browser at the address [http://localhost:60129](http://localhost:60129).

To stop the web server, type Ctrl+C or Cmd+C. 

Note that you can also run the pipeline (update the data) while the web-server is running. 

## Additional notes
Whenever a change is applied to the library files in `libs/` or to the dependencies in `environment.yml`, you must run the following command to update the virtual environment
```bash
docker compose build --no-cache
```
