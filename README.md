# OpenRecombinHunt: An Automated Recombination Analysis Pipeline

OpenRecombinHunt provides an online temporal snapshot of potential recombination events (i.e., single sequences with suspected mosaic structure, with one or two breakpoints).
In this repository, we provide the code of our method and corresponding Web Server, with the potential of unlocking the automatic detection of recombination in viruses along the most current genomic surveillance interests.

#### Motivation
Viruses undergo change affecting their genome by several mechanisms, including point mutation and recombination. With the availability of open databases with large amounts of genome sequences (e.g., NCBI Virus and Nextstrain/pathogens) and the implementation of genomic surveillance systems, the need for light-weight automatic computational methods for  monitoring  continuously updating open data surges.

#### Methodology
OpenRecombinHunt extends our previously published RecombinHunt method ([Alfonsi et al., 2024](https://doi.org/10.1038/s41467-024-47464-5)), which we extensively applied to identify recombinant SARS-CoV-2 lineages, to any virus for which a large corpus of sequences is publicly available. 
Here, we couple RecombinHunt with the HaploCov method ([Chiara et al., 2023](https://doi.org/10.1038/s42003-023-04784-4)), which allows the stratification of any virus in a discrete number of groups, even in the absence of a reference nomenclature, and derives a list of characterizing mutations. 

OpenRecombinHunt provides 1) the holistic framework that exploits the output of HaploCoV as input of RecombinHunt, and 2) the automatic pipeline and Web reporting system that regularly updates datasets from NCBI Virus/Nextstrain sources and identifies novel recombination events, with a substantial contribution with respect to the previously cited works.

#### Use cases 
We apply this framework to openly-accessible datasets of SARS-CoV-2, Respiratory syncytial virus (RSV) A/B, monkeypox, Zika, Yellow Fever, and even hemagglutinin segments of H5N1 Influenza A, reporting several interesting insights.


 --------

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
## System requirements
- a terminal / command prompt
- Docker CLI v2 or [Docker Desktop >=4.0.0](https://docs.docker.com/get-started/get-docker/);
- storage space: 8GB minimum; 16GB recommended

> [!NOTE] Storage space
> The software image alone requires ~4GB without considering any viral data. Because the pipeline downloads all the available sequences for the viruses configured in the `config.yaml` file (see Section [Configuration](#configuration) for more info), the storage requirements depend on the configuration and are generally destined to increase over time. In Docker Desktop, it's always possible to increase/decrease the maximum amount of storage allowed for the installed images. 

## Setup and Installation

**All the commands are intended to be used in a terminal window inside the directory where this file resides.**

Software dependencies are listed in the environment.yml and installed through Conda and PIP within a Docker container (openrecombinhunt-base) that prepares the vitual environment for running the software. To prepare the virtual environment, start Docker and run:
```bash
docker compose build base && docker compose build
```

<small>Notes for software developers: whenever a change to the library files or dependencies (i.e., the `libs/` directory or the `environment.yml` file) is made, rebuild the virtual environment with the option `--no-cache` to apply the change.</small>

## Usage

This package can be executed in two **modes**:
- as a [**pipeline**](#running-the-pipeline) to download and analyse the viral data
- as a **web application** to display the recombinant analyses' results

Which mode to run?
- "I want to use the software for the first time": run the pipeline to build the initial output dataset. Then visualize it through the web application. 
- "I already used the web app, and want to look at it again": you can restart only the web application. 
- "I want updated recombinant analyses on more recent data": if the data source contains more recent data, you can re-run the pipeline to download it and update the previous results. Then, re-run the web app to visualize it. 

### Run OpenRecombinHunt as a Pipeline (updating the data)

#### Pipeline options
The entire pipeline is controlled by the `config/config.yaml` file which defines:
- the paths for the data, results, logs, etc.
- the viruses to include
- for each virus: the data quality filters, download source, method, taxon ID, reference seqeunce, and parameters for HaploCov and RecombinHunt.

#### Running the pipeline
In a terminal, run:
```bash
docker compose run --rm pipeline
```
The above command will open a shell in the virtualized environment. There you can execute the entire pipeline as:

```bash
python src/00_master/pipeline.py --virus yellow-fever   # only for a virus (e.g. yellow fever)
# or
python src/00_master/pipeline.py --virus all            # for all viruses configured in config/config.yaml
```

The pipeline writes the result of download and analyses in the folders 
- `samples/` = sequences to anlyse through RecombinHunt
- `results/` = results of the RecombinHunt analyses
- `environments/` = RecombinHunt's environments (needed for the analyses)

Once the analyses are completed, close the terminal by typing
```bash
exit
```

> [!WARNING] Pipeline execution times
> The pipeline requires from few minutes (for the smallest dataset) to several hours (e.g., > 10 hours for SARS-CoV-2).


### Run OpenRecombinHunt as a Web Application (display the analyses)

In a terminal, run:
```bash
docker compose build frontend && docker compose up frontend
```

The above command will copy the content of the `results/` folder within the virtual environment and start a web-server accessible through a browser at the address [http://localhost:60129](http://localhost:60129).

To stop the web server, press `Ctrl+C` or `Cmd+C` in the same terminal window where you started the web application. 

Note that you can also run the pipeline (update the data) while the web-server is running, but the updates will be visible only after the web application is rebooted. 

## Uninstall 
Open Docker Desktop and delete the images and containers related to OpenRecombinHunt. Then delete this repository. 

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.
