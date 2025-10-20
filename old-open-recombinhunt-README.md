# OpenRecombinHunt

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

#### Repository essentials
This repository will contain:
- The code of the Web application running on our servers
- A user manual to run the application on a local machine:
    - Docker architecture
    - Software download
    - Preparation/format of input
    - Application startup
In this preliminary version, we only report the documentation of our essential module, RecombinHunt, below.

----

# RecombinHunt
RecombinHunt is a Python library implementing a data-driven novel method for identifying contributing lineages and breakpoints in recombinant viral sequences. The method is described in the following manuscript:
> [Data-driven recombination detection in viral genomes](https://doi.org/10.1038/s41467-024-47464-5),
> Tommaso Alfonsi, Anna Bernasconi, Matteo Chiara, Stefano Ceri> 
> Nature Communications 15, 3313 (2024). https://doi.org/10.1038/s41467-024-47464-5
The code is available on a [Zenodo repository](https://doi.org/10.5281/zenodo.8123832.

## Installation
Installation requires Python 3.10.12 and PIP. The software is independent of the operating system.
It is suggested to use a dedicated python environment (e.g., conda, miniconda or venv). Below, it is described how to create one with conda.

#### System requirements
Here we describe how to create a conda environment suitable for the installation of RecombinHunt. If you already know how to create one or want to use a different virtual environment, you can safely skip this subsection.

1. Follow the instructions at https://docs.conda.io/en/latest/miniconda.html# to download and install the latest version of miniconda.
2. Create and activate a dedicated conda environment with Python 3.10.12 and PIP
    ```bash
   $ conda create -n rh_env python=3.10.12 pip
   $ conda activate rh_env
    ```
    
Once the prerequisites are satisfied, move into the RecombinHunt-CoV directory and install the package with:  

```bash
$ pip install recombinhunt-7.0.0-py3-none-any.whl
```

The installation procedure will take ~ 1 minute (depending on internet connection speed) and install the following packages:
```
python 3.10
numpy 1.26.0
pandas 2.1.1
plotly 5.17.0
tqdm 4.66.1
inflect 7.0.0
tabulate 0.9.0
kaleido 0.2.1
jupyter 1.0.0
fastjsonschema 2.20.0
pyarrow 17.0.0
```

## Usage
This package already provides the context information (*environment*) that are needed to evaluate the sequences. (The given *environment* is compressed to save storage space; please unzip it before proceeding). 

You can load any suitable environment as:
```python
from recombinhunt.core.environment import Environment
env = Environment("environments/env_nextstrain_2023_03_30") # <- path to the unzipped environment folder
```

At the core of the package is the *Experiment* class, which analyses a single sequence and detects if the input is a recombination, the contributing 
lineages and the breakpoint position. 
To run an Experiment, you need to provide the *Environment* and a target sequence:

```python
from recombinhunt.core.method import Experiment

target = ['22029_22034', '28248_28253', '28271_28271', '22204_.|GAGCCAGAA', '210_G|T', '241_C|T', '1321_A|C', '3037_C|T', '4181_G|T', '4890_C|T', '6402_C|T', '7124_C|T', '7851_C|T', '8723_A|G', '8986_C|T', '9053_G|T', '10029_C|T', '11201_A|G', '11332_A|G', '14407_C|T', '14408_C|T', '15264_T|C', '15451_G|A', '16466_C|T', '18366_A|G', '19220_C|T', '20032_C|T', '21618_C|G', '21641_G|T', '21846_C|T', '21987_G|A', '22578_G|A', '22673_T|C', '22674_C|T', '22679_T|C', '22686_C|T', '22813_G|T', '22882_T|G', '22898_G|A', '22992_G|A', '22995_C|A', '23013_A|C', '23040_A|G', '23048_G|A', '23055_A|G', '23063_A|T', '23075_T|C', '23202_C|A', '23403_A|G', '23525_C|T', '23599_T|G', '23604_C|A', '23854_C|A', '23948_G|T', '24130_C|A', '24424_A|T', '24469_T|A', '24503_C|T', '25000_C|T', '25667_C|T', '25855_G|T', '26767_T|C', '27638_T|C', '27752_C|T', '27874_C|T', '28461_A|G', '28881_G|T', '28916_G|T', '29402_G|T', '29540_G|A', '29645_G|T', '29742_G|T'] # = list of mutations in the target sequence

experiment = Experiment(environment=env)
experiment.set_target(target)
result = experiment.run()
```

### Output
Results can be displayed by simply calling ```print(result) ```. An example output looks like this:
```json
target length : 69 
designated candidates :  AY.4 + BA.1.15.3 + AY.4  
region details :   1 
                    pos_start_in_t : 1 
                    pos_end_in_t : 25 
                    designated :  AY.4
                  2 
                    pos_start_in_t : 26 
                    pos_end_in_t : 54 
                    designated :  BA.1.15.3
                  3 
                    pos_start_in_t : 55 
                    pos_end_in_t : 69 
                    designated :  AY.4
AIK :  AY.4 : 350.2372865177939 
       BA.1.15.3 : 1452.8035425553069
       AY.4 + BA.1.15.3 + AY.4 : -412.89739429828103  
p_values :  AY.4 + BA.1.15.3 + AY.4 vs AY.4 :  1.93e-166  
            AY.4 + BA.1.15.3 + AY.4 vs BA.1.15.3 :  0.00e+00
```

Likelihood ratio can be visualized with:
```python
from recombinhunt.core.graphics import *

plot_likelihood(result.genome_view, xaxis="changes")
```

![alt text](next_XD_plot.png "Plot of likelihood ratio for case XD")

### Customize the Experiment's parameters
The default behavior of RecombinHunt can be modified by overriding the parameters in the Experiment constructor method: 
```python
Experiment(
  ...
  min_searchable_region_length=3,
  min_candidate_region_length=3,
  min_l2_enclosed_region_length=2,
  alt_candidate_p_value_difference=1e-05,
  alt_candidate_max_pos_distance_t=1)
```

Default values for such parameters are stored in `recombinhunt.core.method.DefaultParams`.

### Ignore specific candidate variants
The knowledge of specific variants/lineages can be ignored, if necessary,  by altering the Environment like so:

```python
# explicitly name the candidates to ignore 
ignored_candidates = ['XBB.1', 'XBB.2']
# or filter the available candidates
all_candidates = base_environment.included_lineages()
ignored_candidates = [l for l in all_candidates if not l.startswith('XBB.')]

# clone an existing environment and remove those candidates (quicker method)
custom_environment = base_environment.copy_with_exclusions(ignored_candidates)
# or create an environment without those canididates 
custom_environment = Environment("environments/env_nextstrain_2023_03_30", ignored_candidates)

# run an experiment...
```

This possibility is helpful if the analysed target sequence is recognised as a descendant of a recombinant variant (e.g., XBB.1), while the desired output should report, instead, the recombinant parental candidates (i.e., BJ.1 + BM.1.1.1)

## Demo

In the ```demo/``` directory, you find the Jupyter Notebook ```recombinant_cases_nextstrain.ipynb```. This 
notebook computes the recombinant cases in Nextstrain dataset using the consensus of all the sequence of 
good quality found for each recombinant case. 

### Demo input
The original nucleotide sequences are
stored in ```demo/demo_input_nextstrain``` - for example, the sequences of recombinant case XD are stored in 
```demo/demo_input_nextstrain/sampels_XD.csv```. For each case, the consensus sequence is computed at runtime.

### Demo output
The notebook produces two files stored in ```demo/demo_input_nextstrain```:
- *summary.md* is a markdown file organising the output of RecombinHunt in a tabular form (one row for each case) and 
comparing the output against the ground truth when available.
- *detail.html* is an HTML file that can be viewed in a browser (internet connection is required to load the 
plotting library). This file contains a more detailed output of RecombinHunt; it includes the plot of the likelihood 
ratio for the candidate lineages contributing to a recombination, and the consensus sequence for each case.

### Expected run time
The demo runs in ~ 1 minute. 

### Instructions
The demo is a Juputer notebook and requires a Jupyter server to run.
The RecombinHunt package will automatically install Jupyter among its dependencies. To start a jupyter server locally, 
open a terminal, move inside this README is located and run:
```bash
$ jupyter notebook
```
In case the browser doesn't open automatically, you can click on the link printed in the terminal (the URL 
will be similar to http://localhost:8888/?token=5f38de823...). Once the browser starts, navigate to the demo directory
and execute every cell of the notebook.

## Data

This package already contain some data in order to ease the testing of the software.
Included data files are: 
-  ```demo/demo_input_nextstrain```: nucleotide sequences of recombinant cases downloaded from Nextstrain. Only the 
sequences satisfying our quality filters were retained.
- ```environments/env_nextstrain_2023_03_30```: information about the probability of nucleotide changes globally and 
for each lineage.
- ```demo/validation_data/alias_key.json```: File of recombinant definitions provided from PANGO GitHub repository.

## Source code

The source code is located in the ```src/``` directory.


 --------

# HaploCov
HaploCoV is a library collecting a set of utilities and methods to identify novel variants of viruses.
In the manuscript:
> [HaploCoV: unsupervised classification and rapid detection of novel emerging variants of SARS-CoV-2](https://doi.org/10.1038/s42003-023-04784-4),
> Matteo Chiara, David S. Horner, Erika Ferrandi, Carmela Gissi, Graziano Pesole> 
> Communication Biology 6, 443 (2023). https://doi.org/10.1038/s42003-023-04784-4
Its use has been described for SARS-CoV-2.
The code is available on a [GitHub repository](https://github.com/matteo14c/HaploCoV).



--------

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.


