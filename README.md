# Monte Carlo HP Protein Folding
This program is a implementation of a Monte Carlo algorithm for protein folding as described in the paper _"A replica exchange Monte Carlo algorithm for protein folding in the HP model"_, by THACHUK C. and _al_. (2007)
## Installation

Clone this repository:

```
git clone git@github.com:
```
Move to the directory:

```
cd MC_protein_folding
```
## Conda Environment
Install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Create the `projet_mc` conda environment:

```
conda env create -f environment_monte_carlo.yml
```

Load the `projet_mc` conda environment:

```
conda activate projet_mc
```
## Running Monte Carlo Protein folding

```
python fold_protein.py -f data/sample2.fasta -n 1000 --display-graph optimal --display-final-frame --initial-conformation random```
