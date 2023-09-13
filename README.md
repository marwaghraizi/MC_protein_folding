# Monte Carlo HP Protein Folding
This program is a implementation of a Monte Carlo algorithm for protein folding as described in the paper _"A replica exchange Monte Carlo algorithm for protein folding in the HP model"_, by THACHUK C. and _al_. (2007)
## Installation

Clone this repository:

```
git clone git@github.com:marwaghraizi/MC_protein_folding
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
python fold_protein.py -f data/sample2.fasta -n 1000 --display-graph optimal --display-final-frame --initial-conformation random
```
### Argument Options

|           argument           | description                                                                                                                            |
| :--------------------------: | :------------------------------------------------------------------------------------------------------------------------------------- |
|          -h, --help          | Displays the help.                                                                                                                     |
|          -f, --file          | Fasta file containing the protein to fold in HP format or classic amino acid format                                                    |
|          -p, --protein       | Protein sequence to fold in HP format or classic amino acid format.                                                                    |
|  -i, --initial-conformation  | Initial protein conformation: linear or randomized. Default: linear.                                                                   |
|      -n,--n-iterations       | Number of iterations.                                                                                                                  |
|      -t, --temperature       | Search temperature.                                                                                                                    |
|    -s, --search-space        | Moves search neighborhood: VSHD (corner, crankshaft and end moves), pull (pull move only), VSHD-pull (all four moves).                 |
|      --probability-pull      | In the case of the VSHD-pull hybrid search space, specify the probability of applying the pull move.                                   |
|      --display-final-frame   | Display final protein frame as directions trace.                                                                                       |
|      --display-graph         | Creates a .PNG with a visualization of the conformation with hydrophobic residues in blue and polar residues in red: The options are final or optimal and the default is the final conformation.                                     |

## Outputs
The program generates a **log file** detailing the protein HP sequence, the initial energy, the final energy and the optimal energy. It also outputs the detailed result of each iteration by outputting the protein in a traced directions format as well as the corresponding energy and the move that caused the new change. In the case of a randomized initial conformation, the program saves an image containing a visual representation of the protein as a directed graph. Finally, the desired output (optimal frame or final frame) is also saved as an image.

### Output example
Command:
```
python fold_protein.py -f data/sample1.fasta -n 3000 --display-graph optimal --display-final-frame --initial-conformation random
```
Output:
```
The HP model of the protein is: HPHPPHHPHPPHPHHPPHPH
Initial conformation of the protein is randomized with a starting energy of -1.
Find the graph representation in initial_randomized_protein.png
Progress: 100%|████████████████████████████████████████████████████████████████████████████████| 3000/3000 [00:02<00:00, 1115.26it/s]
The final frame can be traced as ULULDLULDDRDRDLLLUR
Graph representation is saved to sample1_optimal_display.png
```
