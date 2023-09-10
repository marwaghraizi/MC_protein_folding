import random
import argparse
from protein import Protein
from residue import Residue
from manipulation import Manipulation

random.seed(5)
polar_residues = ["E", "D", "H", "T", "S", "Y", "N", "Q", "R", "K", "H"]
hydrophobic_residues = ["C", "W", "G", "A", "P", "I", "L", "M", "F", "V"]
nb_iterations = 1000

parser = argparse.ArgumentParser(description="Fold Protein")
parser.add_argument('-f', '--file', metavar='protein', type=str, help='Protein File Path')
args = parser.parse_args()

protein_file = args.file
with open(protein_file, "r") as filein:
    sequence = ""
    for line in filein:
        if not line.startswith(">"):
            sequence += line.strip()

HP_sequence = ""
ALL_RESIDUES = []
# set of tuples for hydrophobic contacts

for idx, residue in enumerate(sequence):
    if residue in polar_residues:
        HP_sequence += "P"
        residue_object = Residue("P", idx)
        residue_object.set_coordinates(idx, 0)
    else:
        HP_sequence += "H"
        residue_object = Residue("H", idx)
        residue_object.set_coordinates(idx, 0)
    ALL_RESIDUES.append(residue_object)

initial_protein = Protein(ALL_RESIDUES)
# testing corner
#initial_protein.all_residues[0].coordI = 1
#initial_protein.all_residues[0].coordJ = 1

# testing crankshaft
initial_protein.all_residues[1].coordI = 0
initial_protein.all_residues[1].coordJ = -1
initial_protein.all_residues[2].coordI = 1
initial_protein.all_residues[2].coordJ = -1
initial_protein.all_residues[3].coordI = 1
initial_protein.all_residues[3].coordJ = 0


print(initial_protein.show())

manipulation = Manipulation()
manipulation.add_frame(initial_protein)
manipulation.apply_monte_carlo()
manipulation.show_all_frames()


""" 
occupied_positions = [residue.get_coordinates() for residue in ALL_RESIDUES]
print(occupied_positions)


def is_free(a,b):
    if (a, b) in occupied_positions:
        return False
    else:
        return True


def get_four_neighbors(a,b):
    return (a-1, b), (a+1, b), (a, b-1), (a, b+1)


for i in range(len(ALL_RESIDUES)-1):
    neighbors = get_four_neighbors(*(ALL_RESIDUES[i].get_coordinates()))
    random_neighbor = rd.choice(neighbors)
    isFree = False

    if is_free(*(random_neighbor)):
        isFree = True

    while(not isFree):
        random_neighbor = rd.choice(neighbors)
        if is_free(*(random_neighbor)):
            isFree = True

    ALL_RESIDUES[i+1].set_coordinates(*(random_neighbor))
    occupied_positions = [residue.get_coordinates() for residue in ALL_RESIDUES]







"""
