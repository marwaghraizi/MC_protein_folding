import random
import random as rd
import numpy as np
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
initial_protein.all_residues[0].coordI = 1
initial_protein.all_residues[0].coordJ = 1


print(HP_sequence)
print(initial_protein.show())

manipulation = Manipulation()
manipulation.add_frame(initial_protein)
print(initial_protein)
manipulation.apply_monte_carlo()
manipulation.choose_random_move()
manipulation.show_all_frames()
print(len(manipulation.all_frames))


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

print(occupied_positions)

test_protein = Protein(ALL_RESIDUES)
print(test_protein.calculate_energy())

print(test_protein.show())

hydrophobic_contacts = set()

def calculate_energy():
    hydrophobic_objects = [residue for residue in ALL_RESIDUES if residue.type=="H"]
    print(hydrophobic_objects)

hydrophobic_objects = [residue for residue in ALL_RESIDUES if residue.type=="H"]
print(HP_sequence)
print(len(hydrophobic_objects))


for i in range(1,len(hydrophobic_objects)-1):
    # getting the 4 neighboring positions
    neighbors = get_four_neighbors(*(hydrophobic_objects[i].get_coordinates()))
    # getting the preceding and the next residue index
    preceding_neighbor_coord = (ALL_RESIDUES[hydrophobic_objects[i].index-1]).get_coordinates()
    next_neighbor_coord = (ALL_RESIDUES[hydrophobic_objects[i].index+1]).get_coordinates()
    neighbors2 = list(neighbors)
    neighbors2.remove(preceding_neighbor_coord)
    neighbors2.remove(next_neighbor_coord)
    topo_neighbors = tuple(neighbors2)
    print(topo_neighbors)
    coordinates_dictionary = {}
    for neighbor in topo_neighbors:
        # get the type of the residue whose coordinates match the topo neighbor coordinates
        # need to check if its occupied
        # memory speed tradeoff: create a dictionary with coordinates as keys
        neighbor_type = str([x.type for x in ALL_RESIDUES if x.get_coordinates()==neighbor])
        print(neighbor_type)


"""
