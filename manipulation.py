import math
import random
import random as rd
from protein import Protein
from copy import deepcopy


class Manipulation:
    temperature = 5
    nb_iterations = 0

    def __init__(self):
        # registering all the proteins
        self.all_frames = []

    """
    choose a random amino acid --> choose a random move --> test if it can do the move --> if not choose another move 
    --> if all moves are impossible --> choose another aa if not apply the move --> test the energy --> if favorable -->
    keep it --> if not calculate proba --> if proba does not allow it --> remove from frame --> if not keep it
    
    """
    def add_frame(self, conformation):
        self.all_frames.append(conformation)

    def show_all_frames(self):
        for conformation in self.all_frames:
            print(conformation.show())

    # maybe in protein
    def choose_random_amino_acid(self):
        random_aa = rd.choice(self.all_frames[-1].all_residues)
        return random_aa

    def choose_random_move(self):
        new_protein = self.all_frames[-1].copy_protein()
        #residue = new_protein.all_residues[1]
        residue = self.choose_random_amino_acid()
        decision = False
        if residue.index == 0 or residue.index == len(self.all_frames[-1].all_residues)-1:
            # can extremities do other moves?
            print("end")
            decision = self.end_move(new_protein, residue)
        else:
            # can the choices be functions?
            choice = random.choice(["corner"])
            print(choice)
            if choice == "corner":
                decision = self.corner_move(new_protein, residue)
            if choice == "crankshaft":
                decision = self.crankshaft_move(new_protein, residue)

            # test if the moves can be applied
            # what to do if none of them work? gotta choose another amino acid (recursion?)
        return decision

    # if the move was not energetically favorable and the probability was low --> need to remove it
    def undo_move(self):
        del self.all_frames[-1]

    def test_movement(self):
        current_energy = self.all_frames[-1].calculate_energy()
        previous_energy = self.all_frames[-2].calculate_energy()
        if current_energy <= previous_energy:
            print("Move was successful")
        else:
            energy_difference = current_energy - previous_energy
            probability = math.exp(energy_difference/self.temperature)
            random_number = rd.random()
            if random_number < probability:
                print("Unfavorable move was accepted")
            else:
                print("Move was not accepted")
                self.undo_move()
    # apply the move and return true if it was successful
    # check if a position is empty

    # moves (return true or false)
    def corner_move(self, conformation, res):
        # if i-1 and i+1 share an available position --> move i to this position
        previous_residue = conformation.get_previous_residue(res)
        next_residue = conformation.get_next_residue(res)
        previous_residue_neighbors = conformation.get_empty_topological_positions(previous_residue)
        next_residue_neighbors = conformation.get_empty_topological_positions(next_residue)
        # should account for empty lists??
        for neighbor in previous_residue_neighbors:
            if neighbor in next_residue_neighbors:
                print("YAY")
                newI, newJ = neighbor
                res.set_coordinates(newI, newJ)
                self.add_frame(conformation)
                return True
        print("did not work")
        return False

    def end_move(self, conformation, res):
        # ensure that its an extremity (redundant)
        # get neighbor (second or before last residue)
        # get free topological neighbors of the neighbor
        # there should be two --> randomly choose one
        # for the end move we have two possible positions --> random choice
        if res.index == 0 or res.index == len(conformation.all_residues)-1:
            pass
        else:
            return False

        if res.index == 0:
            neighbor = conformation.all_residues[1]
        else:
            neighbor = conformation.all_residues[-2]

        topological_positions = neighbor.get_topological_neighbors_positions()
        options = []

        # choose topological positions of the neighbor that are not occupied
        for position in topological_positions:
            if position not in conformation.coordinates:
                options.append(position)

        probability = rd.random()
        if probability < 0.5:
            newI, newJ = options[0]
            res.set_coordinates(newI, newJ)
            self.add_frame(conformation)
            return True
        else:
            newI, newJ = options[1]
            res.set_coordinates(newI, newJ)
            self.add_frame(conformation)
            return True

    def crankshaft_move(self, conformation, res):

        return False

    def apply_monte_carlo(self):
        for i in range(self.nb_iterations):
            move_successful = False
            # while can be infinite
            while move_successful is False:
                aa = self.choose_random_amino_acid()
                move_successful = self.choose_random_move()
            self.test_movement()
