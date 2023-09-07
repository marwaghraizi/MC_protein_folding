import math
import random
import random as rd


class Manipulation:
    temperature = 5

    def __init__(self):
        # registering all the proteins
        self.all_frames = []

    """
    choose a random amino acid --> choose a random move --> test if it can do the move --> if not choose another move 
    --> if all moves are impossible --> choose another aa if not apply the move --> test the energy --> if favorable -->
    keep it --> if not calculate proba --> if proba does not allow it --> remove from frame --> if not keep it
    
    Note: may need to create a test function for the possibility of each move but it causes redundancy)
    """
    def add_frame(self, conformation):
        self.all_frames.append(conformation)

    # maybe in protein
    def choose_random_amino_acid(self, conformation):
        random_aa = rd.choice(conformation.all_residues)
        return random_aa

    @staticmethod
    def choose_random_move(self, conformation, residue):

        if residue.index == 0 or residue.index == len(conformation.ALL_RESIDUES)-1:
            # can extremities do other moves?
            return "end"

        else:
            # can the choices be functions?
            choice = random.choice(["corner", "crankshaft", "pull"])
            if choice == "corner":
                self.corner_move(conformation)

            # test if the moves can be applied
            # what to do if none of them work? gotta choose another amino acid (recursion?)
        return False

    def apply_move(self, conformation):
        return 0

    # if the move was not energetically favorable and the probability was low --> need to remove it
    def undo_move(self, move):
        del self.all_frames[-1]
        return ()

    def test_movement(self, move):
        self.apply_move(move)
        current_energy = self.all_frames[-1].calculate_energy()
        previous_energy = self.all_frames[-2].energy
        if current_energy < previous_energy:
            print("Move was successful")
            return()
        else:
            energy_difference = current_energy - previous_energy
            probability = math.exp(energy_difference/self.temperature)
            random_number = rd.random()
            if random_number < probability:
                print("Unfavorable move was accepted")
                return()
            else:
                print("Move was not accepted")
                self.undo_move(move)
                return()

    # apply the move and return true if it was successful
    # check if a position is empty

    # moves (return true or false)
    @staticmethod
    def corner_move(self, conformation, res):
        # if i-1 and i+1 share an available position --> move i to this position
        previous_residue = conformation.get_previous_residue(res)
        next_residue = conformation.get_next_residue(res)
        previous_residue_neighbors = previous_residue.get_topological_neighbors()
        next_residue_neighbors = next_residue.get_topological_neighbors()

        for neighbor in previous_residue_neighbors:
            if neighbor in next_residue_neighbors and neighbor not in conformation.coordinates:
                res.set_coordinates(*neighbor)
                return True
        return False

    @staticmethod
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

        topological_positions = conformation.get_topological_positions(neighbor)
        options = []

        # choose topological positions of the neighbor that are not occupied
        for position in topological_positions:
            if position not in conformation.coordinates:
                options.add(position)

        probability = rd.random()
        if probability < 0.5:
            res.set_coordinates(options[0])
            return True
        else:
            res.set_coordinates(options[1])
            return True
        return False

    @staticmethod
    def crankshaft_move(self, res):

        return False

    @staticmethod
    def choose_random_move(self):
        return False
