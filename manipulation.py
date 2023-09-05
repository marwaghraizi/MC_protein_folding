import math
import random
import random as rd
class Manipulation:
    def __init__(self):
        # registering all the proteins
        all_frames = []
        temperature = 5 # placement provisoire

    def choose_random_amino_acid(self, conformation):
        return()

    def choose_random_move(self,conformation, residue):
        if residue.index == 0 or residue.index == len(conformation.ALL_RESIDUES)-1:
            return "end"
        if residue.index > 0 and residue.index < len(conformation.ALL_RESIDUES)-1:
            choice = random.choice(["corner", "crankshaft", "pull"])
            # test if the moves can be applied
            # what to do if none of them work? gotta choose another amino acid (recursion?)

    def apply_move(self, move):
        return 0

    # if the move was not energetically favorable and the probability was low --> need to remove it
    def undo_move(self, move):
        del self.all_frames[-1]
        return ()

    def test_movement(self, res, move):
        self.apply_move(move)
        current_energy = self.all_frames[-1].calculate_energy()
        previous_energy = self.all_frames[-2].enery
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

