import math
import random as rd
class Manipulation:
    def __init__(self):
        # registering all the proteins
        all_frames = []
        temperature = 5 # placement provisoire



    def apply_move(self, move):
        return 0

    def undo_move(self, move):
        del self.all_frames[-1]
        return ()

    def test_movement(self, res, move):
        self.apply_movement(move)
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

