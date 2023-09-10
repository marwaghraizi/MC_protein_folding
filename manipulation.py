import math
import random
import random as rd
from residue import Residue


class Manipulation:
    temperature = 5
    nb_iterations = 10000

    def __init__(self):
        # registering all the proteins
        self.all_frames = []

    def add_frame(self, conformation):
        self.all_frames.append(conformation)
        # class attribute maybe?
        self.protein_length = int(len(conformation.all_residues))

    def show_all_frames(self):
        for conformation in self.all_frames:
            print(conformation.show())

    # maybe in protein
    def choose_random_amino_acid(self):
        random_aa = rd.choice(self.all_frames[-1].all_residues)
        return random_aa

    def choose_random_move(self, residue):
        new_protein = self.all_frames[-1].copy_protein()
        #residue = new_protein.all_residues[1]
        #residue = self.choose_random_amino_acid()
        decision = False
        if residue.index == 0 or residue.index == len(self.all_frames[-1].all_residues)-1:
            # can extremities do other moves?
            print("end")
            decision = self.end_move(new_protein, residue)
        else:
            # can the choices be functions?
            choice = random.choice([self.crankshaft_move, self.corner_move])
            decision = choice(new_protein, residue)


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
        if previous_residue_neighbors:
            for neighbor in previous_residue_neighbors:
                if neighbor in next_residue_neighbors:
                    #print("YAY")
                    newI, newJ = neighbor
                    res.set_coordinates(newI, newJ)
                    self.add_frame(conformation)
                    return True
            #print("did not work")
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
        if len(conformation.all_residues) < 4:
            return False

        # checking the U shape
        next_residue = conformation.get_next_residue(res)
        if conformation.is_right_angle(res) and conformation.is_right_angle(next_residue):
            # horizontal case
            if res.coordJ == next_residue.coordJ:
                print("TRUE CRANKSHAFT")
                if conformation.is_free(res.coordI, res.coordJ+2) and conformation.is_free(next_residue.coordI,
                                                                                           next_residue.coordJ+2):
                    res.set_coordinates(res.coordI, res.coordJ+2)
                    next_residue.set_coordinates(next_residue.coordI, next_residue.coordJ+2)
                    self.add_frame(conformation)
                    print("YAY CRANKSHAFT")
                    return True
            # vertical case
            elif res.coordI == next_residue.coordI:
                # checking if the 180 flip positions are available
                if conformation.is_free(res.coordI+2, res.coordJ) and conformation.is_free(next_residue.coordI+2,
                                                                                           next_residue.coordJ):
                    res.set_coordinates(res.coordI+2, res.coordJ)
                    next_residue.set_coordinates(next_residue.coordI+2, next_residue.coordJ)
                    self.add_frame(conformation)
                    print("YAY")
                    return True
        return False

    def pull_move(self, conformation, residue):
        # i-1
        previous_residue = conformation.get_previous_residue(residue)
        # i-2
        previous_previous_residue = conformation.get_previous_residue(previous_residue)
        # i+1
        next_residue = conformation.get_next_residue(residue)
        # checking the availability of the C position (same as corner)
        previous_previous_residue_neighbors = conformation.get_empty_topological_positions(previous_previous_residue)
        residue_neighbors = conformation.get_empty_topological_positions(residue)
        decision = False
        C_neighbors = []
        if previous_previous_residue_neighbors:
            for neighbor in previous_previous_residue_neighbors:
                if neighbor in residue_neighbors:
                    decision = True
                    C_coordI, C_coordJ = neighbor
                    temp_residue_C = Residue("H", -1)
                    temp_residue_C.set_coordinates(C_coordI, C_coordJ)
                    C_neighbors = conformation.get_empty_topological_positions(temp_residue_C)

        # checking the availability of the L position (same as corner)
        next_residue_neighbors = conformation.get_empty_topological_positions(next_residue)

        if C_neighbors:
            for neighbor in C_neighbors:
                if neighbor in next_residue_neighbors:
                    decision = True
                    L_coordI, L_coordJ = neighbor
                    residue.set_coordinates(C_coordI, C_coordJ)
                    next_residue.set_coordinates(L_coordI, L_coordJ)
                    self.add_frame(conformation)
                    return decision
        return decision

    def apply_monte_carlo(self):
        for i in range(self.nb_iterations):
            print(i)
            move_successful = False
            amino_acids_used = []
            while move_successful is False and len(amino_acids_used) < self.protein_length:
                aa = self.choose_random_amino_acid()
                if aa not in amino_acids_used:
                    amino_acids_used.append(aa)
                    move_successful = self.choose_random_move(aa)
            self.test_movement()

    def apply_MC_VSHD(self):
        return False

    def apply_MC_VSHD_pull(self):
        return False

    def apply_MC_pull(self):
        return False