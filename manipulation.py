import math
import random
import random as rd
from residue import Residue


class Manipulation:
    temperature = 5
    nb_iterations = 1

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
        for idx, conformation in enumerate(self.all_frames):
            print(f"-------- Frame {idx} --------")
            print(conformation.grid_show())
            print(f"-------------------------")

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
                    new_i, new_j = neighbor
                    res.set_coordinates(new_i, new_j)
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
        if res.index == 0 or res.index == len(conformation.all_residues) - 1:
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

        # case of no available empty positions
        if len(options) == 0:
            return False

        probability = rd.random()
        # if the probability is favorable or if we have only one available option
        if probability < 0.5 or len(options) == 1:
            new_i, new_j = options[0]
            res.set_coordinates(new_i, new_j)
            self.add_frame(conformation)
            return True
        else:
            new_i, new_j = options[1]
            res.set_coordinates(new_i, new_j)
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
                # print("TRUE CRANKSHAFT")
                # checking if the 180 flip positions are available
                if conformation.is_free(res.coordI, res.coordJ + 2) and conformation.is_free(next_residue.coordI,
                                                                                             next_residue.coordJ + 2):
                    res.set_coordinates(res.coordI, res.coordJ + 2)
                    next_residue.set_coordinates(next_residue.coordI, next_residue.coordJ+2)
                    self.add_frame(conformation)
                    print("YAY CRANKSHAFT")
                    return True
            # vertical case
            elif res.coordI == next_residue.coordI:
                if conformation.is_free(res.coordI + 2, res.coordJ) and conformation.is_free(next_residue.coordI + 2,
                                                                                             next_residue.coordJ):
                    res.set_coordinates(res.coordI + 2, res.coordJ)
                    next_residue.set_coordinates(next_residue.coordI + 2, next_residue.coordJ)
                    self.add_frame(conformation)
                    print("YAY CRANKSHAFT")
                    return True
        return False

    def pull_move(self, _conformation, _residue):
        conformation = _conformation.copy_protein()
        prev_residue = conformation.get_previous_residue(_residue)
        residue = conformation.get_next_residue(prev_residue)
        residue_i, residue_j = residue.get_coordinates()
        # to be generalized
        c_positions = (residue_i, residue_j + 1)
        l_positions = (residue_i + 1, residue_j + 1)

        # if i-1 is in C aka a corner move
        if c_positions == prev_residue.get_coordinates():
            residue.set_coordinates(*l_positions)
            self.add_frame(conformation)
            return True

            # C is not empty
        if not conformation.is_free(*c_positions):
            return False

        prev_residue.set_coordinates(*c_positions)
        residue.set_coordinates(*l_positions)

        residue_i_minus_2 = conformation.get_previous_residue(prev_residue)
        # conformation is considered valid if i-2 is next to C which is now occupied by i-1
        if conformation.are_adjacent(residue_i_minus_2.get_coordinates(), c_positions):
            self.add_frame(conformation)
            return True

        j = residue.index - 2
        curr_conformation = conformation
        while j >= 0:
            curr_conformation = curr_conformation.copy_protein()
            curr_residue = curr_conformation.get_residue_at_idx(j)
            # the coordinates that have just been vacated ?
            updated_coordinates = _conformation.get_residue_at_idx(j + 2).get_coordinates()
            curr_residue.set_coordinates(*updated_coordinates)
            # break early if valid conformation is found
            residue_i_minus_2 = curr_conformation.get_residue_at_idx(residue.index - 2)
            if curr_conformation.are_adjacent(residue_i_minus_2.get_coordinates(), c_positions) and \
                    curr_conformation.is_valid():
                break
            j -= 1
        self.add_frame(curr_conformation)
        return True

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