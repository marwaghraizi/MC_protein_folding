import math
import random
import random as rd


class Manipulation:

    def __init__(self, n_iterations, temperature):
        self.all_frames = []
        self.frame_caused_by_move = []
        self.n_iterations = n_iterations
        self.temperature = temperature
        self.all_energies = []
        self.optimal_frame = None

    def add_frame(self, conformation, cause='UNKNOWN'):
        self.all_frames.append(conformation)
        self.frame_caused_by_move.append(cause)
        # class attribute maybe?
        self.protein_length = int(len(conformation.all_residues))

    def show_all_frames(self):
        for conformation in self.all_frames:
            print(conformation.show())
        for idx, conformation in enumerate(self.all_frames):
            print(f"-------- Frame {idx} - Reason {self.frame_caused_by_move[idx]} --------")
            print(conformation.grid_show())
            print(f"-------------------------")

    @staticmethod
    def choose_random_amino_acid(self, protein):
        random_aa = rd.choice(protein.all_residues)
        return random_aa

    def choose_random_move(self, new_protein, residue, search_space="VSHD-pull", pull_proba=0.5):
        #new_protein = self.all_frames[-1].copy_protein()
        #residue = new_protein.all_residues[1]
        #residue = self.choose_random_amino_acid()

        if search_space == "VSHD":
            if residue.index == 0 or residue.index == len(self.all_frames[-1].all_residues)-1:
                decision = self.end_move(new_protein, residue)
            else:
                choice = random.choice([self.crankshaft_move, self.corner_move])
                decision = choice(new_protein, residue)

        if search_space == "VSHD-pull":

            if residue.index == 0 or residue.index == len(self.all_frames[-1].all_residues)-1:
                decision = self.end_move(new_protein, residue)
            else:
                probability = random.random()
                if probability < pull_proba:
                    choice = self.pull_move(new_protein, residue)
                else:
                    choice = random.choice([self.crankshaft_move, self.corner_move])
                decision = choice(new_protein, residue)

        if search_space == "pull":
            decision = self.pull_move(new_protein, residue)

        return decision

    # if the move was not energetically favorable and the probability was low --> need to remove it
    def undo_move(self):
        del self.all_frames[-1]
        del self.frame_caused_by_move[-1]

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
        print()
        # if i-1 and i+1 share an available position --> move i to this position
        print("corner is happening")
        print(conformation.coordinates)
        print("residue of interest:")
        print(res)
        print("its index is:" + str(res.index))
        print(f"its coordinates are: {res.get_coordinates()}")
        print(conformation.show())

        previous_residue = conformation.get_previous_residue(res)
        next_residue = conformation.get_next_residue(res)
        previous_residue_neighbors = conformation.get_empty_topological_positions(previous_residue)
        next_residue_neighbors = conformation.get_empty_topological_positions(next_residue)

        if previous_residue_neighbors:
            for neighbor in previous_residue_neighbors:
                if neighbor in next_residue_neighbors:
                    print(neighbor)
                    new_i, new_j = neighbor
                    conformation.set_coordinates(res, new_i, new_j)
                    self.add_frame(conformation, 'corner move')
                    print(conformation.show())
                    return True
            print("corner does not work here")
            return False

    def end_move(self, conformation, res):
        print()
        print("end is happening")
        # ensure that its an extremity (redundant)
        # get neighbor (second or before last residue)
        # get free topological neighbors of the neighbor
        # there should be two --> randomly choose one
        # for the end move we have two possible positions --> random choice
        # print(conformation.show())
        # print(conformation.coordinates)
        # print("index " + str(res.index))
        if res.index == 0 or res.index == len(conformation.all_residues) - 1:
            pass
        else:
            return False

        if res.index == 0:
            neighbor = conformation.all_residues[1]
        else:
            neighbor = conformation.all_residues[-2]
        #print("neighbor:")
        #print(neighbor.get_coordinates())

        topological_positions = conformation.get_empty_topological_positions(neighbor)
        options = topological_positions
        #print("options are:")
        #print(options)

        # choose topological positions of the neighbor that are not occupied
        #for position in topological_positions:
            #if position not in conformation.coordinates:
                #options.append(position)

        # case of no available empty positions
        if len(options) == 0:
            return False

        probability = rd.random()
        # if the probability is favorable or if we have only one available option
        if probability < 0.5 or len(options) == 1:
            new_i, new_j = options[0]
            #print(new_i, new_j)
            conformation.set_coordinates(res, new_i, new_j)
            self.add_frame(conformation, 'end move w/ 1 option or probability')
            #print("end")
            #print(conformation.coordinates)
            #print(conformation.show())
            print("end move over")
            return True
        else:
            new_i, new_j = options[1]
            #print(new_i, new_j)
            conformation.set_coordinates(res, new_i, new_j)
            self.add_frame(conformation, 'end move w/ 2 option')
            #print(conformation.coordinates)
            #print("end")
            #print(conformation.show())
            print("end move over")
            return True

    def crankshaft_move(self, conformation, res):
        print()
        if len(conformation.all_residues) < 4:
            return False
        print("crankshaft is happening")
        print("residue of interest:")
        print(res)
        print("its index is:" + str(res.index))
        print(f"its coordinates are: {res.get_coordinates()}")
        print(conformation.show())
        # checking the U shape
        next_residue = conformation.get_next_residue(res)
        if conformation.is_right_angle(res) and conformation.is_right_angle(next_residue):
            # horizontal case
            if res.coordJ == next_residue.coordJ:
                # print("TRUE CRANKSHAFT")
                # checking if the 180 flip positions are available
                if conformation.is_free(res.coordI, res.coordJ + 2) and conformation.is_free(next_residue.coordI,
                                                                                             next_residue.coordJ + 2):
                    residue_i_plus_2 = conformation.get_next_residue(next_residue)
                    prev_residue = conformation.get_previous_residue(res)
                    if conformation.are_adjacent(residue_i_plus_2.get_coordinates(), prev_residue.get_coordinates()):
                        print(f"i+2 coordinates {residue_i_plus_2.get_coordinates()}")
                        print(f"i-1 coordinates {prev_residue.get_coordinates()}")
                        # flip up
                        if res.coordJ < prev_residue.coordJ:
                            conformation.set_coordinates(res, res.coordI, res.coordJ + 2)
                            conformation.set_coordinates(next_residue, next_residue.coordI, next_residue.coordJ + 2)
                        else:
                            # flip down
                            conformation.set_coordinates(res, res.coordI, res.coordJ - 2)
                            conformation.set_coordinates(next_residue, next_residue.coordI, next_residue.coordJ - 2)
                        self.add_frame(conformation, 'crankshaft horizontal')
                        print(conformation.coordinates)
                        for res in conformation.all_residues:
                            print(f"residue {res.index} {res} with coordinates {res.get_coordinates()}")
                        print(conformation.show())
                        return True
            # vertical case
            elif res.coordI == next_residue.coordI:
                if conformation.is_free(res.coordI + 2, res.coordJ) and conformation.is_free(next_residue.coordI + 2,
                                                                                             next_residue.coordJ):
                    residue_i_plus_2 = conformation.get_next_residue(next_residue)
                    prev_residue = conformation.get_previous_residue(res)
                    if conformation.are_adjacent(residue_i_plus_2.get_coordinates(), prev_residue.get_coordinates()):
                        # flip right
                        if res.coordI < prev_residue.coordI:
                            conformation.set_coordinates(res, res.coordI + 2, res.coordJ)
                            conformation.set_coordinates(next_residue, next_residue.coordI + 2, next_residue.coordJ + 2)
                        else:
                            # flip left
                            conformation.set_coordinates(res, res.coordI - 2, res.coordJ)
                            conformation.set_coordinates(next_residue, next_residue.coordI - 2, next_residue.coordJ)
                        self.add_frame(conformation, 'crankshaft vertical')
                        print(conformation.coordinates)
                        for res in conformation.all_residues:
                            print(f"residue {res.index} {res} with coordinates {res.get_coordinates()}")
                        print(conformation.show())
                        return True
        print("crankshaft does not work here")
        return False

    @staticmethod
    def get_adjacent(position):
        x, y = position
        return (x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)

    def pull_move(self, _conformation, _residue):
        #print("pull is happening")
        #print(conformation.coordinates)
        #print("residue of interest:")
        #print(residue)
        #print("its index is:" + str(residue.index))
        #print(f"its coordinates are: {residue.get_coordinates()}")
        #print(conformation.show())
        if _residue.index == 0 or _residue.index == len(self.all_frames[-1].all_residues) - 1:
            return False
        conformation = _conformation.copy_protein()
        prev_residue = conformation.get_previous_residue(_residue)
        residue = conformation.get_next_residue(prev_residue)
        next_residue = conformation.get_next_residue(residue)

        residue_coords = residue.get_coordinates()
        residue_i, residue_j = residue_coords

        # L: empty lattice position which is adjacent to i+1 and diagonally adjacent to i
        l_position_options = []
        for position in conformation.get_empty_topological_positions(next_residue):
            if conformation.are_diagonally_adjacent(position, residue_coords):
                l_position_options.append(position)

        # C mutually adjacent to L and i
        c_positions = None
        if l_position_options:
            for position in l_position_options:
                neighbors_of_l = self.get_adjacent(position)
                for neighbor in neighbors_of_l:
                    if neighbor in residue.get_topological_neighbors_positions():
                        if neighbor in conformation.coordinates:
                            if neighbor == prev_residue.get_coordinates():
                                l_positions = position
                                c_positions = neighbor
                        else:
                            l_positions = position
                            c_positions = neighbor
        else:
            # if no L positions are present?
            return False

        #print(l_positions)
        #print(c_positions)

        if not c_positions:
            return False

        #c_positions = (residue_i, residue_j + 1)
        #l_positions = (residue_i + 1, residue_j + 1)

        # if i-1 is in C aka a corner move
        if c_positions == prev_residue.get_coordinates():
            conformation.set_coordinates(residue, *l_positions)
            self.add_frame(conformation, 'pull w/ i-1 in C')
            #print(conformation.show())
            return True

        # C is not empty
        if not conformation.is_free(*c_positions):
            return False

        conformation.set_coordinates(prev_residue, *c_positions)
        conformation.set_coordinates(residue, *l_positions)

        residue_i_minus_2 = conformation.get_previous_residue(prev_residue)
        # conformation is considered valid if i-2 is next to C which is now occupied by i-1
        if conformation.are_adjacent(residue_i_minus_2.get_coordinates(), c_positions):
            self.add_frame(conformation,  'pull w/ 1 step')
            #print(conformation.show())
            return True

        j = residue.index - 2
        curr_conformation = conformation
        while j >= 0:
            curr_conformation = curr_conformation.copy_protein()
            curr_residue = curr_conformation.get_residue_at_idx(j)
            # the coordinates that have just been vacated ?
            updated_coordinates = _conformation.get_residue_at_idx(j + 2).get_coordinates()
            curr_conformation.set_coordinates(curr_residue, *updated_coordinates)
            # break early if valid conformation is found
            residue_i_minus_2 = curr_conformation.get_residue_at_idx(residue.index - 2)
            if curr_conformation.are_adjacent(residue_i_minus_2.get_coordinates(), c_positions) and \
                    curr_conformation.is_valid():
                break
            j -= 1
        self.add_frame(curr_conformation, 'pull w/ n steps')
        print(conformation.show())
        return True

    def apply_monte_carlo(self, search_space="VSHD-pull"):
        for i in range(self.n_iterations):
            print()
            print("-----------ITERATION------------: " + str(i))
            move_successful = False
            amino_acids_used = []
            while move_successful is False and len(amino_acids_used) < self.protein_length:
                new_protein = self.all_frames[-1].copy_protein()
                aa = self.choose_random_amino_acid(new_protein)
                if aa not in amino_acids_used:
                    amino_acids_used.append(aa)
                    move_successful = self.choose_random_move(new_protein, aa, search_space)
            self.test_movement()
