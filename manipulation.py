import math
import random
from tqdm import tqdm


class Manipulation:
    """
        A class to apply the Monte Carlo protein folding algorithm.

        ...

        Attributes
        ----------
        all_frames : list
            list of all the protein objects forming the protein.
        moves : list
            list of applied moves in each iteration to create a new conformation
        n_iterations : int
            energy of the protein.
        temperature : float
            search temperature.

        Methods
        -------
        add_frame(conformation, cause='UNKNOWN'):
            Adds the new frame (protein conformation) and the move causing it
            to the corresponding attributes.
        get_adjacent(position):
            Returns tuples of coordinates adjacent to a given position.
        choose_random_amino_acid(protein):
            Returns a random amino acid in a given protein.
        choose_random_move(new_protein, residue, search_space="VSHD-pull",
                            pull_probability=0.5)
            Chooses a random move and applies it if possible.
        end_move(conformation, res):
            Applies the end move if possible and adds a new frame.
        corner_move(conformation, res):
            Applies the corner move if possible and adds a new frame.
        crankshaft_move(conformation, res):
            Applies the crankshaft move if possible and adds a new frame.
        pull_move(conformation, res)
            Applies the pull move if possible and adds a new frame.
        undo_move()
            Deletes the last frame
        test_movement()
            Tests the energy of last protein added to the frames and removes it
            if it is unfavorable with a probability
        apply_monte_carlo(search_space="VSHD-pull", pull_probability=0.5)
            Applies the Monte Carlo search for n iterations.
        get_optimal_frame()
            Returns the protein with the lowest energy.
        show_all_frames()
            Prints all frames and the move causing it in the direction text
            visualization.
        """

    def __init__(self, n_iterations, temperature):
        """
        Constructs all the necessary attributes for the manipulation object.

        Parameters
        ----------
            n_iterations : int
                Number of search iterations.
            temperature : float
                Search temperature
        """
        self.n_iterations = n_iterations
        self.temperature = temperature
        self.all_frames = []
        self.moves = []
        self.all_energies = []
        self.optimal_frame = None

    def add_frame(self, conformation, cause='UNKNOWN'):
        """Adds the new frame (protein conformation) and the move causing it to
        the corresponding attributes."""
        self.all_frames.append(conformation)
        self.moves.append(cause)

    @staticmethod
    def get_adjacent(position):
        """Returns tuples of coordinates adjacent to a given position."""
        x, y = position
        return (x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)

    @staticmethod
    def choose_random_amino_acid(protein):
        """Returns a random amino acid in a given protein."""
        random_aa = random.choice(protein.all_residues)
        return random_aa

    def choose_random_move(self, new_protein, residue, search_space="VSHD-pull",
                           pull_probability=0.5):
        """
        Choose a random move to execute and returns a boolean.

        Parameters
        ----------
        new_protein : Protein
            Protein to apply the move on.
        residue : Residue
            Residue of interest to be moved.
        search_space : str
            Search neighborhood: VSHD-pull, VSHD or pull
        pull_probability : float
            Probability of choosing the pull move in the hybrid VSHD-pull search
            neighborhood

        Returns
        -------
        bool
            Indicating the success or failure of the randomly chosen move.
        """
        decision = False
        if search_space == "VSHD":
            if residue.index == 0 or \
                    residue.index == len(self.all_frames[-1].all_residues)-1:
                decision = self.end_move(new_protein, residue)
            else:
                choice = random.choice([self.crankshaft_move, self.corner_move])
                decision = choice(new_protein, residue)

        if search_space == "VSHD-pull":

            if residue.index == 0 or \
                    residue.index == len(self.all_frames[-1].all_residues)-1:
                decision = self.end_move(new_protein, residue)
            else:
                probability = random.random()
                if probability < pull_probability:
                    decision = self.pull_move(new_protein, residue)
                else:
                    choice = random.choice([self.crankshaft_move,
                                            self.corner_move])
                    decision = choice(new_protein, residue)

        if search_space == "pull":
            decision = self.pull_move(new_protein, residue)

        return decision

    def end_move(self, conformation, res):
        """Applies the end move of a given residue in a given conformation. If it is possible, the new conformation
        is added to the frames and a true is returned else a false is returned. """
        # checking that the residue is an extremity
        if res.index == 0 or res.index == len(conformation.all_residues) - 1:
            pass
        else:
            return False

        # extracting corresponding connected residue
        if res.index == 0:
            neighbor = conformation.all_residues[1]
        else:
            neighbor = conformation.all_residues[-2]

        # extracting the available neighboring positions of the connected residue
        options = conformation.get_empty_topological_positions(neighbor)

        # case of no available empty positions
        if len(options) == 0:
            return False

        probability = random.random()
        # if the probability is favorable or if we have only one available option
        if probability < 0.5 or len(options) == 1:
            new_i, new_j = options[0]
            conformation.set_coordinates(res, new_i, new_j)
            self.add_frame(conformation, 'end move w/ 1 option or probability')
            return True
        else:
            new_i, new_j = options[1]
            conformation.set_coordinates(res, new_i, new_j)
            self.add_frame(conformation, 'end move w/ 2 option')
            return True

    def corner_move(self, conformation, res):
        """Applies the corner move of a given residue in a given conformation. If it is possible, the new conformation
        is added to the frames and a true is returned else a false is returned. """
        # extracting connecting residues and their empty neighboring positions
        previous_residue = conformation.get_previous_residue(res)
        next_residue = conformation.get_next_residue(res)
        previous_residue_neighbors = conformation.get_empty_topological_positions(previous_residue)
        next_residue_neighbors = conformation.get_empty_topological_positions(next_residue)

        # if the previous (i-1) residue has empty neighbors
        if previous_residue_neighbors:
            for neighbor in previous_residue_neighbors:
                # check if it is also the neighbor of the next (i+1) residue
                if neighbor in next_residue_neighbors:
                    new_i, new_j = neighbor
                    conformation.set_coordinates(res, new_i, new_j)
                    self.add_frame(conformation, 'corner move')
                    return True
        return False

    def crankshaft_move(self, conformation, res):
        """Applies the crankshaft move of a given residue in a given
        conformation. If it is possible, the new conformation is added
        to the frames and a true is returned else a false is returned. """
        # Crankshaft requires a U shape of at least 4 amino acids
        if len(conformation.all_residues) < 4:
            return False

        next_residue = conformation.get_next_residue(res)
        # checking if residues i and i+1 constitute right angles
        if conformation.is_right_angle(res) and \
                conformation.is_right_angle(next_residue):
            # horizontal case
            if res.coordJ == next_residue.coordJ:
                residue_i_plus_2 = conformation.get_next_residue(next_residue)
                prev_residue = conformation.get_previous_residue(res)
                # checking if i+2 and i-1 are adjacent to avoid Z shapes
                if conformation.are_adjacent(residue_i_plus_2.get_coordinates(),
                                             prev_residue.get_coordinates()):
                    # flip up
                    if res.coordJ < prev_residue.coordJ:
                        # checking if the 180 flip positions are available
                        if conformation.is_free(res.coordI, res.coordJ + 2) \
                                and conformation.is_free(next_residue.coordI,
                                                         next_residue.coordJ + 2):
                            conformation.set_coordinates(res, res.coordI,
                                                         res.coordJ + 2)
                            conformation.set_coordinates(next_residue,
                                                         next_residue.coordI,
                                                         next_residue.coordJ + 2)
                            self.add_frame(conformation, 'crankshaft horizontal')
                            return True
                    # flip down
                    else:
                        if conformation.is_free(res.coordI, res.coordJ - 2) \
                                and conformation.is_free(next_residue.coordI,
                                                         next_residue.coordJ - 2):
                            conformation.set_coordinates(res, res.coordI,
                                                         res.coordJ - 2)
                            conformation.set_coordinates(next_residue,
                                                         next_residue.coordI,
                                                         next_residue.coordJ - 2)
                            self.add_frame(conformation, 'crankshaft horizontal')
                            return True
            # vertical case
            elif res.coordI == next_residue.coordI:
                residue_i_plus_2 = conformation.get_next_residue(next_residue)
                prev_residue = conformation.get_previous_residue(res)
                if conformation.are_adjacent(residue_i_plus_2.get_coordinates(),
                                             prev_residue.get_coordinates()):
                    # flip right
                    if res.coordI < prev_residue.coordI:
                        if conformation.is_free(res.coordI + 2, res.coordJ) and \
                                conformation.is_free(next_residue.coordI + 2,
                                                     next_residue.coordJ):
                            conformation.set_coordinates(res, res.coordI + 2,
                                                         res.coordJ)
                            conformation.set_coordinates(next_residue,
                                                         next_residue.coordI + 2,
                                                         next_residue.coordJ)
                            self.add_frame(conformation, 'crankshaft vertical')
                            return True
                    # flip left
                    else:
                        if conformation.is_free(res.coordI - 2, res.coordJ) and \
                                conformation.is_free(next_residue.coordI - 2,
                                                     next_residue.coordJ):
                            conformation.set_coordinates(res, res.coordI - 2,
                                                         res.coordJ)
                            conformation.set_coordinates(next_residue,
                                                         next_residue.coordI - 2,
                                                         next_residue.coordJ)
                            self.add_frame(conformation, 'crankshaft vertical')
        return False

    def pull_move(self, _conformation, _residue):
        if _residue.index == 0 \
                or _residue.index == len(self.all_frames[-1].all_residues) - 1:
            return False
        conformation = _conformation.copy_protein()
        prev_residue = conformation.get_previous_residue(_residue)
        residue = conformation.get_next_residue(prev_residue)
        next_residue = conformation.get_next_residue(residue)

        residue_coordinates = residue.get_coordinates()

        # L: empty lattice position which is adjacent to i+1 and
        # diagonally adjacent to i
        l_position_options = []
        for position in conformation.get_empty_topological_positions(next_residue):
            if conformation.are_diagonally_adjacent(position,
                                                    residue_coordinates):
                l_position_options.append(position)

        # C mutually adjacent to L and i
        c_positions = tuple()
        l_positions = tuple()
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
            # if no L positions are present
            return False

        # if c positions do not exist
        if not c_positions:
            return False

        # if i-1 is in C --> a corner move
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
        # conformation is considered valid if i-2 is next to C
        # which is now occupied by i-1
        if conformation.are_adjacent(residue_i_minus_2.get_coordinates(),
                                     c_positions):
            self.add_frame(conformation,  'pull w/ 1 step')
            return True

        j = residue.index - 2
        curr_conformation = conformation
        while j >= 0:
            curr_conformation = curr_conformation.copy_protein()
            curr_residue = curr_conformation.get_residue_at_idx(j)
            # the coordinates that have just been vacated
            updated_coordinates = _conformation.get_residue_at_idx(j + 2).get_coordinates()
            curr_conformation.set_coordinates(curr_residue, *updated_coordinates)
            # break early if valid conformation is found
            residue_i_minus_2 = curr_conformation.get_residue_at_idx(residue.index - 2)
            if curr_conformation.are_adjacent(residue_i_minus_2.get_coordinates(),
                                              c_positions) and curr_conformation.is_valid():
                break
            j -= 1
        self.add_frame(curr_conformation, 'pull w/ n steps')
        return True

    def undo_move(self):
        """Deletes the last added frame and move causing it in the case of it being energetically unfavorable"""
        del self.all_frames[-1]
        del self.moves[-1]

    def test_movement(self):
        """Tests the last protein added to the frames. If its energy increases with respect to the previous one,
        the move is considered successful. Otherwise, a probability is calculated based on the temperature to
        accept or reject an energetically unfavorable move."""
        current_energy = self.all_frames[-1].calculate_energy()
        previous_energy = self.all_frames[-2].calculate_energy()
        # if the added protein has a lower energy than the previous one,
        # frames are kept as is
        if current_energy <= previous_energy:
            #print("Move was successful")
            pass
        else:
            # if the protein has a higher energy than the previous one,
            # it is accepted based on random probability.
            energy_difference = current_energy - previous_energy
            probability = math.exp(energy_difference/self.temperature)
            random_number = random.random()
            if random_number < probability:
                #print("Unfavorable move was accepted")
                pass
            else:
                #print("Move was not accepted")
                self.undo_move()

    def apply_monte_carlo(self, search_space="VSHD-pull", pull_probability=0.5):
        """Applies the search whereby each iteration includes choosing an amino
        acid, a random movement, testing if the amino acid can undergo the
        movement, if not choosing another amino acid until a successful move
        that is subsequently tested for its energy."""
        for _ in tqdm(range(self.n_iterations), desc="Progress"):
            move_successful = False
            # list to keep track of amino acids in the case of none of
            # them being able to apply a move
            amino_acids_used = []
            # if the randomly chosen move fails choose another amino acid
            while move_successful is False \
                    and len(amino_acids_used) < len(self.all_frames[0].all_residues):
                # deep copy the protein before choosing the amino acid
                # and applying the move
                new_protein = self.all_frames[-1].copy_protein()
                aa = self.choose_random_amino_acid(new_protein)
                if aa not in amino_acids_used:
                    amino_acids_used.append(aa)
                    move_successful = self.choose_random_move(new_protein, aa,
                                                              search_space,
                                                              pull_probability)
            # once a move is successful, it has to be tested energy wise
            self.test_movement()

    def get_optimal_frame(self):
        """Returns the protein with the lowest energy."""
        min_energy = 0
        optimal_conformation = None
        for protein in self.all_frames:
            self.all_energies.append(protein.calculate_energy())
            curr_energy = protein.calculate_energy()
            if curr_energy <= min_energy:
                min_energy = curr_energy
                optimal_conformation = protein
        return optimal_conformation

    def show_all_frames(self):
        """Prints all frames and the move causing it in the direction
        text visualization."""
        for idx, conformation in enumerate(self.all_frames):
            print(f"-------- Frame {idx} - Caused by {self.moves[idx]} -------")
            print(conformation.show())
            #print(conformation.grid_show())
