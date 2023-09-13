from copy import deepcopy


class Protein:
    """
        A class to represent a protein.

        ...

        Attributes
        ----------
        all_residues : list
            list of all the residue objects forming the protein.
        coordinates : dict
            dictionary with a tuple of coordinates as the key and corresponding
            residue as the value.
        energy : int
            energy of the protein.

        Methods
        -------
        set_coordinates(res, new_x, new_y):
            Updates a residues coordinates as well as the conformation
            coordinates and its energy.
        is_free(x, y):
            Checks if a given position is empty.
        get_residue_at_idx(idx):
            Returns Residue object at a given index.
        get_next_residue(res)
            Returns the following residue of a given residue.
        get_previous_residue(res):
            Returns the previous residue of a given residue.
        get_connected_residues(res):
            Returns a list of covalently connected residues to a given residue.
        get_occupied_topological_neighbors(res):
            Returns the residues surrounding a given residue and not covalently
            connected to it.
        get_empty_topological_positions(res):
            Returns the coordinates of the empty positions surrounding a given
            residue.
        is_right_angle(res):
            Checks if a given residue forms a right angle with its previous and
            next residue.
        are_adjacent(first_residue_coordinates, second_residue_coordinates):
            Checks if two residues are adjacent to each other.
        are_diagonally_adjacent(first_residue_coordinates,
                                second_residue_coordinates):
            Checks if two residues are diagonally adjacent to each other.
        is_valid():
            Checks if a conformation is valid if all residues are adjacent to
            their previous residue.
        copy_protein():
            Returns a deep copied version of the protein.
        calculate_energy():
            Returns the overall energy of a protein obtained
            from all hydrophobic contacts.
        describe_direction(x1, y1, x2, y2):
            Returns movement direction from one position to another.
        show():
            Returns a text representation that traces
            a protein with direction symbols.
        grid_show():
            Returns a command line graphical visualisation of a protein.
        graph_show(file_name):
            Creates a .PNG file with the graphical representation
            of the protein.
        """
    def __init__(self, residues):
        """
        Constructs all the necessary attributes for the residue object.

        Parameters
        ----------
            residues : list
                list of Residue objects
        """
        self.all_residues = residues
        self.coordinates = {}
        for residue in residues:
            self.coordinates[(residue.coordI, residue.coordJ)] = residue
        self.energy = self.calculate_energy()

    def set_coordinates(self, residue, new_x, new_y):
        """Updates the coordinate of the residue, the overall conformation
        and the energy given a residue and new coordinates"""
        coordinates = (new_x, new_y)
        old_coordinates = residue.get_coordinates()
        residue.set_coordinates(new_x, new_y)
        del self.coordinates[old_coordinates]
        self.coordinates[coordinates] = residue
        self.energy = self.calculate_energy()

    def is_free(self, x, y):
        """Checks if a given position is empty."""
        if (x, y) in self.coordinates:
            return False
        else:
            return True

    def get_residue_at_idx(self, index):
        """Returns residue object at a given index"""
        return self.all_residues[index]

    def get_next_residue(self, res):
        """Returns the following residue (C' end) of a given residue."""
        return self.all_residues[res.index+1]

    def get_previous_residue(self, res):
        """Returns the previous residue (N' end) of a given residue."""
        return self.all_residues[res.index-1]

    def get_connected_residues(self, res):
        """Returns a list of covalently connected residues to a
        given residue."""
        connected_residues = []
        # first residue
        if res.index == 0:
            connected_residues.append(self.all_residues[res.index+1])
        # last residue
        elif res.index == len(self.all_residues)-1:
            connected_residues.append(self.all_residues[res.index-1])
        # middle residue
        else:
            connected_residues.append(self.all_residues[res.index - 1])
            connected_residues.append(self.all_residues[res.index + 1])

        return connected_residues

    def get_occupied_topological_neighbors(self, res):
        """Returns the residues surrounding a given residue and not covalently
        connected to it

        Parameters
        ----------
        res : Residue object
            residue of interest.

        Returns
        -------
        list
            list of neighboring but not connected residues of a residue.
        """
        # get connected residues
        connected_residues = self.get_connected_residues(res)
        connected_residues_positions = {}

        for connected_residue in connected_residues:
            connected_residues_positions[(connected_residue.coordI, connected_residue.coordJ)] = connected_residue

        topological_positions = res.get_topological_neighbors_positions()
        # topological neighbors are in the topological positions and are
        # not connected to the residue of interest
        occupied_topological_neighbors = []

        for position in topological_positions:
            # if the position is occupied and is not a connected amino acid
            if (position in self.coordinates) and (position not in connected_residues_positions):
                occupied_topological_neighbors.append(self.coordinates[position])

        return occupied_topological_neighbors

    def get_empty_topological_positions(self, res):
        """Returns the coordinates of the empty positions surrounding a given
        residue

        Parameters
        ----------
        res : Residue object
            residue of interest.

        Returns
        -------
        list
            list of tuples of coordinates of empty neighboring
            positions of a residue.
        """
        topological_positions = res.get_topological_neighbors_positions()
        empty_topological_positions = []
        for position in topological_positions:
            if position not in self.coordinates:
                empty_topological_positions.append(position)
        return empty_topological_positions

    def is_right_angle(self, res):
        """Checks if a given residue forms a right angle with its
        previous and next residue."""
        if res.index == 0 or res.index == len(self.all_residues) - 1:
            return False
        x_res, y_res = res.get_coordinates()
        previous_res = self.get_previous_residue(res)
        x_prev, y_prev = previous_res.get_coordinates()
        next_res = self.get_next_residue(res)
        x_next, y_next = next_res.get_coordinates()
        # vectors formed by residue and its previous one as well as
        # residue and its following one
        vector1 = (x_next - x_res, y_next - y_res)
        vector2 = (x_prev - x_res, y_prev - y_res)
        # dot product is the sum of the product of corresponding
        # elements in two vectors
        dot_product = vector1[0] * vector2[0] + vector1[1] * vector2[1]
        # if the dot product == 0 --> vectors are perpendicular
        return dot_product == 0

    @staticmethod
    def are_adjacent(first_residue_coordinates, second_residue_coordinates):
        """Checks if two residues are adjacent to each other."""
        i1, j1 = first_residue_coordinates
        i2, j2 = second_residue_coordinates

        if abs(i1 - i2) == 1 and j1 == j2:
            return True

        if i1 == i2 and abs(j1 - j2) == 1:
            return True

        return False

    @staticmethod
    def are_diagonally_adjacent(first_residue_coordinates, second_residue_coordinates):
        """Checks if two residues are diagonally adjacent to each other."""
        i1, j1 = first_residue_coordinates
        i2, j2 = second_residue_coordinates
        return abs(i1 - i2) == 1 and abs(j1 - j2) == 1

    def is_valid(self):
        """Checks if a conformation is valid if all residues are adjacent
        to their previous residue."""
        idx = 1
        for residue in self.all_residues[1:]:
            if residue.get_coordinates() not in self.all_residues[idx - 1].get_topological_neighbors_positions():
                return False
            idx += 1
        return True

    def copy_protein(self):
        """Returns a deep copied version of the protein."""
        copied_residues = []
        for residue in self.all_residues:
            copied_residue = deepcopy(residue)
            copied_residues.append(copied_residue)
        copied_conformation = Protein(copied_residues)
        return copied_conformation

    def calculate_energy(self):
        """Returns the overall energy of a protein obtained from
        all hydrophobic contacts."""
        hydrophobic_contacts = set()
        for res in self.all_residues:
            if res.HP_type == "H":
                occupied_neighbors = self.get_occupied_topological_neighbors(res)
                # if occupied neighbors exist
                if occupied_neighbors:
                    # for every neighbor check type and if the interaction
                    # has not been added
                    for neighbor in occupied_neighbors:
                        if neighbor.HP_type == "H" and (neighbor.index, res.index) not in hydrophobic_contacts:
                            hydrophobic_contacts.add((res.index, neighbor.index))
        # energy is equal to the number of hydrophobic contacts
        return -len(hydrophobic_contacts)

    @staticmethod
    def describe_direction(x1, y1, x2, y2):
        """Returns movement direction (U for up, D for down, L for left
        and R for right from one residue to another given two sets of
        coordinates. """
        i1, j1 = x1, y1
        i2, j2 = x2, y2
        if j1 < j2:
            return "U"
        elif j1 > j2:
            return "D"
        elif i1 < i2:
            return "R"
        elif i1 > i2:
            return "L"

    def show(self):
        """Returns a text representation that traces a protein
        with direction symbols."""
        representation = ""
        for i in range(len(self.all_residues) - 1):
            i1, j1 = self.all_residues[i].get_coordinates()
            i2, j2 = self.all_residues[i + 1].get_coordinates()
            direction = self.describe_direction(i1, j1, i2, j2)
            if direction:
                representation += direction

        return representation

    def grid_show(self):
        """Returns a command line graphical visualisation of a protein, note
        that it becomes unreadable for proteins with 10 or more residues"""
        representation = ''
        coordinates = list(map(lambda res: res.get_coordinates(), self.all_residues))
        if len(coordinates) == 0:
            return representation

        min_x = min(coord[0] for coord in coordinates)
        max_x = max(coord[0] for coord in coordinates)
        min_y = min(coord[1] for coord in coordinates)
        max_y = max(coord[1] for coord in coordinates)

        grid_width = max_x - min_x + 1
        grid_height = max_y - min_y + 1

        grid = [[' ' for _ in range(grid_width)] for _ in range(grid_height)]

        count = 0

        for coord in coordinates:
            x, y = coord
            grid_y = y - min_y
            grid_x = x - min_x
            grid[grid_y][grid_x] = str(count)
            count += 1

        for row in grid:
            representation += ''.join(row) + '\n'

        return representation

    def graph_show(self, file_name):
        """Creates the DOT source & Visualization PNG file to
        represent the directed graph"""
        from graphviz import Digraph

        g = Digraph('G', engine="neato", filename=file_name, format='png')
        g.attr(size=str(len(self.all_residues)))

        # creating the nodes by looping over the residues and taking into
        # account their HP_type and their coordinates
        for residue in self.all_residues:
            coordinates = residue.get_coordinates()
            color = 'blue' if residue.HP_type == 'H' else 'red'
            g.node(f"{str(residue.index)}", pos=f"{coordinates[0]},{coordinates[1]}!", fillcolor=color, style='filled')
        # creating the edges directed from the previous residue
        # to the current one
        for i in range(1, len(self.all_residues)):
            g.edge(f"{str(i - 1)}", f"{str(i)}")
        g.render()
