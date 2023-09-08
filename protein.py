from copy import deepcopy


class Protein:

    def __init__(self, residues):
        self.all_residues = residues
        self.coordinates = {}
        # will it be updated?
        for residue in residues:
            self.coordinates[(residue.coordI, residue.coordJ)] = residue
        # is this legal?
        self.energy = self.calculate_energy()

    @staticmethod
    def describe_direction(x1, y1, x2, y2):
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

        return None

    def show(self):
        representation = ""
        for i in range(len(self.all_residues)-1):
            i1, j1 = self.all_residues[i].get_coordinates()
            i2, j2 = self.all_residues[i+1].get_coordinates()
            direction = self.describe_direction(i1, j1, i2, j2)
            if direction:
                representation += direction

        return representation

    def is_free(self, x, y):
        if (x, y) in self.coordinates:
            return False
        else:
            return True

    def copy_protein(self):
        copied_residues = []
        for residue in self.all_residues:
            copied_residue = deepcopy(residue)
            copied_residues.append(copied_residue)
        copied_conformation = Protein(copied_residues)
        return copied_conformation

    def get_connected_residues(self, res):
        connected_residues = []
        if res.index == 0:
            connected_residues.append(self.all_residues[res.index+1])
        elif res.index == len(self.all_residues)-1:
            connected_residues.append(self.all_residues[res.index-1])
        else:
            connected_residues.append(self.all_residues[res.index - 1])
            connected_residues.append(self.all_residues[res.index + 1])

        return connected_residues

    def get_next_residue(self, res):
        return self.all_residues[res.index+1]

    def get_previous_residue(self, res):
        return self.all_residues[res.index-1]

    # returns occupied topological positions that are not connected to the residue
    # to be cleaned?
    def get_occupied_topological_neighbors(self, res):
        # get connected residues (obj)
        connected_residues = self.get_connected_residues(res)
        connected_residues_positions = {}

        for connected_residue in connected_residues:
            connected_residues_positions[(connected_residue.coordI, connected_residue.coordJ)] = connected_residue

        topological_positions = res.get_topological_neighbors_positions()
        # topological neighbors are in the topological positions and are not connected to the residue of interest
        occupied_topological_neighbors = []

        for position in topological_positions:
            # if the position is occupied and is not a connected amino acid
            if (position in self.coordinates) and (position not in connected_residues_positions):
                occupied_topological_neighbors.append(self.coordinates[position])

        return occupied_topological_neighbors

    # returns topological positions that do not contain any residues
    def get_empty_topological_positions(self, res):
        topological_positions = res.get_topological_neighbors_positions()
        empty_topological_positions = []
        for position in topological_positions:
            if position not in self.coordinates:
                empty_topological_positions.append(position)
        return empty_topological_positions

    def calculate_energy(self):
        hydrophobic_contacts = set()
        for res in self.all_residues:
            if res.type == "H":
                occupied_neighbors = self.get_occupied_topological_neighbors(res)
                # if occupied neighbors exist
                if occupied_neighbors:
                    #neighbors = self.get_occupied_topological_neighbors(res)
                    for neighbor in occupied_neighbors:
                        if neighbor.type == "H" and (neighbor.index, res.index) not in hydrophobic_contacts:
                            hydrophobic_contacts.add((res.index, neighbor.index))
        # if i do this then no need to attribute the method in the constructor so which is it?
        self.energy = -len(hydrophobic_contacts)
        return -len(hydrophobic_contacts)
