import Residue
import Move

class Protein:
    def __init__(self, residues):
        self.ALL_RESIDUES = residues
        self.energy = int()
        self.coordinates = {}
        # will it be updated?
        for residue in residues:
            self.coordinates[(residue.coordI, residue.coordJ)] = residue

    def get_occupied_topological_neighbors(self, res):
        #returns topological neighbors
        connected_residues = self.get_connected_residues(res)
        topological_positions = (res.coordI - 1, res.coordJ), (res.coordI+1, res.coordJ), \
                                (res.coordI, res.coordJ - 1), (res.coordI, res.coordJ + 1)
        # topological neighbors are in the topological positions and are not connected to the residue of interest
        occupied_topological_neighbors = (res for x in self.ALL_RESIDUES if x not in connected_residues
                                          and (x.coordI, x.coordJ) in topological_positions)
        return occupied_topological_neighbors

    def get_topological_neighbors(self, res):
        topological_positions = [(res.coordI - 1, res.coordJ), (res.coordI+1, res.coordJ),
                                 (res.coordI, res.coordJ - 1), (res.coordI, res.coordJ + 1)]
        return topological_positions

    def get_connected_residues(self, res):
        connected_residues = [] # should it be a ?
        if res.index == 0:
            connected_residues.append(self.ALL_RESIDUES[res.index+1])
        elif res.index == len(self.ALL_RESIDUES)-1:
            connected_residues.append(self.ALL_RESIDUES[res.index-1])
        else:
            connected_residues.append([self.ALL_RESIDUES[res.index - 1], self.ALL_RESIDUES[res.index + 1]])
        return connected_residues

    def get_next_residue(self, res):
        return self.ALL_RESIDUES[res.index+1]

    def get_previous_residue(self, res):
        return self.ALL_RESIDUES[res.index-1]

    def calculate_energy(self):
        hydrophobic_contacts = set()
        for res in self.ALL_RESIDUES:
            if res.type == "H":
                neighbors = self.get_topological_neighbors(res)
                for neighbor in neighbors:
                    if neighbor.type == "H" and (neighbor.index, res.index) not in hydrophobic_contacts:
                        hydrophobic_contacts.add((res.index, neighbor.index))
        self.energy = -len(hydrophobic_contacts)
        return()
