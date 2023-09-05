import protein
import residue


class Move:
    def is_available(self, conformation, coordinates):
        if coordinates in conformation.coordinates:
            return False
        else:
            return True

    def corner_move(self, conformation, res):
        # if i-1 and i+1 share an available position --> move i to this position
        previous_residue = conformation.get_previous_residue(res)
        next_residue = conformation.get_next_residue(res)
        previous_residue_neighbors = previous_residue.get_topological_neighbors()
        next_residue_neighbors = next_residue.get_topological_neighbors()

        for neighbor in previous_residue_neighbors:
            if neighbor in next_residue_neighbors and neighbor not in conformation.coordinates:
                res.set_coordinates(*neighbor)

        return()

    def end_move(self, conformation, res):

        return Protein

    def crankshaft_move(self, res):
        return Protein

    def choose_random_move(self):
        return Move

