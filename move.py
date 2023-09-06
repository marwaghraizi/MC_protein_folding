from protein import Protein
import random as rd

class Move:
    # apply the move and return true if it was successful
    # check if a position is empty
    def is_available(self, conformation, coordinates):
        if coordinates in conformation.coordinates:
            return False
        else:
            return True

    # moves (not sure what to return)
    def corner_move(self, conformation, res):
        # if i-1 and i+1 share an available position --> move i to this position
        previous_residue = conformation.get_previous_residue(res)
        next_residue = conformation.get_next_residue(res)
        previous_residue_neighbors = previous_residue.get_topological_neighbors()
        next_residue_neighbors = next_residue.get_topological_neighbors()

        for neighbor in previous_residue_neighbors:
            if neighbor in next_residue_neighbors and neighbor not in conformation.coordinates:
                res.set_coordinates(*neighbor)
                return True
        return False

    def end_move(self, conformation, res):
        # ensure that its an extremity
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

        topological_positions = conformation.get_topological_positions(neighbor)
        options = []

        # choose topological positions of the neighbor that are not occupied
        for position in topological_positions:
            if position not in conformation.coordinates:
                options.add(position)

        probability = rd.random()
        if probability < 0.5:
            res.set_coordinates(options[0])
            return True
        else:
            res.set_coordinates(options[1])
            return True
        return False

    def crankshaft_move(self, res):

        return False

    def choose_random_move(self):
        return False

