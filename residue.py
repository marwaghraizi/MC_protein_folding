class Residue:
    """
        A class to represent an amino acid as a residue within a protein.

        ...

        Attributes
        ----------
        index : int
            index of the residue within the protein conformation.
        HP_type : char
            polarity of the residue: H for hydrophobic and P for polar.
        coordI : int
            i coordinate of the residue's position (equivalent to the x coordinate)
        coordJ: int
            j coordinate of the residue's position (equivalent to the y coordinate)

        Methods
        -------
        get_coordinates():
            Retrieves the residue's coordinates as tuple.
        set_coordinates(new_i, new_j):
            Updates the residue's coordinates.
        get_topological_neighbors_positions():
            Retrieves the coordinates of the 4 adjacent topological positions around the residue.
        """

    def __init__(self, hp_type, index):
        """
        Constructs all the necessary attributes for the person object.

        Parameters
        ----------
            hp_type : char
                first name of the person
            index : int
                family name of the person
        """
        self.HP_type = hp_type
        self.index = index
        self.coordI = None
        self.coordJ = None

    def get_coordinates(self):
        """Returns the coordinates of a residue as a tuple."""
        return self.coordI, self.coordJ

    def set_coordinates(self, new_i, new_j):
        """Updates the coordinates of a residue given new coordinates"""
        self.coordI = new_i
        self.coordJ = new_j

    def get_topological_neighbors_positions(self):
        """Returns the four adjacent positions of a residue as a tuple of tuples"""
        return (self.coordI - 1, self.coordJ), (self.coordI + 1, self.coordJ), (self.coordI,  self.coordJ - 1), \
               (self.coordI, self.coordJ + 1)
