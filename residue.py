import protein


class Residue:
    def __init__(self, type, index):
        self.type = type
        self.index = index
        self.coordI = None
        self.coordJ = None

    def get_coordinates(self):
        return self.coordI, self.coordJ

    def set_coordinates(self, new_coordI, new_coordJ):
        self.coordI = new_coordI
        self.coordJ = new_coordJ

    def get_topological_neighbors(self):
        return (self.coordI - 1, self.coordJ), (self.coordI + 1, self.coordJ), (self.coordI,  self.coordJ - 1), (self.coordI, self.coordJ + 1)


