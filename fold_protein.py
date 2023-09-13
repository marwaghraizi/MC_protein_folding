import random
import argparse
from protein import Protein
from residue import Residue
from manipulation import Manipulation

if __name__ == '__main__':
    #random.seed(5)
    polar_residues = ["E", "D", "H", "T", "S", "Y", "N", "Q", "R", "K", "H"]
    hydrophobic_residues = ["C", "W", "G", "A", "P", "I", "L", "M", "F", "V"]

    parser = argparse.ArgumentParser(description="Fold HP Protein")
    parser.add_argument('-f', '--file', type=str, default='test.fasta',
                        help='Protein File Path')
    parser.add_argument('-p', '--protein',
                        help="input protein sequence in classic or HP format")
    parser.add_argument('-i', '--initial-conformation',
                        choices=["linear", "random"], default="linear",
                        help='initial conformation of the protein: linear '
                             'or randomized placements')
    parser.add_argument('-n', '--n-iterations', type=int, default=1000,
                        help="number of search iterations")
    parser.add_argument('-t', '--temperature', type=float, default=100.0,
                        help='search temperature: higher temperatures '
                             'increases the probability of accepting '
                             'energetically unfavorable moves')
    parser.add_argument('-s', '--search-space',
                        choices=["VSHD", "VSHD-pull", 'pull'],
                        default='VSHD-pull')
    parser.add_argument('--probability-pull', type=float, default=0.5)
    parser.add_argument('--display-grid-all', action='store_true',
                        help='print out all of the frames')
    parser.add_argument('--display-grid', action='store_true',
                        help='print out final frame')
    parser.add_argument('--display-graph', choices=["final", "linear"],
                        default="final",
                        help='create png of final or optimal frame.')
    #parser.set_defaults(display_grid=False)
    args = parser.parse_args()
    protein_file = args.file
    initial_confirmation = args.initial_conformation
    n_iterations = args.n_iterations
    temperature = args.temperature
    search_space = args.search_space

    with open(protein_file, "r") as filein:
        sequence = ""
        for line in filein:
            if not line.startswith(">"):
                sequence += line.strip()

    HP_sequence = ""
    ALL_RESIDUES = []

    for idx, residue in enumerate(sequence):
        if residue in polar_residues:
            HP_sequence += "P"
            residue_object = Residue("P", idx)
            residue_object.set_coordinates(idx, 0)
        else:
            HP_sequence += "H"
            residue_object = Residue("H", idx)
            residue_object.set_coordinates(idx, 0)
        ALL_RESIDUES.append(residue_object)

    if initial_confirmation == 'random':
        # occupied_positions = [residue.get_coordinates() for residue in ALL_RESIDUES]
        # print(occupied_positions)
        #
        # for i in range(len(ALL_RESIDUES) - 1):
        #     neighbors = get_four_neighbors(*(ALL_RESIDUES[i].get_coordinates()))
        #     random_neighbor = rd.choice(neighbors)
        #     isFree = False
        #
        #     if is_free(*(random_neighbor)):
        #         isFree = True
        #
        #     while (not isFree):
        #         random_neighbor = rd.choice(neighbors)
        #         if is_free(*(random_neighbor)):
        #             isFree = True
        #
        #     ALL_RESIDUES[i + 1].set_coordinates(*(random_neighbor))
        #     occupied_positions = [residue.get_coordinates() for residue in ALL_RESIDUES]
        pass
    initial_protein = Protein(ALL_RESIDUES)
    manipulation = Manipulation(n_iterations, temperature)
    manipulation.add_frame(initial_protein, 'initial')
    manipulation.apply_monte_carlo(search_space)
    final_frame = manipulation.all_frames[-1]
    if args.display_grid_all:
        manipulation.show_all_frames()
    if args.display_grid:
        print(final_frame.show())
        print(final_frame.grid_show())
    if args.display_graph == "final":
        final_frame.graph_show(f"{protein_file}_display")
    print(f"----- Final Energy {final_frame.calculate_energy()} -----")
    for residue in final_frame.all_residues:
        print(f"residue {residue.index} with coordinates {residue.get_coordinates()}")

    # if the sequence was given in HP format keep it as is if not (if not sequence.strip("HP") == "" --> translate it
    # initial conformation: linear or random (random can start with res 0 on (0,0) and place the rest randomly relative to each other
    # test cases from paper:
    ## HPHPPHHPHPPHPHHPPHPH
    ## HHPPHPPHPPHPPHPPHPPHPPHH
    ## PPHPPHHPPPPHHPPPPHHPPPPHH
    ## PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP
    ## PPHPPHHPPHHPPPPPHHHHHHHHHHPPPPPPPPHHPPHHHPPHHHHH

    # throw error if amino acid sequence contains invalid amino acid symbols (example: X)
    # throw error if n-iterations less than 1
    # throw error if fasta file in wrong format
    # throw error if sequence length < 2
    # input log file name to be created

    # initial_protein = Protein(ALL_RESIDUES)

    # testing corner
    # initial_protein.all_residues[0].coordI = 1
    # initial_protein.all_residues[0].coordJ = 1

    # testing crankshaft
    # initial_protein.all_residues[1].coordI = 0
    # initial_protein.all_residues[1].coordJ = -1
    # initial_protein.all_residues[2].coordI = 1
    # initial_protein.all_residues[2].coordJ = -1
    # initial_protein.all_residues[3].coordI = 1
    # initial_protein.all_residues[3].coordJ = 0

    # Figure 3 (a) testing that pull can do the corner move
    #     ALL_RESIDUES[0].coordJ = 3
    #     ALL_RESIDUES[1].coordI = 0
    #     ALL_RESIDUES[1].coordJ = 2
    #     ALL_RESIDUES[2].coordI = 0
    #     ALL_RESIDUES[2].coordJ = 1
    #     ALL_RESIDUES[3].coordI = 1
    #     ALL_RESIDUES[3].coordJ = 1

    # Figure 3 (b)
    #     ALL_RESIDUES[0].coordJ = 1
    #     ALL_RESIDUES[1].coordI = 0
    #     ALL_RESIDUES[2].coordI = 0
    #     ALL_RESIDUES[3].coordI = 2
    #     ALL_RESIDUES[2].coordI = 1

    # Figure 3 (c)
    # FOR THIS ONE, CHANGE test.fasta to 'AEAEAEAEA' (9 residue)

    # ALL_RESIDUES[0].coordI = 3
    # ALL_RESIDUES[1].coordI = 2
    # ALL_RESIDUES[2].coordI = 1
    # ALL_RESIDUES[3].coordI = 1
    # ALL_RESIDUES[3].coordJ = 1
    # ALL_RESIDUES[4].coordI = 0
    # ALL_RESIDUES[4].coordJ = 1
    # ALL_RESIDUES[5].coordI = 0
    # ALL_RESIDUES[5].coordJ = 2
    # ALL_RESIDUES[6].coordI = 1
    # ALL_RESIDUES[6].coordJ = 2
    # ALL_RESIDUES[7].coordI = 2
    # ALL_RESIDUES[7].coordJ = 2
    # ALL_RESIDUES[8].coordI = 3
    # ALL_RESIDUES[8].coordJ = 2

    # print("INITIAL PROTEIN COORDINATES:")
    # print(initial_protein.coordinates)
    # print()
    # manipulation = Manipulation()
    # manipulation.add_frame(initial_protein, 'initial')
    # manipulation.apply_monte_carlo()
    # # figure b
    # # manipulation.pull_move(manipulation.all_frames[0].copy_protein(), initial_protein.all_residues[2])
    # # figure 3 c
    # # manipulation.pull_move(manipulation.all_frames[0].copy_protein(), initial_protein.all_residues[7])
    #
    # manipulation.show_all_frames()