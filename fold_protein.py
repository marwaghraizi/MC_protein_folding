import argparse
import random
import sys

from protein import Protein
from residue import Residue
from manipulation import Manipulation


def translate_to_HP(_sequence):
    polar_residues = ["E", "D", "H", "T", "S", "Y", "N", "Q", "R", "K", "H"]
    hydrophobic_residues = ["C", "W", "G", "A", "P", "I", "L", "M", "F", "V"]
    hp_sequence = ""
    # if the sequence is already in HP format
    if _sequence.strip("HP") == "":
        print(f"The HP model of the protein is: {_sequence}")
        return _sequence
    else:
        for aa in _sequence:
            if aa in hydrophobic_residues:
                hp_sequence += "H"
            elif aa in polar_residues:
                hp_sequence += "P"
            else:
                sys.exit("You entered an invalid amino acid one letter symbol."
                         " Exiting the program.")
        print(f"The HP model of the protein is: {hp_sequence}")
    return hp_sequence


def create_list_of_residue_objects(HP_sequence):
    residues = []
    for idx, res in enumerate(HP_sequence):
        residue_object = Residue(res, idx)
        residue_object.set_coordinates(idx, 0)
        residues.append(residue_object)
    return residues


def randomize_protein(_protein):
    new_protein = _protein.copy_protein()
    for i in range(len(new_protein.all_residues)-1):
        res = new_protein.all_residues[i]
        next_res = new_protein.all_residues[i+1]
        neighbors = new_protein.get_empty_topological_positions(res)
        if neighbors:
            new_coordinates = random.choice(neighbors)
            new_protein.set_coordinates(next_res, *new_coordinates)
        else:
            return randomize_protein(_protein)
    return new_protein


if __name__ == '__main__':
    random.seed(5)
    parser = argparse.ArgumentParser(description="Fold HP Protein")
    parser.add_argument('-f', '--file', type=str,
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
    parser.add_argument('--display-final-frame', action='store_true',
                        help='Display final frame as directions trace.')
    parser.add_argument('--display-graph', choices=["final", "optimal"],
                        default="final",
                        help='create png of final or optimal frame.')
    args = parser.parse_args()

    initial_conformation = args.initial_conformation
    n_iterations = args.n_iterations
    temperature = args.temperature
    search_space = args.search_space

    sequence = ""
    # if the input is a fasta file, read the sequence and transform it to HP if
    # needed
    if args.file:
        with open(args.file, "r") as filein:
            read_sequence = ""
            for line in filein:
                if not line.startswith(">"):
                    read_sequence += line.strip()
            sequence = translate_to_HP(read_sequence)

    # if the sequence is a command line input, transform it to HP if needed
    if args.protein:
        sequence = translate_to_HP(args.protein)

    # if no sequence was provided, throw error and exit
    if sequence == "":
        sys.exit("No input protein was provided. Exiting.")

    # create a list of residue objects from input sequence
    all_residues = create_list_of_residue_objects(sequence)

    # create initial protein which is linear
    initial_linear_protein = Protein(all_residues)

    # if the initial conformation desired is random
    if initial_conformation == 'random':
        protein = randomize_protein(initial_linear_protein)
        protein.graph_show("initial_randomized_protein")
        starting_energy = protein.calculate_energy()
        print(f"Initial conformation of the protein is randomized with a "
              f"starting energy of {starting_energy}.\nFind the graph "
              f"representation in initial_randomized_protein.png")
    else:
        protein = initial_linear_protein
        print("Initial conformation of the protein is linear with a starting"
              "energy of 0")

    # create manipulation object
    manipulation = Manipulation(n_iterations, temperature)
    # add the initial protein as the first frame
    manipulation.add_frame(protein, 'initial input')
    # apply the monte carlo protein folding algorithm in the desired search
    # neighborhood
    manipulation.apply_monte_carlo(search_space)
    # retrieve final frame
    final_frame = manipulation.all_frames[-1]
    # retrieve optimal frame
    optimal_protein = manipulation.get_optimal_frame()

    # shows the final frame as traced directions
    if args.display_final_frame:
        final_frame = manipulation.all_frames[-1]
        print(f"The final frame can be traced as {final_frame.show()}")

    if args.display_graph == "final":
        if args.file:
            file_name = args.file.rsplit(".")[0]
            final_frame.graph_show(f"{file_name}_final_display")
            print(f"Graph representation is saved to {file_name}_"
                  f"final_display.png")
        elif args.protein:
            final_frame.graph_show("protein_display")
            print(f"Graph representation is saved to "
                  f"protein_final_display.png")
    elif args.display_graph == "optimal":
        if args.file:
            file_name = args.file.rsplit(".")[0]
            optimal_protein.graph_show(f"{file_name}_optimal_display")
            print(f"Graph representation is saved to {file_name}_"
                  f"optimal_display.png")
        elif args.protein:
            optimal_protein.graph_show("protein_display")
            print(f"Graph representation is saved to "
                  f"protein_optimal_display.png")

    with open("log.txt", "w") as filout:
        filout.write(f"The HP model of the protein is: {sequence}\n")
        filout.write(f"The starting energy is {protein.calculate_energy()}\n")
        filout.write(f"The final energy is {final_frame.calculate_energy()}\n")
        filout.write(f"The optimal energy is "
                     f"{optimal_protein.calculate_energy()}\n")
        filout.write("-----------------------------")
        for idx, conformation in enumerate(manipulation.all_frames):
            energy = conformation.calculate_energy()
            print(f"-------- Frame {idx} - Caused by {manipulation.moves[idx]} "
                  f"- Energy: {energy} -------", file=filout)
            print(conformation.show(), file=filout)
        filout.close()










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