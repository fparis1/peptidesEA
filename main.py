import random
import numpy as np
from peptides import Peptide
from genetic_algorithm import GeneticAlgorithm


def compute_properties(sequence_custom):
    # Function to compute properties of the peptide sequence
    molecular_weight_custom = sequence_custom.molecular_weight()
    isoelectric_point_custom = sequence_custom.isoelectric_point()
    hydrophobicity_custom = sequence_custom.hydrophobicity()
    charge_custom = sequence_custom.charge()
    aliphatic_index_custom = sequence_custom.aliphatic_index()
    instability_index_custom = sequence_custom.instability_index()
    boman_custom = sequence_custom.boman()
    hydrophobic_moment_custom = sequence_custom.hydrophobic_moment()
    return (molecular_weight_custom, isoelectric_point_custom, hydrophobicity_custom, charge_custom,
            aliphatic_index_custom, instability_index_custom, boman_custom, hydrophobic_moment_custom)


def generate_random_sequence(length):
    # Function to generate a random peptide sequence
    return Peptide(sequence="".join(random.choices("ACDEFGHIKLMNPQRSTVWY", k=length)))


if __name__ == "__main__":
    # Switch to choose how to generate target values
    mode = input("Choose mode: (1) Random peptide sequence, (2) Manually define values, (3) Read from file: ")

    if mode == '1':
        # Generate a random peptide sequence and compute its properties
        random_sequence = generate_random_sequence(np.random.randint(15, 50))  # Example length of 50
        target_values = compute_properties(random_sequence)
        target_molecular_weight, target_isoelectric_point, target_hydrophobicity, target_charge, \
            target_aliphatic_index, target_instability_index, target_boman, target_hydrophobic_moment = target_values

        initial_sequence = random_sequence.sequence

    elif mode == '2':
        # Manually define target values
        target_molecular_weight = 3000
        target_isoelectric_point = 9.5
        target_hydrophobicity = 0.4
        target_charge = 5
        target_aliphatic_index = 90
        target_instability_index = 40
        target_boman = 1.0
        target_hydrophobic_moment = 0.6

        initial_sequence = None  # No initial sequence for manually defined values

    elif mode == '3':
        # Read peptide sequence from file and compute its properties
        filename = 'sequence_input.txt'
        try:
            with open(filename, 'r') as file:
                sequence = file.read().strip()
            target_values = compute_properties(Peptide(sequence=sequence))
            (target_molecular_weight, target_isoelectric_point, target_hydrophobicity, target_charge,
             target_aliphatic_index, target_instability_index, target_boman,
             target_hydrophobic_moment) = target_values

            initial_sequence = sequence
        except FileNotFoundError:
            print(f"File '{filename}' not found.")
            exit(1)

    else:
        print("Invalid mode selected.")
        exit(1)

    # Parameters for the genetic algorithm
    population_size = 80
    offspring_size = 50
    num_generations = 300
    mutation_probability = 0.3
    tournament_size = 3

    # Initialize the genetic algorithm with the target values
    genetic_algorithm = GeneticAlgorithm(target_molecular_weight, target_isoelectric_point, target_hydrophobicity,
                                         target_charge, target_aliphatic_index, target_instability_index,
                                         target_boman, target_hydrophobic_moment,
                                         population_size, offspring_size, num_generations, mutation_probability,
                                         tournament_size)

    # Run the optimization
    best_sequence, molecular_weight, isoelectric_point, hydrophobicity, charge, aliphatic_index, instability_index, \
        boman, hydrophobic_moment = genetic_algorithm.optimize()

    # Print initial sequence and its properties
    print("Mode Selected:", mode)
    if mode == '1' or mode == '3':
        print("\nInitial Peptide:")
        print("Sequence:", initial_sequence)
    elif mode == '2':
        print("\nManually Defined Target Values:")

    print("Target Molecular Weight:", target_molecular_weight)
    print("Target Isoelectric Point:", target_isoelectric_point)
    print("Target Hydrophobicity:", target_hydrophobicity)
    print("Target Charge:", target_charge)
    print("Target Aliphatic Index:", target_aliphatic_index)
    print("Target Instability Index:", target_instability_index)
    print("Target Boman:", target_boman)
    print("Target Hydrophobic Moment:", target_hydrophobic_moment)

    # Print the optimized peptide properties
    print("\nOptimized Peptide:")
    print("Sequence:", best_sequence)
    print("Molecular Weight:", molecular_weight)
    print("Isoelectric Point:", isoelectric_point)
    print("Hydrophobicity:", hydrophobicity)
    print("Charge:", charge)
    print("Aliphatic Index:", aliphatic_index)
    print("Instability Index:", instability_index)
    print("Boman:", boman)
    print("Hydrophobic Moment:", hydrophobic_moment)

    genetic_algorithm.plot_scores()

    # Generate radar chart comparing initial targets and final values
    final_values = [molecular_weight, isoelectric_point, hydrophobicity, charge, aliphatic_index,
                    instability_index, boman, hydrophobic_moment]
    target_values = [target_molecular_weight, target_isoelectric_point, target_hydrophobicity, target_charge,
                     target_aliphatic_index, target_instability_index, target_boman, target_hydrophobic_moment]

    genetic_algorithm.create_radar_chart(target_values, final_values)
