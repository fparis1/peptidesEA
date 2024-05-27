from genetic_algorithm import GeneticAlgorithm

if __name__ == "__main__":
    # Define the target values
    target_molecular_weight = 3000
    target_isoelectric_point = 9.5
    target_hydrophobicity = 0.4
    target_charge = 5
    target_aliphatic_index = 90
    target_instability_index = 40
    target_boman = 1.0
    target_hydrophobic_moment = 0.6

    # Parameters for the genetic algorithm
    population_size = 50
    offspring_size = 30
    num_generations = 500
    mutation_probability = 0.1

    # Initialize the genetic algorithm with the target values
    genetic_algorithm = GeneticAlgorithm(target_molecular_weight, target_isoelectric_point, target_hydrophobicity,
                                         target_charge, target_aliphatic_index, target_instability_index,
                                         target_boman, target_hydrophobic_moment,
                                         population_size, offspring_size, num_generations, mutation_probability)

    # Run the optimization
    best_sequence, molecular_weight, isoelectric_point, hydrophobicity, charge, aliphatic_index, instability_index, \
        boman, hydrophobic_moment = genetic_algorithm.optimize()

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
