from genetic_algorithm import GeneticAlgorithm

if __name__ == "__main__":
    target_molecular_weight = 1500
    target_isoelectric_point = 10.0
    target_hydrophobicity = 0.5
    target_charge = 5
    target_aliphatic_index = 50
    target_instability_index = 40
    target_boman = 2
    target_hydrophobic_moment = 0.5
    target_mass_to_charge_ratio = 800

    population_size = 100
    offspring_size = 100
    num_generations = 100
    mutation_probability = 0.1

    genetic_algorithm = GeneticAlgorithm(target_molecular_weight, target_isoelectric_point, target_hydrophobicity,
                                         target_charge, target_aliphatic_index, target_instability_index,
                                         target_boman, target_hydrophobic_moment, target_mass_to_charge_ratio,
                                         population_size, offspring_size, num_generations, mutation_probability)
    best_sequence, molecular_weight, isoelectric_point, hydrophobicity, charge, aliphatic_index, instability_index, \
        boman, hydrophobic_moment, mass_to_charge_ratio = genetic_algorithm.optimize()

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
    print("Mass-to-Charge Ratio:", mass_to_charge_ratio)

