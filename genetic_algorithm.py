import random

import matplotlib.pyplot as plt
import numpy as np
from peptides import Peptide


class GeneticAlgorithm:
    class Solution:
        def __init__(self, sequence, molecular_weight, isoelectric_point, hydrophobicity,
                     charge, aliphatic_index, instability_index, boman, hydrophobic_moment):
            self.sequence = sequence
            self.molecular_weight = molecular_weight
            self.isoelectric_point = isoelectric_point
            self.hydrophobicity = hydrophobicity
            self.charge = charge
            self.aliphatic_index = aliphatic_index
            self.instability_index = instability_index
            self.boman = boman
            self.hydrophobic_moment = hydrophobic_moment

    def __init__(self, target_molecular_weight, target_isoelectric_point, target_hydrophobicity,
                 target_charge, target_aliphatic_index, target_instability_index,
                 target_boman, target_hydrophobic_moment,
                 population_size, offspring_size, num_generations, mutation_probability, tournament_size):
        self.target_molecular_weight = target_molecular_weight
        self.target_isoelectric_point = target_isoelectric_point
        self.target_hydrophobicity = target_hydrophobicity
        self.target_charge = target_charge
        self.target_aliphatic_index = target_aliphatic_index
        self.target_instability_index = target_instability_index
        self.target_boman = target_boman
        self.target_hydrophobic_moment = target_hydrophobic_moment
        self.population_size = population_size
        self.offspring_size = offspring_size
        self.num_generations = num_generations
        self.mutation_probability = mutation_probability
        self.tournament_size = tournament_size
        self.best_scores = []
        self.worst_scores = []

    def calculate_fitness(self, peptide):
        mw_diff = abs(peptide.molecular_weight() - self.target_molecular_weight)
        p_i_diff = abs(peptide.isoelectric_point() - self.target_isoelectric_point)
        hydrophobicity = peptide.hydrophobicity(scale="Eisenberg")
        hydro_diff = abs(hydrophobicity - self.target_hydrophobicity)
        charge_diff = abs(peptide.charge(pH=7.0, pKscale="EMBOSS") - self.target_charge)
        aliphatic_index_diff = abs(self.target_aliphatic_index - peptide.aliphatic_index())
        instability_index_diff = abs(peptide.instability_index() - self.target_instability_index)
        boman_diff = abs(peptide.boman() - self.target_boman)
        hydrophobic_moment_diff = abs(peptide.hydrophobic_moment() - self.target_hydrophobic_moment)

        # Weights for each property
        weight_mw = 10
        weight_pi = 100
        weight_hydro = 1000 if hydro_diff < 0.5 else 0
        weight_charge = 50
        weight_aliphatic = 50
        weight_instability = 50
        weight_boman = 200
        weight_hydrophobic_moment = 200

        # Properties to use MAE
        mae_fitness = (weight_mw * mw_diff +
                       weight_hydro * hydro_diff +
                       weight_aliphatic * aliphatic_index_diff +
                       weight_instability * instability_index_diff)

        # Properties to use MSE
        mse_fitness = (weight_pi * (p_i_diff ** 2) +
                       weight_charge * (charge_diff ** 2) +
                       weight_boman * (boman_diff ** 2) +
                       weight_hydrophobic_moment * (hydrophobic_moment_diff ** 2))

        # Combining MAE and MSE components
        fitness = mae_fitness + mse_fitness

        return fitness

    def generate_random_peptide(self):
        return Peptide(sequence="".join(random.choices("ACDEFGHIKLMNPQRSTVWY", k=np.random.randint(15, 50))))

    def generate_random_population(self):
        population = []

        for _ in range(self.population_size):
            peptide = self.generate_random_peptide()
            fitness = self.calculate_fitness(peptide)
            solution = self.Solution(peptide.sequence, peptide.molecular_weight(), peptide.isoelectric_point(),
                                     peptide.hydrophobicity(), peptide.charge(pH=7.0, pKscale="EMBOSS"),
                                     peptide.aliphatic_index(), peptide.instability_index(),
                                     peptide.boman(), peptide.hydrophobic_moment())
            population.append((solution, fitness))

        return population

    def crossover(self, parent1, parent2):
        min_length = min(len(parent1.sequence), len(parent2.sequence))
        crossover_point = random.randint(1, min_length - 1)
        child_sequence = parent1.sequence[:crossover_point] + parent2.sequence[crossover_point:]
        return Peptide(sequence=child_sequence)

    def mutate(self, peptide):
        if np.random.rand() <= self.mutation_probability:
            mutated_sequence = list(peptide.sequence)
            mutation_type = np.random.choice(['substitution', 'addition', 'deletion', 'swap'])

            if mutation_type == 'substitution':
                mutation_point = random.randint(0, len(mutated_sequence) - 1)
                original_amino_acid = mutated_sequence[mutation_point]

                new_amino_acid = original_amino_acid
                while new_amino_acid == original_amino_acid:
                    new_amino_acid = random.choice("ACDEFGHIKLMNPQRSTVWY")

                mutated_sequence[mutation_point] = new_amino_acid

            elif mutation_type == 'addition':
                addition_point = random.randint(0, len(mutated_sequence))
                new_amino_acid = random.choice("ACDEFGHIKLMNPQRSTVWY")
                mutated_sequence.insert(addition_point, new_amino_acid)

            elif mutation_type == 'deletion':
                deletion_point = random.randint(0, len(mutated_sequence) - 1)
                del mutated_sequence[deletion_point]

            elif mutation_type == 'swap':
                swap_points = random.sample(range(len(mutated_sequence)), 2)
                while mutated_sequence[swap_points[0]] == mutated_sequence[swap_points[1]]:
                    swap_points = random.sample(range(len(mutated_sequence)), 2)
                mutated_sequence[swap_points[0]], mutated_sequence[swap_points[1]] = mutated_sequence[
                    swap_points[1]], mutated_sequence[swap_points[0]]

            return Peptide(sequence="".join(mutated_sequence))
        else:
            return peptide

    def tournament_selection(self, population):
        # Randomly select tournament_size individuals from the population
        tournament = random.sample(population, self.tournament_size)

        # Select the best individual from the tournament
        best_individual = min(tournament, key=lambda solution: solution[1])

        return best_individual[0]

    def generate_offspring(self, population):
        offspring = []

        for _ in range(self.offspring_size):
            parent1 = self.tournament_selection(population)
            parent2 = self.tournament_selection(population)

            child = self.crossover(parent1, parent2)
            child = self.mutate(child)

            fitness = self.calculate_fitness(child)
            offspring.append((self.Solution(child.sequence, child.molecular_weight(), child.isoelectric_point(),
                                            child.hydrophobicity(), child.charge(pH=7.0, pKscale="EMBOSS"),
                                            child.aliphatic_index(), child.instability_index(),
                                            child.boman(), child.hydrophobic_moment()), fitness))

        return offspring

    def next_generation(self, population, offspring):
        combined_population = population + offspring
        combined_population.sort(key=lambda x: x[1])
        return combined_population[:self.population_size]

    def optimize(self):
        population = self.generate_random_population()

        for generation in range(self.num_generations):
            offspring = self.generate_offspring(population)
            population = self.next_generation(population, offspring)

            best_fitness = population[0][1]
            worst_fitness = population[-1][1]
            self.best_scores.append(best_fitness)
            self.worst_scores.append(worst_fitness)
            print("Generation:", generation + 1, "/", self.num_generations)

            if population[0][1] == 0:
                print("Target properties achieved!")
                break

        return (population[0][0].sequence, population[0][0].molecular_weight, population[0][0].isoelectric_point,
                population[0][0].hydrophobicity, population[0][0].charge, population[0][0].aliphatic_index,
                population[0][0].instability_index, population[0][0].boman,
                population[0][0].hydrophobic_moment, self.best_scores, self.worst_scores)