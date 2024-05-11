#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import random
from peptides import Peptide


class GeneticAlgorithm:
    class Solution:
        def __init__(self, sequence, molecular_weight, isoelectric_point, hydrophobicity):
            self.sequence = sequence
            self.molecular_weight = molecular_weight
            self.isoelectric_point = isoelectric_point
            self.hydrophobicity = hydrophobicity

    def __init__(self, target_molecular_weight, target_isoelectric_point, target_hydrophobicity,
                 population_size, offspring_size, num_generations, mutation_probability):
        self.target_molecular_weight = target_molecular_weight
        self.target_isoelectric_point = target_isoelectric_point
        self.target_hydrophobicity = target_hydrophobicity
        self.population_size = population_size
        self.offspring_size = offspring_size
        self.num_generations = num_generations
        self.mutation_probability = mutation_probability

    def calculate_fitness(self, peptide):
        mw_diff = abs(peptide.molecular_weight() - self.target_molecular_weight)
        pI_diff = abs(peptide.isoelectric_point() - self.target_isoelectric_point)
        hydro_diff = abs(peptide.hydrophobicity() - self.target_hydrophobicity)
        return mw_diff + pI_diff + hydro_diff

    def generate_random_peptide(self):
        return Peptide(sequence="".join(random.choices("ACDEFGHIKLMNPQRSTVWY", k=10)))

    def generate_random_population(self):
        population = []

        for _ in range(self.population_size):
            peptide = self.generate_random_peptide()
            fitness = self.calculate_fitness(peptide)
            solution = self.Solution(peptide.sequence, peptide.molecular_weight(), peptide.isoelectric_point(),
                                     peptide.hydrophobicity())
            population.append((solution, fitness))

        return population

    def crossover(self, parent1, parent2):
        crossover_point = random.randint(1, len(parent1.sequence) - 1)
        child_sequence = parent1.sequence[:crossover_point] + parent2.sequence[crossover_point:]
        return Peptide(sequence=child_sequence)

    def mutate(self, peptide):
        mutated_sequence = list(peptide.sequence)
        mutation_point = random.randint(0, len(mutated_sequence) - 1)
        mutated_sequence[mutation_point] = random.choice("ACDEFGHIKLMNPQRSTVWY")
        return Peptide(sequence="".join(mutated_sequence))

    def generate_offspring(self, population):
        offspring = []

        for _ in range(self.offspring_size):
            parent1 = random.choices(population, weights=[1 / fitness for (_, fitness) in population])[0][0]
            parent2 = random.choices(population, weights=[1 / fitness for (_, fitness) in population])[0][0]

            child = self.crossover(parent1, parent2)
            child = self.mutate(child)

            fitness = self.calculate_fitness(child)
            offspring.append((self.Solution(child.sequence, child.molecular_weight(), child.isoelectric_point(),
                                            child.hydrophobicity()), fitness))

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

            best_solution = population[0][0]
            print(f"Generation {generation + 1}: Best Peptide - {best_solution.sequence}")
            print("Molecular Weight:", best_solution.molecular_weight)
            print("Isoelectric Point:", best_solution.isoelectric_point)
            print("Hydrophobicity:", best_solution.hydrophobicity)
            print("Fitness:", population[0][1])
            print()

            if population[0][1] == 0:
                print("Target properties achieved!")
                break

        return population[0][0].sequence, population[0][0].molecular_weight, population[0][0].isoelectric_point, \
        population[0][0].hydrophobicity


if __name__ == "__main__":
    target_molecular_weight = 1500
    target_isoelectric_point = 7.0
    target_hydrophobicity = -0.5
    population_size = 100
    offspring_size = 100
    num_generations = 50
    mutation_probability = 0.1

    genetic_algorithm = GeneticAlgorithm(target_molecular_weight, target_isoelectric_point, target_hydrophobicity,
                                         population_size, offspring_size, num_generations, mutation_probability)
    best_sequence, molecular_weight, isoelectric_point, hydrophobicity = genetic_algorithm.optimize()

    print("\nOptimized Peptide:")
    print("Sequence:", best_sequence)
    print("Molecular Weight:", molecular_weight)
    print("Isoelectric Point:", isoelectric_point)
    print("Hydrophobicity:", hydrophobicity)
