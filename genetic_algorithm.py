import random
from peptides import Peptide
import numpy as np


class GeneticAlgorithm:
    class Solution:
        def __init__(self, sequence, molecular_weight, isoelectric_point, hydrophobicity,
                     charge, aliphatic_index, instability_index, boman, hydrophobic_moment, mass_to_charge_ratio):
            self.sequence = sequence
            self.molecular_weight = molecular_weight
            self.isoelectric_point = isoelectric_point
            self.hydrophobicity = hydrophobicity
            self.charge = charge
            self.aliphatic_index = aliphatic_index
            self.instability_index = instability_index
            self.boman = boman
            self.hydrophobic_moment = hydrophobic_moment
            self.mass_to_charge_ratio = mass_to_charge_ratio

    def __init__(self, target_molecular_weight, target_isoelectric_point, target_hydrophobicity,
                 target_charge, target_aliphatic_index, target_instability_index,
                 target_boman, target_hydrophobic_moment, target_mass_to_charge_ratio,
                 population_size, offspring_size, num_generations, mutation_probability):
        self.target_molecular_weight = target_molecular_weight
        self.target_isoelectric_point = target_isoelectric_point
        self.target_hydrophobicity = target_hydrophobicity
        self.target_charge = target_charge
        self.target_aliphatic_index = target_aliphatic_index
        self.target_instability_index = target_instability_index
        self.target_boman = target_boman
        self.target_hydrophobic_moment = target_hydrophobic_moment
        self.target_mass_to_charge_ratio = target_mass_to_charge_ratio
        self.population_size = population_size
        self.offspring_size = offspring_size
        self.num_generations = num_generations
        self.mutation_probability = mutation_probability

    def calculate_fitness(self, peptide):
        mw_diff = abs(peptide.molecular_weight() - self.target_molecular_weight)
        p_i_diff = abs(peptide.isoelectric_point() - self.target_isoelectric_point)
        hydrophobicity = peptide.hydrophobicity(scale="Eisenberg")
        hydro_diff = abs(hydrophobicity - self.target_hydrophobicity)

        # Adjusting fitness based on hydrophobicity difference
        if hydro_diff <= 0.5:
            hydro_fitness = 0
        else:
            hydro_fitness = 1000 * hydro_diff

        charge_diff = abs(peptide.charge(pH=7.0, pKscale="EMBOSS") - self.target_charge)
        aliphatic_index_diff = abs(peptide.aliphatic_index() - self.target_aliphatic_index)
        instability_index_diff = abs(peptide.instability_index() - self.target_instability_index)
        boman_diff = abs(peptide.boman() - self.target_boman)
        hydrophobic_moment_diff = abs(peptide.hydrophobic_moment() - self.target_hydrophobic_moment)
        mass_to_charge_ratio_diff = abs(peptide.mz() - self.target_mass_to_charge_ratio)

        # Assigning more weight to hydrophobicity, hydrophobic moment, charge, and Boman
        fitness = (10 * mw_diff + 100 * p_i_diff + hydro_fitness + 50 * charge_diff +
                   50 * aliphatic_index_diff + 50 * instability_index_diff +
                   200 * boman_diff + 200 * hydrophobic_moment_diff +
                   10 * mass_to_charge_ratio_diff)

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
                                     peptide.boman(), peptide.hydrophobic_moment(), peptide.mz())
            population.append((solution, fitness))

        return population

    def crossover(self, parent1, parent2):
        min_length = min(len(parent1.sequence), len(parent2.sequence))
        crossover_point = random.randint(1, min_length - 1)
        child_sequence = parent1.sequence[:crossover_point] + parent2.sequence[crossover_point:]
        return Peptide(sequence=child_sequence)

    def mutate(self, peptide):
        mutated_sequence = list(peptide.sequence)
        mutation_point = random.randint(0, len(mutated_sequence) - 1)
        mutated_sequence[mutation_point] = random.choice("ACDEFGHIKLMNPQRSTVWY")
        return Peptide(sequence="".join(mutated_sequence))

    def selection(self, population):
        fitness_scores = [1 / (solution[1] + 1) for solution in population]

        # Create a roulette wheel.
        roulette_wheel = np.cumsum(
            fitness_scores / np.sum(fitness_scores)
        )

        random_number = np.random.rand()
        for index, score in enumerate(roulette_wheel):
            if random_number <= score:
                return population[index][0]

    def generate_offspring(self, population):
        offspring = []

        for _ in range(self.offspring_size):
            parent1 = self.selection(population)
            parent2 = self.selection(population)

            child = self.crossover(parent1, parent2)
            child = self.mutate(child)

            fitness = self.calculate_fitness(child)
            offspring.append((self.Solution(child.sequence, child.molecular_weight(), child.isoelectric_point(),
                                            child.hydrophobicity(), abs(child.charge(pH=7.0, pKscale="EMBOSS")),
                                            abs(child.aliphatic_index()), abs(child.instability_index()),
                                            abs(child.boman()), abs(child.hydrophobic_moment()), abs(child.mz())),
                              fitness))

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
            print("Generation:", generation + 1, "/", self.num_generations)

            if population[0][1] == 0:
                print("Target properties achieved!")
                break

        return (population[0][0].sequence, population[0][0].molecular_weight, population[0][0].isoelectric_point,
                population[0][0].hydrophobicity, population[0][0].charge, population[0][0].aliphatic_index,
                population[0][0].instability_index, population[0][0].boman,
                population[0][0].hydrophobic_moment, population[0][0].mass_to_charge_ratio)
