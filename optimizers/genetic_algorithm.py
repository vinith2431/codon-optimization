"""
Genetic Algorithm based codon optimization.

The chromosome represents synonymous codon choices
for the amino-acid sequence encoded by the input DNA.
"""

import random
import math

from core.sequence import get_codons
from core.verification import CODON_TO_AA
from core.cai import calculate_cai
from data.codon_tables import get_codon_table

from constraints.gc import gc_content
from constraints.local_gc import local_gc_in_range
from constraints.repeats import has_repeats
from constraints.homopolymer import has_homopolymer
from constraints.restriction_sites import contains_restriction_site


class GeneticAlgorithmOptimizer:

    def __init__(
        self,
        sequence,
        organism,
        population_size=50,
        generations=100,
        mutation_rate=0.05,
        crossover_rate=0.8,
        elite_size=2,
        gc_min=None,
        gc_max=None,
        local_gc_window=None,
        local_gc_min=None,
        local_gc_max=None,
        repeat_length=None,
        max_homopolymer=5,
        forbidden_sites=None,
        seed=None
    ):
        """
        Initialize the GA optimizer.
        """

        self.sequence = sequence.upper()
        self.organism = organism

        self.population_size = population_size
        self.generations = generations

        self.mutation_rate = mutation_rate
        self.crossover_rate = crossover_rate

        self.elite_size = elite_size

        self.gc_min = gc_min
        self.gc_max = gc_max

        self.local_gc_window = local_gc_window
        self.local_gc_min = local_gc_min
        self.local_gc_max = local_gc_max

        self.repeat_length = repeat_length

        self.max_homopolymer = max_homopolymer

        self.forbidden_sites = (
            forbidden_sites
            if forbidden_sites is not None
            else []
        )

        if seed is not None:
            random.seed(seed)

        self.original_codons = get_codons(
            self.sequence
        )

        self.codon_table = get_codon_table(
            self.organism
        )

        self.amino_acids = [
            CODON_TO_AA[codon]
            for codon in self.original_codons
        ]

    # --------------------------------------------------
    # Chromosome creation
    # --------------------------------------------------

    def create_random_chromosome(self):
        """
        Create a chromosome by randomly selecting a
        synonymous codon at every amino-acid position.
        """

        chromosome = []

        for amino_acid in self.amino_acids:

            choices = list(
                self.codon_table[amino_acid].keys()
            )

            chromosome.append(
                random.choice(choices)
            )

        return chromosome

    def chromosome_to_sequence(
        self,
        chromosome
    ):
        """
        Convert chromosome to DNA sequence.
        """

        return "".join(chromosome)

    # --------------------------------------------------
    # Population
    # --------------------------------------------------

    def initialize_population(self):
        """
        Create the initial population.
        """

        return [
            self.create_random_chromosome()
            for _ in range(self.population_size)
        ]

    # --------------------------------------------------
    # Objective
    # --------------------------------------------------

    def codon_usage_score(
        self,
        chromosome
    ):
        """
        Calculate average codon usage frequency.
        """

        frequencies = []

        for codon, amino_acid in zip(
            chromosome,
            self.amino_acids
        ):

            frequency = self.codon_table[
                amino_acid
            ][codon]

            frequencies.append(
                frequency
            )

        if not frequencies:
            return 0.0

        return sum(frequencies) / len(frequencies)

    # --------------------------------------------------
    # Constraint penalty
    # --------------------------------------------------

    def calculate_penalty(
        self,
        sequence
    ):
        """
        Calculate penalties for constraint violations.
        """

        penalty = 0.0

        # Global GC
        gc = gc_content(sequence)

        if (
            self.gc_min is not None
            and gc < self.gc_min
        ):
            penalty += (
                self.gc_min - gc
            ) * 100

        if (
            self.gc_max is not None
            and gc > self.gc_max
        ):
            penalty += (
                gc - self.gc_max
            ) * 100

        # Local GC
        if self.local_gc_window is not None:

            valid = local_gc_in_range(
                sequence,
                self.local_gc_window,
                self.local_gc_min,
                self.local_gc_max
            )

            if not valid:
                penalty += 10.0

        # Repeats
        if self.repeat_length is not None:

            if has_repeats(
                sequence,
                self.repeat_length
            ):
                penalty += 10.0

        # Homopolymer
        if has_homopolymer(
            sequence,
            self.max_homopolymer
        ):
            penalty += 10.0

        # Restriction sites
        if contains_restriction_site(
            sequence,
            self.forbidden_sites
        ):
            penalty += 20.0

        return penalty

    # --------------------------------------------------
    # Fitness
    # --------------------------------------------------

    def fitness(
        self,
        chromosome
    ):
        """
        Calculate chromosome fitness.

        Higher fitness is better.
        """

        sequence = self.chromosome_to_sequence(
            chromosome
        )

        usage_score = self.codon_usage_score(
            chromosome
        )

        penalty = self.calculate_penalty(
            sequence
        )

        return usage_score - penalty

    # --------------------------------------------------
    # Selection
    # --------------------------------------------------

    def tournament_selection(
        self,
        population,
        tournament_size=3
    ):
        """
        Select one chromosome using tournament selection.
        """

        participants = random.sample(
            population,
            min(
                tournament_size,
                len(population)
            )
        )

        return max(
            participants,
            key=self.fitness
        ).copy()

    # --------------------------------------------------
    # Crossover
    # --------------------------------------------------

    def crossover(
        self,
        parent1,
        parent2
    ):
        """
        Single-point crossover.
        """

        if (
            random.random()
            > self.crossover_rate
        ):
            return (
                parent1.copy(),
                parent2.copy()
            )

        if len(parent1) < 2:

            return (
                parent1.copy(),
                parent2.copy()
            )

        point = random.randint(
            1,
            len(parent1) - 1
        )

        child1 = (
            parent1[:point]
            + parent2[point:]
        )

        child2 = (
            parent2[:point]
            + parent1[point:]
        )

        return child1, child2

    # --------------------------------------------------
    # Mutation
    # --------------------------------------------------

    def mutate(
        self,
        chromosome
    ):
        """
        Mutate codon choices while preserving amino acids.
        """

        chromosome = chromosome.copy()

        for i, amino_acid in enumerate(
            self.amino_acids
        ):

            if random.random() > self.mutation_rate:
                continue

            choices = list(
                self.codon_table[
                    amino_acid
                ].keys()
            )

            if len(choices) <= 1:
                continue

            current = chromosome[i]

            alternatives = [
                codon
                for codon in choices
                if codon != current
            ]

            chromosome[i] = random.choice(
                alternatives
            )

        return chromosome

    # --------------------------------------------------
    # Optimize
    # --------------------------------------------------

    def optimize(self):
        """
        Run the genetic algorithm.
        """

        population = (
            self.initialize_population()
        )

        best_chromosome = None
        best_fitness = -math.inf

        for generation in range(
            self.generations
        ):

            population.sort(
                key=self.fitness,
                reverse=True
            )

            current_best = population[0]
            current_fitness = self.fitness(
                current_best
            )

            if current_fitness > best_fitness:

                best_fitness = current_fitness
                best_chromosome = (
                    current_best.copy()
                )

            new_population = []

            # Elitism
            elite_count = min(
                self.elite_size,
                len(population)
            )

            for chromosome in population[
                :elite_count
            ]:

                new_population.append(
                    chromosome.copy()
                )

            # Generate remaining population
            while len(new_population) < self.population_size:

                parent1 = self.tournament_selection(
                    population
                )

                parent2 = self.tournament_selection(
                    population
                )

                child1, child2 = self.crossover(
                    parent1,
                    parent2
                )

                child1 = self.mutate(
                    child1
                )

                child2 = self.mutate(
                    child2
                )

                new_population.append(
                    child1
                )

                if (
                    len(new_population)
                    < self.population_size
                ):
                    new_population.append(
                        child2
                    )

            population = new_population

        optimized_sequence = (
            self.chromosome_to_sequence(
                best_chromosome
            )
        )

        return optimized_sequence


def optimize_ga(
    sequence,
    organism,
    **kwargs
):
    """
    Convenience function for GA optimization.
    """

    optimizer = GeneticAlgorithmOptimizer(
        sequence=sequence,
        organism=organism,
        **kwargs
    )

    return optimizer.optimize()