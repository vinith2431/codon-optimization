"""
Genetic Algorithm for synonymous codon optimization.

This version is compatible with the current project structure.

Objectives:
    - Maximize CAI
    - Respect global GC constraints
    - Optionally minimize codon changes

The GA returns an optimized DNA string because the existing
test suite expects optimize_ga() to return a string.
"""

from __future__ import annotations

import random
import time
from typing import List, Optional, Tuple

from core.cai import calculate_cai
from core.sequence import get_codons, gc_content
from core.verification import (
    CODON_TO_AA,
    translate_dna,
    verify_translation,
)
from data.codon_tables import get_codon_table


# ============================================================
# Helper functions
# ============================================================

def split_codons(sequence: str) -> List[str]:
    """Split DNA sequence into codons."""
    return get_codons(sequence)


def count_changes(
    original_codons: List[str],
    candidate_codons: List[str],
) -> int:
    """Count changed codon positions."""
    return sum(
        original != candidate
        for original, candidate
        in zip(original_codons, candidate_codons)
    )


def edit_rate(
    original_codons: List[str],
    candidate_codons: List[str],
) -> float:
    """Calculate fraction of changed codons."""
    if not original_codons:
        return 0.0

    return (
        count_changes(original_codons, candidate_codons)
        / len(original_codons)
    )


def gc_penalty(
    gc: float,
    gc_min: Optional[float],
    gc_max: Optional[float],
) -> float:
    """
    Calculate the amount by which GC violates the requested range.
    """
    penalty = 0.0

    if gc_min is not None and gc < gc_min:
        penalty += gc_min - gc

    if gc_max is not None and gc > gc_max:
        penalty += gc - gc_max

    return penalty


def is_gc_feasible(
    sequence: str,
    gc_min: Optional[float],
    gc_max: Optional[float],
) -> bool:
    """Check whether sequence satisfies global GC constraint."""
    gc = gc_content(sequence)

    if gc_min is not None and gc < gc_min:
        return False

    if gc_max is not None and gc > gc_max:
        return False

    return True


# ============================================================
# Synonymous codons
# ============================================================

def build_synonymous_codons(
    original_codons: List[str],
    codon_table: dict,
) -> List[List[str]]:
    """
    Build possible synonymous codons for each position.

    Stop codons are kept unchanged.
    """
    choices_per_position = []

    for codon in original_codons:

        codon = codon.upper()

        if codon not in CODON_TO_AA:
            raise ValueError(f"Unknown codon: {codon}")

        amino_acid = CODON_TO_AA[codon]

        # Preserve stop codons.
        if amino_acid == "*":
            choices_per_position.append([codon])
            continue

        if amino_acid not in codon_table:
            raise ValueError(
                f"Amino acid '{amino_acid}' "
                f"not found in codon table."
            )

        choices = list(codon_table[amino_acid].keys())

        if not choices:
            raise ValueError(
                f"No synonymous codons available for "
                f"amino acid '{amino_acid}'."
            )

        choices_per_position.append(choices)

    return choices_per_position


# ============================================================
# Population initialization
# ============================================================

def initialize_population(
    original_codons: List[str],
    synonymous: List[List[str]],
    population_size: int,
    rng: random.Random,
) -> List[List[str]]:
    """Create initial GA population."""

    if population_size < 1:
        raise ValueError(
            "population_size must be at least 1."
        )

    population = [
        original_codons.copy()
    ]

    while len(population) < population_size:

        individual = []

        for position, choices in enumerate(synonymous):

            original = original_codons[position]

            # Keep some original codons for diversity.
            if (
                original in choices
                and rng.random() < 0.30
            ):
                individual.append(original)
            else:
                individual.append(
                    rng.choice(choices)
                )

        population.append(individual)

    return population


# ============================================================
# Fitness
# ============================================================

def calculate_fitness(
    codons: List[str],
    original_codons: List[str],
    organism: str,
    gc_min: Optional[float] = None,
    gc_max: Optional[float] = None,
    cai_weight: float = 1.0,
    gc_weight: float = 10.0,
    edit_weight: float = 0.0,
) -> float:
    """
    Calculate GA fitness.

    F =
        CAI contribution
        - GC constraint penalty
        - edit-rate penalty
    """

    sequence = "".join(codons)

    cai = calculate_cai(
        sequence,
        organism,
    )

    gc = gc_content(sequence)

    violation = gc_penalty(
        gc,
        gc_min,
        gc_max,
    )

    changes = count_changes(
        original_codons,
        codons,
    )

    rate = (
        changes / len(original_codons)
        if original_codons
        else 0.0
    )

    # Strong penalty for violating global GC.
    HARD_CONSTRAINT_PENALTY = 1000.0

    return (
        cai_weight * cai
        - gc_weight
        * HARD_CONSTRAINT_PENALTY
        * violation
        - edit_weight * rate
    )


# ============================================================
# Selection
# ============================================================

def tournament_selection(
    population: List[List[str]],
    fitnesses: List[float],
    rng: random.Random,
    tournament_size: int = 3,
) -> List[str]:
    """Select an individual using tournament selection."""

    size = min(
        tournament_size,
        len(population),
    )

    indices = rng.sample(
        range(len(population)),
        size,
    )

    winner = max(
        indices,
        key=lambda index: fitnesses[index],
    )

    return population[winner].copy()


# ============================================================
# Crossover
# ============================================================

def crossover(
    parent1: List[str],
    parent2: List[str],
    rng: random.Random,
) -> Tuple[List[str], List[str]]:
    """Single-point crossover."""

    if len(parent1) != len(parent2):
        raise ValueError(
            "Parents must have equal length."
        )

    if len(parent1) < 2:
        return (
            parent1.copy(),
            parent2.copy(),
        )

    point = rng.randint(
        1,
        len(parent1) - 1,
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


# ============================================================
# Mutation
# ============================================================

def mutate(
    individual: List[str],
    synonymous: List[List[str]],
    mutation_rate: float,
    rng: random.Random,
) -> List[str]:
    """Replace codons with synonymous alternatives."""

    result = individual.copy()

    for position, choices in enumerate(synonymous):

        if rng.random() >= mutation_rate:
            continue

        alternatives = [
            codon
            for codon in choices
            if codon != result[position]
        ]

        if alternatives:
            result[position] = rng.choice(
                alternatives
            )

    return result


# ============================================================
# Genetic Algorithm
# ============================================================

def optimize_ga(
    dna_sequence: str,
    organism: str = "ecoli",

    # Global GC
    gc_min: Optional[float] = None,
    gc_max: Optional[float] = None,

    # Existing main.py parameters
    local_gc_window: Optional[int] = None,
    local_gc_min: Optional[float] = None,
    local_gc_max: Optional[float] = None,
    repeat_length: Optional[int] = None,
    max_homopolymer: Optional[int] = None,
    forbidden_sites: Optional[List[str]] = None,

    # GA parameters
    population_size: int = 100,
    generations: int = 100,
    mutation_rate: float = 0.02,
    crossover_rate: float = 0.80,
    tournament_size: int = 3,
    elitism: int = 2,

    # Fitness weights
    cai_weight: float = 1.0,
    gc_weight: float = 10.0,
    edit_weight: float = 0.0,

    seed: Optional[int] = None,
) -> str:
    """
    Optimize DNA using a Genetic Algorithm.

    The additional constraint parameters are accepted so that
    main.py can call the GA without a TypeError.

    They are currently reserved for the next constraint-wiring
    stage. Global GC is currently used directly in fitness.
    """

    start_time = time.perf_counter()

    # Prevent unused-parameter warnings while keeping API
    # compatible with main.py.
    _ = (
        local_gc_window,
        local_gc_min,
        local_gc_max,
        repeat_length,
        max_homopolymer,
        forbidden_sites,
    )

    # --------------------------------------------------------
    # Validate parameters
    # --------------------------------------------------------

    if population_size < 1:
        raise ValueError(
            "population_size must be at least 1."
        )

    if generations < 0:
        raise ValueError(
            "generations cannot be negative."
        )

    if not 0.0 <= mutation_rate <= 1.0:
        raise ValueError(
            "mutation_rate must be between 0 and 1."
        )

    if not 0.0 <= crossover_rate <= 1.0:
        raise ValueError(
            "crossover_rate must be between 0 and 1."
        )

    if gc_min is not None and gc_max is not None:

        if gc_min > gc_max:
            raise ValueError(
                "gc_min cannot be greater than gc_max."
            )

    # --------------------------------------------------------
    # Prepare sequence
    # --------------------------------------------------------

    dna_sequence = dna_sequence.upper().strip()

    original_codons = split_codons(
        dna_sequence
    )

    original_protein = translate_dna(
        dna_sequence
    )

    codon_table = get_codon_table(
        organism
    )

    synonymous = build_synonymous_codons(
        original_codons,
        codon_table,
    )

    rng = random.Random(seed)

    # --------------------------------------------------------
    # Initial population
    # --------------------------------------------------------

    population = initialize_population(
        original_codons,
        synonymous,
        population_size,
        rng,
    )

    # --------------------------------------------------------
    # Best solution
    # --------------------------------------------------------

    best_individual = original_codons.copy()

    best_fitness = calculate_fitness(
        best_individual,
        original_codons,
        organism,
        gc_min,
        gc_max,
        cai_weight,
        gc_weight,
        edit_weight,
    )

    # --------------------------------------------------------
    # Evolution
    # --------------------------------------------------------

    for _generation in range(generations):

        fitnesses = [
            calculate_fitness(
                individual,
                original_codons,
                organism,
                gc_min,
                gc_max,
                cai_weight,
                gc_weight,
                edit_weight,
            )
            for individual in population
        ]

        current_best_index = max(
            range(len(population)),
            key=lambda index: fitnesses[index],
        )

        current_best_fitness = fitnesses[
            current_best_index
        ]

        if current_best_fitness > best_fitness:

            best_fitness = current_best_fitness

            best_individual = population[
                current_best_index
            ].copy()

        # ----------------------------------------------------
        # Elitism
        # ----------------------------------------------------

        elite_count = min(
            max(elitism, 0),
            population_size,
        )

        ranked_indices = sorted(
            range(len(population)),
            key=lambda index: fitnesses[index],
            reverse=True,
        )

        new_population = [
            population[index].copy()
            for index in ranked_indices[:elite_count]
        ]

        # ----------------------------------------------------
        # Generate children
        # ----------------------------------------------------

        while len(new_population) < population_size:

            parent1 = tournament_selection(
                population,
                fitnesses,
                rng,
                tournament_size,
            )

            parent2 = tournament_selection(
                population,
                fitnesses,
                rng,
                tournament_size,
            )

            # Crossover
            if rng.random() < crossover_rate:

                child1, child2 = crossover(
                    parent1,
                    parent2,
                    rng,
                )

            else:

                child1 = parent1.copy()
                child2 = parent2.copy()

            # Mutation
            child1 = mutate(
                child1,
                synonymous,
                mutation_rate,
                rng,
            )

            child2 = mutate(
                child2,
                synonymous,
                mutation_rate,
                rng,
            )

            new_population.append(child1)

            if len(new_population) < population_size:
                new_population.append(child2)

        population = new_population

    # --------------------------------------------------------
    # Final population evaluation
    # --------------------------------------------------------

    final_fitnesses = [
        calculate_fitness(
            individual,
            original_codons,
            organism,
            gc_min,
            gc_max,
            cai_weight,
            gc_weight,
            edit_weight,
        )
        for individual in population
    ]

    final_best_index = max(
        range(len(population)),
        key=lambda index: final_fitnesses[index],
    )

    if final_fitnesses[final_best_index] > best_fitness:

        best_individual = population[
            final_best_index
        ].copy()

    optimized_sequence = "".join(
        best_individual
    )

    # --------------------------------------------------------
    # Protein verification
    # --------------------------------------------------------

    verification_result = verify_translation(
        dna_sequence,
        optimized_sequence,
    )

    # verify_translation() returns a dictionary.
    verification_passed = bool(
        verification_result["valid"]
    )

    if not verification_passed:
        raise RuntimeError(
            "GA produced a sequence that does not "
            "preserve the original protein."
        )

    # --------------------------------------------------------
    # GC feasibility
    # --------------------------------------------------------

    if not is_gc_feasible(
        optimized_sequence,
        gc_min,
        gc_max,
    ):

        feasible_candidates = [
            individual
            for individual in population
            if is_gc_feasible(
                "".join(individual),
                gc_min,
                gc_max,
            )
        ]

        if feasible_candidates:

            best_feasible = max(
                feasible_candidates,
                key=lambda individual:
                    calculate_fitness(
                        individual,
                        original_codons,
                        organism,
                        gc_min,
                        gc_max,
                        cai_weight,
                        gc_weight,
                        edit_weight,
                    ),
            )

            optimized_sequence = "".join(
                best_feasible
            )

        else:
            raise RuntimeError(
                "GA did not find a sequence satisfying "
                "the requested GC constraint."
            )

    # --------------------------------------------------------
    # Final verification
    # --------------------------------------------------------

    final_protein = translate_dna(
        optimized_sequence
    )

    if final_protein != original_protein:
        raise RuntimeError(
            "Protein sequence was not preserved."
        )

    _runtime = time.perf_counter() - start_time
    _ = _runtime

    # Existing tests expect a DNA string.
    return optimized_sequence


# ============================================================
# Compatibility alias
# ============================================================

def optimize_sequence_ga(
    dna_sequence: str,
    organism: str = "ecoli",
    **kwargs,
) -> str:
    """Compatibility wrapper."""
    return optimize_ga(
        dna_sequence,
        organism=organism,
        **kwargs,
    )