"""
NSGA-II (Non-dominated Sorting Genetic Algorithm II) for multi-objective codon optimization.

Objectives:
1. Maximize CAI (minimized as -CAI)
2. Minimize deviation from target GC (|actual_GC - target_GC|)
3. Minimize edit rate (fraction of altered codons from original sequence)

Constraints handled via Deb's constrained-dominance principle:
- Feasible solutions always dominate infeasible solutions.
- Among infeasible solutions, the one with smaller constraint violation dominates.
- Among feasible solutions, Pareto dominance across the 3 objectives is applied.
"""

from __future__ import annotations
import math
import random
from typing import Any, Dict, List, Optional, Tuple, Union

from core.cai import calculate_cai
from core.sequence import clean_dna, gc_content, get_codons
from core.similarity import edit_rate
from core.verification import CODON_TO_AA, translate_dna, verify_translation
from core.validator import assess_candidate, DEFAULT_CONSTRAINTS
from constraints.local_gc import local_gc_violations
from constraints.homopolymer import find_homopolymers
from constraints.repeats import find_repeats
from constraints.restriction_sites import find_restriction_sites
from data.codon_tables import get_codon_table


def build_synonymous_choices(original_codons: List[str], codon_table: dict) -> List[List[str]]:
    """Build list of synonymous codon choices for each position."""
    choices = []
    for codon in original_codons:
        aa = CODON_TO_AA[codon]
        if aa == "*":
            choices.append([codon])
        else:
            syn = list(codon_table[aa].keys())
            choices.append(syn if syn else [codon])
    return choices


class Individual:
    def __init__(self, codons: List[str]):
        self.codons = list(codons)
        self.sequence = "".join(self.codons)
        # Objectives (all formulated for minimization):
        # f1: -CAI
        # f2: |GC - target_GC|
        # f3: edit_rate
        self.objectives: List[float] = [0.0, 0.0, 0.0]
        self.cai: float = 0.0
        self.gc_dev: float = 0.0
        self.edit: float = 0.0
        # Constraint violation & assessment
        self.feasible: bool = False
        self.constraint_violation: float = 0.0
        self.assessment: Dict[str, Any] = {}
        # NSGA-II metadata
        self.rank: int = 0
        self.crowding_distance: float = 0.0
        self.domination_count: int = 0
        self.dominated_set: List[Individual] = []

    def evaluate(
        self,
        original_sequence: str,
        original_codons: List[str],
        organism: str,
        target_gc: float,
        constraints: Optional[Dict[str, Any]] = None,
    ):
        self.cai = calculate_cai(self.sequence, organism)
        actual_gc = gc_content(self.sequence)
        self.gc_dev = abs(actual_gc - target_gc)
        self.edit = edit_rate(original_sequence, self.sequence)

        # Objectives to minimize
        self.objectives = [-self.cai, self.gc_dev, self.edit]

        # Quantitative constraint violation calculation
        cfg = dict(DEFAULT_CONSTRAINTS)
        if constraints:
            cfg.update(constraints)

        violation = 0.0
        # Global GC violation
        gc_min = cfg.get("gc_min")
        gc_max = cfg.get("gc_max")
        if gc_min is not None and actual_gc < gc_min:
            violation += gc_min - actual_gc
        elif gc_max is not None and actual_gc > gc_max:
            violation += actual_gc - gc_max

        # 5' GC violation
        p5_min = cfg.get("five_prime_min")
        p5_max = cfg.get("five_prime_max")
        p5_win = cfg.get("five_prime_window", 60)
        if p5_min is not None or p5_max is not None:
            p5 = gc_content(self.sequence[: min(len(self.sequence), p5_win)])
            if p5_min is not None and p5 < p5_min:
                violation += p5_min - p5
            elif p5_max is not None and p5 > p5_max:
                violation += p5 - p5_max

        # Local GC violation (using optimized rolling window)
        loc_min = cfg.get("local_min")
        loc_max = cfg.get("local_max")
        loc_win = cfg.get("local_window", 60)
        if (loc_min is not None or loc_max is not None) and len(self.sequence) >= loc_win:
            v_windows = local_gc_violations(self.sequence, loc_win, loc_min, loc_max)
            if v_windows:
                violation += len(v_windows) * 0.05

        # Homopolymer violation
        max_hp = cfg.get("max_homopolymer")
        if max_hp is not None and max_hp > 0:
            hps = find_homopolymers(self.sequence, max_hp)
            if hps:
                violation += len(hps) * 0.1

        # Repeats violation
        rep_len = cfg.get("repeat_length")
        if rep_len is not None and rep_len > 0:
            reps = find_repeats(self.sequence, rep_len)
            if reps:
                violation += len(reps) * 0.1

        # Restriction site violation
        forbidden = cfg.get("forbidden_sites")
        if forbidden:
            check_rc = cfg.get("check_reverse_complement", True)
            sites = find_restriction_sites(self.sequence, forbidden, check_reverse_complement=check_rc)
            if sites:
                violation += len(sites) * 0.2

        # Independent validator is authoritative for feasibility.
        self.assessment = assess_candidate(
            original_sequence, self.sequence, constraints
        )

        # Penalize every failed validator check so infeasible candidates
        # cannot be treated as feasible during constrained selection.
        checks = (
            "protein_preserved",
            "global_gc_pass",
            "local_gc_pass",
            "five_prime_gc_pass",
            "homopolymer_pass",
            "repeat_pass",
            "restriction_pass",
        )
        failed_checks = sum(
            not self.assessment.get(check, False) for check in checks
        )
        violation += float(failed_checks)

        self.constraint_violation = violation
        self.feasible = bool(self.assessment["feasible"])


def constrained_dominates(ind_a: Individual, ind_b: Individual) -> bool:
    """
    Deb's Constrained Domination comparison.
    A dominates B if:
    1. A is feasible and B is infeasible, OR
    2. Both are infeasible, and A has smaller constraint violation, OR
    3. Both are feasible, and A Pareto-dominates B in objectives (<= on all, < on at least one).
    """
    if ind_a.feasible and not ind_b.feasible:
        return True
    if not ind_a.feasible and ind_b.feasible:
        return False
    if not ind_a.feasible and not ind_b.feasible:
        return ind_a.constraint_violation < ind_b.constraint_violation

    # Both are feasible: Pareto dominance (minimization)
    better_in_any = False
    for a_val, b_val in zip(ind_a.objectives, ind_b.objectives):
        if a_val > b_val:
            return False
        if a_val < b_val:
            better_in_any = True
    return better_in_any


def fast_non_dominated_sort(population: List[Individual]) -> List[List[Individual]]:
    """Partition population into non-dominated Pareto fronts F_1, F_2, ..."""
    fronts: List[List[Individual]] = [[]]

    for p in population:
        p.dominated_set = []
        p.domination_count = 0
        for q in population:
            if constrained_dominates(p, q):
                p.dominated_set.append(q)
            elif constrained_dominates(q, p):
                p.domination_count += 1

        if p.domination_count == 0:
            p.rank = 1
            fronts[0].append(p)

    i = 0
    while len(fronts[i]) > 0:
        next_front: List[Individual] = []
        for p in fronts[i]:
            for q in p.dominated_set:
                q.domination_count -= 1
                if q.domination_count == 0:
                    q.rank = i + 2
                    next_front.append(q)
        i += 1
        fronts.append(next_front)

    if len(fronts[-1]) == 0:
        fronts.pop()

    return fronts


def calculate_crowding_distance(front: List[Individual]):
    """Calculate crowding distance for individuals in a single Pareto front."""
    size = len(front)
    if size == 0:
        return
    for ind in front:
        ind.crowding_distance = 0.0

    if size <= 2:
        for ind in front:
            ind.crowding_distance = float("inf")
        return

    num_objectives = len(front[0].objectives)

    for m in range(num_objectives):
        front.sort(key=lambda ind: ind.objectives[m])
        front[0].crowding_distance = float("inf")
        front[-1].crowding_distance = float("inf")

        f_min = front[0].objectives[m]
        f_max = front[-1].objectives[m]
        norm = f_max - f_min

        if norm == 0.0:
            continue

        for i in range(1, size - 1):
            if front[i].crowding_distance != float("inf"):
                front[i].crowding_distance += (
                    (front[i + 1].objectives[m] - front[i - 1].objectives[m]) / norm
                )


def crowded_comparison(ind_a: Individual, ind_b: Individual) -> int:
    """
    Crowded comparison operator:
    Returns -1 if a is preferred to b, 1 if b is preferred to a, 0 if tied.
    """
    if ind_a.rank < ind_b.rank:
        return -1
    if ind_a.rank > ind_b.rank:
        return 1
    if ind_a.crowding_distance > ind_b.crowding_distance:
        return -1
    if ind_a.crowding_distance < ind_b.crowding_distance:
        return 1
    return 0


def tournament_select(population: List[Individual], rng: random.Random) -> Individual:
    """Binary tournament selection based on rank and crowding distance."""
    p1, p2 = rng.sample(population, 2)
    cmp = crowded_comparison(p1, p2)
    return p1 if cmp <= 0 else p2


def crossover(p1: Individual, p2: Individual, rng: random.Random) -> Tuple[Individual, Individual]:
    """Single-point synonymous crossover."""
    if len(p1.codons) < 2:
        return Individual(p1.codons), Individual(p2.codons)
    pt = rng.randint(1, len(p1.codons) - 1)
    c1 = p1.codons[:pt] + p2.codons[pt:]
    c2 = p2.codons[:pt] + p1.codons[pt:]
    return Individual(c1), Individual(c2)


def mutate(
    ind: Individual,
    synonymous: List[List[str]],
    mutation_rate: float,
    rng: random.Random,
) -> Individual:
    """Mutate individual via synonymous codon substitutions."""
    mutated = ind.codons.copy()
    for pos, choices in enumerate(synonymous):
        if rng.random() < mutation_rate and len(choices) > 1:
            alternatives = [c for c in choices if c != mutated[pos]]
            if alternatives:
                mutated[pos] = rng.choice(alternatives)
    return Individual(mutated)


def optimize_nsga2(
    sequence: str,
    organism: str = "ecoli",
    target_gc: float = 0.50,
    constraints: Optional[Dict[str, Any]] = None,
    population_size: int = 50,
    generations: int = 50,
    crossover_rate: float = 0.80,
    mutation_rate: float = 0.03,
    seed: Optional[int] = 42,
) -> List[Dict[str, Any]]:
    """
    Run NSGA-II multi-objective codon optimization.

    Returns
    -------
    list of dict
        Non-dominated Pareto candidates with:
        - sequence: DNA sequence
        - cai: Codon Adaptation Index
        - gc_deviation: Deviation from target GC
        - edit_rate: DNA edit rate
        - feasible: Independent validation feasibility
        - rank: Pareto rank (1 is non-dominated front)
        - crowding_distance: Diversity measure
        - failure_reasons: Reasons if infeasible
    """
    clean_seq = clean_dna(sequence)
    original_codons = get_codons(clean_seq)
    table = get_codon_table(organism)
    synonymous = build_synonymous_choices(original_codons, table)
    rng = random.Random(seed)

    # 1. Initialize population
    population: List[Individual] = [Individual(original_codons)]
    while len(population) < population_size:
        ind_codons = []
        for pos, choices in enumerate(synonymous):
            orig_c = original_codons[pos]
            if orig_c in choices and rng.random() < 0.25:
                ind_codons.append(orig_c)
            else:
                ind_codons.append(rng.choice(choices))
        population.append(Individual(ind_codons))

    # Evaluate initial population
    for ind in population:
        ind.evaluate(clean_seq, original_codons, organism, target_gc, constraints)

    # Initial ranking
    fronts = fast_non_dominated_sort(population)
    for f in fronts:
        calculate_crowding_distance(f)

    # 2. Main NSGA-II evolution loop
    for gen in range(1, generations + 1):
        offspring: List[Individual] = []

        while len(offspring) < population_size:
            p1 = tournament_select(population, rng)
            p2 = tournament_select(population, rng)

            if rng.random() < crossover_rate:
                c1, c2 = crossover(p1, p2, rng)
            else:
                c1, c2 = Individual(p1.codons), Individual(p2.codons)

            c1 = mutate(c1, synonymous, mutation_rate, rng)
            c2 = mutate(c2, synonymous, mutation_rate, rng)

            c1.evaluate(clean_seq, original_codons, organism, target_gc, constraints)
            c2.evaluate(clean_seq, original_codons, organism, target_gc, constraints)

            offspring.append(c1)
            if len(offspring) < population_size:
                offspring.append(c2)

        # 3. Elitist selection: Combine parents + offspring (size 2N -> N)
        combined = population + offspring
        fronts = fast_non_dominated_sort(combined)

        new_population: List[Individual] = []
        for front in fronts:
            calculate_crowding_distance(front)
            if len(new_population) + len(front) <= population_size:
                new_population.extend(front)
            else:
                # Sort front by crowding distance descending and take remainder
                front.sort(key=lambda ind: ind.crowding_distance, reverse=True)
                needed = population_size - len(new_population)
                new_population.extend(front[:needed])
                break

        population = new_population

    # 4. Final non-dominated sort of resulting population
    final_fronts = fast_non_dominated_sort(population)
    front_1 = final_fronts[0] if final_fronts else []
    calculate_crowding_distance(front_1)

    # Format result records, deduplicating identical sequences
    seen_seqs = set()
    results = []
    for ind in front_1:
        if ind.sequence in seen_seqs:
            continue
        seen_seqs.add(ind.sequence)
        ind.assessment = assess_candidate(clean_seq, ind.sequence, constraints)
        ind.feasible = ind.assessment["feasible"]
        results.append({
            "sequence": ind.sequence,
            "cai": round(ind.cai, 6),
            "gc_deviation": round(ind.gc_dev, 6),
            "gc_content": round(gc_content(ind.sequence), 6),
            "edit_rate": round(ind.edit, 6),
            "feasible": ind.feasible,
            "pareto_rank": ind.rank,
            "crowding_distance": (
                round(ind.crowding_distance, 4)
                if ind.crowding_distance != float("inf")
                else "inf"
            ),
            "failure_reasons": ind.assessment.get("failure_reasons", []),
            "seed": seed,
            "generations": generations,
        })

    # Sort candidates by CAI descending for ease of inspection
    results.sort(key=lambda x: x["cai"], reverse=True)
    return results
