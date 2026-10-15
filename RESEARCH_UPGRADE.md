# Codon Optimization Research Upgrade & Implementation Architecture

## 1. System Overview

This project implements a rigorous, bioinformatics-grounded codon optimization and multi-objective sequence redesign system for synthetic biology and recombinant protein expression. The system redesigns coding DNA sequences (CDS) for target hosts (such as *Escherichia coli*) by selecting synonymous codons that enhance translation efficiency (Codon Adaptation Index, CAI) while strictly preserving the translated amino acid sequence and satisfying multi-scale biological constraints.

### Key Architectural Principle
**Decoupled Independent Validation**: An optimizer is never permitted to declare its own candidate feasible merely because its internal objective function or heuristic fitness score is high. Every optimizer—Greedy, MILP, single/multi-start Genetic Algorithm, and NSGA-II—outputs candidate sequences that are evaluated by the same independent verification and validation layer (`core/validator.py`). Feasibility is evaluated before ranking by CAI.

---

## 2. Directory & Module Responsibilities

| Module / File | Responsibility |
| :--- | :--- |
| `core/sequence.py` | DNA alphabet cleaning, validation, codon splitting, GC computation, reverse complement. |
| `core/verification.py` | Protein translation (`translate_dna`), synonymous verification (`verify_translation`, `protein_preserved`). |
| `core/validator.py` | **Centralized Independent Validation Layer**. Evaluates global GC, 5' GC, local sliding GC, homopolymers, repeats, and forbidden restriction sites. Emits booleans, quantitative measurements, and failure explanations. |
| `core/expression_prediction.py` | Feature extraction pipeline (CAI, GC, 5' GC, CBP, 6-mer frequency, homopolymers, length), baseline predictor, and regularized regression for proteomic/expression study. |
| `constraints/gc.py` | Global GC calculation, bounds checking, and violation penalties. |
| `constraints/local_gc.py` | Sliding-window GC profiling, range checks, and window violation identification. |
| `constraints/homopolymer.py` | Homopolymer run detection and violation counts. |
| `constraints/repeats.py` | Repeated $k$-mer subsequence detection ($k \ge 8, 10, 12$). |
| `constraints/restriction_sites.py` | Restriction site identification supporting palindromic and non-palindromic sites with reverse complement detection. Standard cloning enzyme catalog (`EcoRI`, `BamHI`, `HindIII`, `NdeI`, `XhoI`, `BsaI`, `NcoI`). |
| `optimizers/greedy.py` | Fast unconstrained baseline selecting the highest-frequency synonymous codon per position. |
| `optimizers/milp.py` | Mixed Integer Linear Programming optimizer (PuLP / CBC) supporting global GC bounds, 5' GC bounds, local GC windows, homopolymer cuts, and restriction-site elimination cuts. |
| `optimizers/genetic_algorithm.py` | Single/multi-start genetic algorithm optimizing CAI, global GC, target GC, local GC, 5' GC, homopolymers, repeats, restriction sites, and codon-pair bias. |
| `optimizers/nsga2.py` | **True Multi-Objective NSGA-II Optimizer**. Implements fast non-dominated sorting, crowding distance assignment, binary tournament selection, synonymous crossover/mutation, and elitist selection across 3 conflicting objectives: maximize CAI, minimize GC deviation, and minimize edit rate. |
| `optimizers/pareto.py` | Weighted-GA Pareto approximation, Pareto dominance helpers (`dominates`, `filter_non_dominated`), and unified interface. |
| `experiments/run_optimizer_comparison.py` | Fair comparative benchmark across Greedy, MILP, MILP_constrained, GA_single, GA_multistart, and NSGA-II. |
| `experiments/run_constraint_sensitivity.py` | Systematic evaluation of constraint stringency (repeats 8–12 bp, 5' GC bounds, local GC windows). |
| `experiments/run_fitness_weight_sensitivity.py` | Weight sensitivity analysis (CBP, 5' GC, $k$-mer complexity, combined weights). |
| `experiments/run_pareto_benchmark.py` | Comparative Pareto benchmark evaluating NSGA-II vs. Weighted-GA trade-offs. |
| `experiments/generate_plots.py` | Automated research figure generator creating all 8 publication-quality plots. |
| `experiments/run_research_suite.py` | Integrated master runner executing the complete research experiment suite. |
| `tests/` | 84 regression tests covering constraints, baselines, solvers, reproducibility, NSGA-II, and expression modeling. |

---

## 3. Mathematical Formulations & Constraints

### 3.1 Codon Adaptation Index (CAI)
For a gene $s$ composed of $L$ codons $(c_1, c_2, \dots, c_L)$:
$$\text{CAI}(s) = \left( \prod_{k=1}^L w_{c_k} \right)^{1/L} = \exp\left( \frac{1}{L} \sum_{k=1}^L \ln w_{c_k} \right)$$
where $w_{c_k} = \frac{f(c_k)}{\max_{c \in \text{syn}(AA_k)} f(c)}$ is the relative adaptiveness of codon $c_k$ within the target host codon usage table.

### 3.2 Global & 5' Region GC Content
$$GC_{\text{global}}(s) = \frac{\#G + \#C}{\text{len}(s)}, \quad GC_{5'}(s) = \frac{\#G + \#C \text{ in } s[0:w]}{w} \quad (w = 60 \text{ bp})$$
Feasible intervals:
$$0.40 \le GC_{\text{global}}(s) \le 0.60, \quad 0.40 \le GC_{5'}(s) \le 0.50$$
Lower 5' GC facilitates mRNA secondary structure unwinding at the ribosome binding site, promoting efficient translation initiation.

### 3.3 Sequence Constraints Conjunction
Feasibility is strictly binary:
$$\text{Feasible}(s) = \bigwedge_{i=1}^n \text{Constraint}_i(s)$$
Specifically:
1. $\text{ProteinPreserved}(s) = \text{True}$
2. $GC_{\text{global}}(s) \in [0.40, 0.60]$
3. $GC_{5'}(s) \in [0.40, 0.50]$ (first 60 bp)
4. $\forall \text{ window } w \text{ of 60 bp}: GC(w) \in [0.30, 0.70]$
5. $\text{LongestHomopolymer}(s) \le 6 \text{ bp}$
6. $\text{MaxRepeatedSubstring}(s) < 12 \text{ bp}$
7. $\text{RestrictionSites}(s) \cap \{\text{EcoRI}, \text{BamHI}, \text{HindIII}, \text{NdeI}, \text{XhoI}\} = \emptyset$

---

## 4. Multi-Objective Pareto Optimization (NSGA-II)

The NSGA-II module (`optimizers/nsga2.py`) formulates codon redesign as a true multi-objective discrete search:
- **Objective 1 (Maximize)**: $\text{CAI}(s)$ $\implies$ Minimize $f_1(s) = -\text{CAI}(s)$
- **Objective 2 (Minimize)**: $f_2(s) = |GC_{\text{global}}(s) - GC_{\text{target}}|$
- **Objective 3 (Minimize)**: $f_3(s) = \text{EditRate}(s_{\text{orig}}, s) = \frac{\sum_{k=1}^L \mathbb{I}[c_k \ne c_k^{\text{orig}}]}{L}$

### Constrained Domination Principle (Deb)
Solution $A$ constrained-dominates $B$ ($A \prec_C B$) if:
1. $A$ is feasible and $B$ is infeasible.
2. Both $A$ and $B$ are infeasible, but $\text{Violation}(A) < \text{Violation}(B)$.
3. Both $A$ and $B$ are feasible, and $A$ Pareto-dominates $B$:
   $$\forall m \in \{1, 2, 3\}: f_m(A) \le f_m(B) \quad \text{and} \quad \exists m \in \{1, 2, 3\}: f_m(A) < f_m(B)$$

### Elitist Selection & Crowding Distance
At each generation, parent population $P_t$ and offspring $Q_t$ ($|P_t| = |Q_t| = N$) are combined ($2N$), partitioned into non-dominated fronts $F_1, F_2, \dots$, and sorted using crowding distance:
$$I[i]_{\text{distance}} = \sum_{m=1}^M \frac{f_m(i+1) - f_m(i-1)}{f_m^{\max} - f_m^{\min} + \epsilon}$$
This preserves candidate diversity along the entire trade-off curve rather than collapsing into a single weight-biased point.

---

## 5. Experimental Verification & Execution

### 5.1 Verification Commands
```powershell
# 1. Verify dependencies
python -m pip install -r requirements.txt

# 2. Run full regression test suite (84 tests)
python -m pytest -q

# 3. Run constraint-aware multi-start GA smoke test
python -m experiments.run_constraint_aware_multistart --smoke

# 4. Run NSGA-II vs. Weighted-GA Pareto smoke test
python -m experiments.run_pareto_benchmark --smoke

# 5. Run full integrated research suite (smoke mode)
python -m experiments.run_research_suite --smoke

# 6. Generate all 8 research plots from CSV results
python -m experiments.generate_plots
```

All verification cycles pass cleanly without warnings or errors.
