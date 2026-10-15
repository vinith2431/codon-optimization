# Bioinformatics Codon Optimization & Multi-Objective Sequence Redesign System

A comprehensive bioinformatics framework for codon optimization and multi-objective coding sequence redesign for recombinant protein expression. The system redesigns coding DNA sequences for target hosts while strictly preserving the translated amino acid sequence and satisfying multi-scale biological constraints.

---

## Key Features

- **Decoupled Independent Validation**: Candidates from all optimizers are assessed by a centralized independent validation layer (`core/validator.py`). Optimizers cannot declare themselves valid.
- **Algorithms Compared Fairly**:
  - **Greedy Baseline**: Fast unconstrained selection of highest-frequency synonymous codons.
  - **MILP (PuLP / CBC)**: Mathematical integer programming supporting global GC bounds, 5' GC bounds, local GC windows, and linear cuts for restriction sites and homopolymers.
  - **Single-Start & Multi-Start GA**: Evolutionary search with configurable penalties for GC content, local windows, 5' initiation regions, homopolymer runs, repeated subsequences, and codon-pair bias.
  - **NSGA-II Multi-Objective Optimization**: True Pareto search exploring trade-offs across CAI maximization, GC target deviation minimization, and sequence edit rate minimization.
- **Biotechnology Sequence Constraints**:
  - Global GC content ($0.40 \le GC \le 0.60$)
  - 5' initiation region GC ($0.40 \le GC_{5'} \le 0.50$ over the first 60 bp)
  - Local sliding-window GC ($0.30 \le GC_{\text{local}} \le 0.70$ over 60 bp windows)
  - Homopolymer maximum length ($\le 6$ consecutive identical nucleotides)
  - Repeated subsequence prohibition (no repeated $\ge 12$-mers)
  - Cloning restriction site avoidance (EcoRI, BamHI, HindIII, NdeI, XhoI, BsaI)
- **Expression Prediction Framework**: Dedicated machine learning module evaluating feature extraction and linear vs. baseline models, distinguishing CAI from wet-lab experimental protein expression.
- **Automated Research Reporting**: End-to-end benchmark suite generating detailed per-run CSVs, aggregate summary tables, validation failure logs, and 8 publication-ready figures.

---

## Directory Structure

```text
├── core/
│   ├── sequence.py              # Sequence cleaning, validation, codon splitting, GC
│   ├── verification.py          # Translation and protein preservation checks
│   ├── validator.py             # Centralized independent validation layer
│   ├── cai.py                   # Codon Adaptation Index computation
│   ├── codon_pair.py            # Codon-pair frequency calculations
│   └── expression_prediction.py  # Feature extraction & expression modeling
├── constraints/
│   ├── gc.py                    # Global GC content bounds
│   ├── local_gc.py              # Sliding-window local GC profiling
│   ├── homopolymer.py           # Homopolymer run detection
│   ├── repeats.py               # Subsequence repeat identification
│   └── restriction_sites.py     # Restriction enzyme recognition & avoidance
├── optimizers/
│   ├── greedy.py                # Greedy baseline optimizer
│   ├── milp.py                  # MILP mathematical optimization
│   ├── genetic_algorithm.py     # Single and multi-start Genetic Algorithm
│   ├── nsga2.py                 # NSGA-II multi-objective Pareto optimizer
│   └── pareto.py                # Weighted-GA Pareto approximation & helpers
├── experiments/
│   ├── run_optimizer_comparison.py       # Fair comparative benchmark
│   ├── run_constraint_sensitivity.py     # Constraint stringency sensitivity study
│   ├── run_fitness_weight_sensitivity.py # Multi-objective weight sensitivity study
│   ├── run_pareto_benchmark.py           # NSGA-II vs. Weighted-GA Pareto benchmark
│   ├── run_constraint_aware_multistart.py# Multi-start GA feasibility benchmark
│   ├── run_research_suite.py             # Master integrated research suite
│   └── generate_plots.py                 # Automated research figure generation
├── data/
│   ├── codon_tables.py          # Host codon usage frequency tables
│   ├── codon_pair/              # Experimentally observed codon-pair tables
│   └── genes/fasta/             # 10 benchmark reference genes (GFP, Luciferase, etc.)
├── results/
│   ├── csv/                     # Output CSVs (summaries, details, failures)
│   └── plots/                   # 8 generated research plots
├── tests/                       # Complete regression test suite (84 tests)
├── requirements.txt             # Project dependencies
└── RESEARCH_UPGRADE.md          # Implementation details and architectural reference
```

---

## Installation & Setup

Ensure Python 3.9+ is installed. In PowerShell or your terminal, run:

```powershell
python -m pip install -r requirements.txt
```

Verify PuLP solver and CBC availability:

```powershell
python -c "import pulp; print('PuLP:', pulp.__version__); print('CBC solver:', pulp.PULP_CBC_CMD().available())"
```

---

## Verification & Testing

Execute the full regression suite:

```powershell
python -m pytest -q
```
All 84 tests should pass.

---

## Running Experiments

### 1. Constraint-Aware Multi-Start Benchmark
```powershell
# Fast smoke test
python -m experiments.run_constraint_aware_multistart --smoke

# Full benchmark across all 10 genes and 5 seeds
python -m experiments.run_constraint_aware_multistart
```

### 2. Pareto Multi-Objective Benchmark (NSGA-II)
```powershell
# Fast smoke test
python -m experiments.run_pareto_benchmark --smoke

# Full benchmark
python -m experiments.run_pareto_benchmark
```

### 3. Fair Optimizer Comparison
```powershell
python -m experiments.run_optimizer_comparison --smoke
```

### 4. Constraint & Fitness Weight Sensitivity Studies
```powershell
python -m experiments.run_constraint_sensitivity --smoke
python -m experiments.run_fitness_weight_sensitivity --smoke
```

### 5. Master Integrated Research Suite
Runs all benchmarks, aggregates failure records, and automatically outputs all 8 research figures:
```powershell
# Smoke mode
python -m experiments.run_research_suite --smoke

# Full research run
python -m experiments.run_research_suite
```

### 6. Plot Generation
Regenerate all 8 research figures at any time from the CSV results:
```powershell
python -m experiments.generate_plots
```
Plots are saved to `results/plots/`.
