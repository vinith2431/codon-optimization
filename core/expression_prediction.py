"""
Protein Expression Prediction Framework (Research Extension).

IMPORTANT RESEARCH DISTINCTION:
Codon Adaptation Index (CAI) measures codon usage conformity to a host reference table.
CAI is NOT an experimental measurement of protein expression.
Direct prediction of experimental protein expression requires labeled wet-lab data
(e.g., mass-spectrometry protein abundance, GFP fluorescence arbitrary units, Western blot).

This module provides:
1. Sequence feature extraction for expression modeling
2. Baseline mean predictor
3. Regularized linear regression model (Ridge)
4. Evaluation metrics (MAE, RMSE, R^2) comparing trained models against baseline
"""

from typing import Any, Dict, List, Optional, Tuple
import math
from collections import Counter
from core.cai import calculate_cai
from core.sequence import clean_dna, gc_content
from core.validator import longest_homopolymer_run


def extract_sequence_features(
    sequence: str,
    organism: str = "ecoli",
    codon_pair_table: Optional[dict] = None,
) -> Dict[str, float]:
    """
    Extract multi-scale biological features from a coding DNA sequence.

    Features:
    - cai: Codon Adaptation Index
    - global_gc: Global GC content
    - five_prime_gc: GC content of the first 60 nucleotides
    - max_6mer_freq: Maximum frequency of any 6-mer (sequence complexity)
    - longest_homopolymer: Maximum consecutive run of identical bases
    - sequence_length: Nucleotide length
    - cbp_score: Codon-pair bias score (if table provided)
    """
    seq = clean_dna(sequence)
    length = len(seq)

    # 1. CAI
    cai_val = calculate_cai(seq, organism)

    # 2. Global GC
    gc_val = gc_content(seq)

    # 3. 5' GC (first 60 nt)
    p5_val = gc_content(seq[: min(60, length)])

    # 4. K-mer frequency (k=6)
    if length >= 6:
        kmers = [seq[i : i + 6] for i in range(length - 5)]
        counts = Counter(kmers)
        max_kmer_freq = max(counts.values()) / len(kmers)
    else:
        max_kmer_freq = 0.0

    # 5. Longest homopolymer
    longest_hp = float(longest_homopolymer_run(seq))

    # 6. Codon-pair score
    cbp_score = 0.0
    if codon_pair_table and length >= 6:
        codons = [seq[i : i + 3] for i in range(0, length, 3)]
        pairs = [(codons[i], codons[i + 1]) for i in range(len(codons) - 1)]
        freqs = [
            float(codon_pair_table.get(p, codon_pair_table.get(f"{p[0]}-{p[1]}", 0.0)))
            for p in pairs
        ]
        if freqs:
            cbp_score = sum(freqs) / len(freqs)

    return {
        "cai": cai_val,
        "global_gc": gc_val,
        "five_prime_gc": p5_val,
        "max_6mer_freq": max_kmer_freq,
        "longest_homopolymer": longest_hp,
        "sequence_length": float(length),
        "cbp_score": cbp_score,
    }


class BaselineMeanModel:
    """Baseline model predicting the training set mean expression."""

    def __init__(self):
        self.mean_val: float = 0.0

    def fit(self, X: List[List[float]], y: List[float]):
        self.mean_val = sum(y) / len(y) if y else 0.0

    def predict(self, X: List[List[float]]) -> List[float]:
        return [self.mean_val for _ in X]


class RidgeRegressionModel:
    """Simple Ridge regression implemented via normal equations (X^T X + alpha I)^-1 X^T y."""

    def __init__(self, alpha: float = 1.0):
        self.alpha = alpha
        self.weights: List[float] = []
        self.intercept: float = 0.0
        self.mean_X: List[float] = []
        self.std_X: List[float] = []

    def fit(self, X: List[List[float]], y: List[float]):
        n_samples = len(X)
        if n_samples == 0:
            return
        n_features = len(X[0])

        # Standardize features
        self.mean_X = [sum(X[i][j] for i in range(n_samples)) / n_samples for j in range(n_features)]
        self.std_X = []
        for j in range(n_features):
            var = sum((X[i][j] - self.mean_X[j]) ** 2 for i in range(n_samples)) / n_samples
            self.std_X.append(math.sqrt(var) if var > 1e-9 else 1.0)

        mean_y = sum(y) / n_samples
        y_centered = [yi - mean_y for yi in y]

        X_norm = [
            [(X[i][j] - self.mean_X[j]) / self.std_X[j] for j in range(n_features)]
            for i in range(n_samples)
        ]

        # Solve (X^T X + alpha I) w = X^T y using gradient descent or direct solver
        # For small n_features (7 features), gradient descent with momentum converges rapidly
        w = [0.0] * n_features
        lr = 0.01
        for epoch in range(1000):
            # Compute predictions
            preds = [sum(X_norm[i][j] * w[j] for j in range(n_features)) for i in range(n_samples)]
            # Gradients
            for j in range(n_features):
                grad = (
                    -2.0
                    * sum(X_norm[i][j] * (y_centered[i] - preds[i]) for i in range(n_samples))
                    / n_samples
                    + 2.0 * self.alpha * w[j]
                )
                w[j] -= lr * grad

        self.weights = w
        self.intercept = mean_y

    def predict(self, X: List[List[float]]) -> List[float]:
        preds = []
        for row in X:
            norm_row = [(row[j] - self.mean_X[j]) / self.std_X[j] for j in range(len(row))]
            pred = self.intercept + sum(norm_row[j] * self.weights[j] for j in range(len(row)))
            preds.append(pred)
        return preds


def evaluate_predictions(y_true: List[float], y_pred: List[float]) -> Dict[str, float]:
    """Calculate MAE, RMSE, and R^2 score."""
    if not y_true or len(y_true) != len(y_pred):
        return {"mae": 0.0, "rmse": 0.0, "r2": 0.0}

    n = len(y_true)
    mae = sum(abs(t - p) for t, p in zip(y_true, y_pred)) / n
    rmse = math.sqrt(sum((t - p) ** 2 for t, p in zip(y_true, y_pred)) / n)

    mean_t = sum(y_true) / n
    ss_tot = sum((t - mean_t) ** 2 for t in y_true)
    ss_res = sum((t - p) ** 2 for t, p in zip(y_true, y_pred))

    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-9 else 0.0

    return {
        "mae": round(mae, 4),
        "rmse": round(rmse, 4),
        "r2": round(r2, 4),
    }
