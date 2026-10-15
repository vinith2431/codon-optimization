"""
Tests for Expression Prediction framework.
"""

from core.expression_prediction import (
    extract_sequence_features,
    BaselineMeanModel,
    RidgeRegressionModel,
    evaluate_predictions,
)


def test_extract_sequence_features():
    seq = "ATGGCCGCGTTCTGA"
    feats = extract_sequence_features(seq, "ecoli")
    assert "cai" in feats
    assert "global_gc" in feats
    assert "five_prime_gc" in feats
    assert "max_6mer_freq" in feats
    assert "longest_homopolymer" in feats
    assert "sequence_length" in feats
    assert feats["sequence_length"] == 15.0


def test_baseline_and_ridge_models():
    # Synthetic dataset
    X = [[1.0, 0.5], [2.0, 0.6], [3.0, 0.7], [4.0, 0.8]]
    y = [10.0, 20.0, 30.0, 40.0]

    # Baseline mean
    base = BaselineMeanModel()
    base.fit(X, y)
    preds_base = base.predict(X)
    assert all(p == 25.0 for p in preds_base)
    metrics_base = evaluate_predictions(y, preds_base)
    assert metrics_base["r2"] == 0.0

    # Ridge
    ridge = RidgeRegressionModel(alpha=0.1)
    ridge.fit(X, y)
    preds_ridge = ridge.predict(X)
    metrics_ridge = evaluate_predictions(y, preds_ridge)
    # Ridge should achieve lower RMSE than baseline on correlated linear data
    assert metrics_ridge["rmse"] < metrics_base["rmse"]
