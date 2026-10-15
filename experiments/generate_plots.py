"""
Automated Plot Generation for Codon Optimization Research Report.

Generates the 8 recommended research figures from saved experiment CSVs:
1. Feasibility rate by optimizer
2. Mean feasible CAI by optimizer
3. Runtime by optimizer
4. CAI versus edit rate
5. CAI versus GC deviation
6. Feasibility rate versus constraint threshold
7. CAI versus fitness weight
8. Pareto scatter plot showing non-dominated trade-offs
"""

from pathlib import Path
import csv
import sys
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CSV_DIR = ROOT / "results" / "csv"
PLOTS_DIR = ROOT / "results" / "plots"


def read_csv(path: Path) -> list:
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


def plot_optimizer_metrics():
    summary_path = CSV_DIR / "optimizer_summary.csv"
    detailed_path = CSV_DIR / "optimizer_comparison.csv"
    summary = read_csv(summary_path)
    detailed = read_csv(detailed_path)

    if not summary:
        print("optimizer_summary.csv not found; skipping optimizer summary plots.")
        return

    methods = [r["method"] for r in summary]
    feas_rates = [float(r["feasibility_rate_percent"]) for r in summary]
    runtimes = [float(r["mean_runtime_seconds"]) for r in summary]

    feasible_cais = []
    for r in summary:
        val = r.get("mean_feasible_cai", "None")
        feasible_cais.append(float(val) if val not in ("None", "") else 0.0)

    # 1. Feasibility rate by optimizer
    plt.figure(figsize=(8, 4.5))
    bars = plt.bar(methods, feas_rates, color="#2b5c8f", edgecolor="black")
    plt.title("Constraint Feasibility Rate by Optimizer", fontsize=13, fontweight="bold")
    plt.ylabel("Feasible Rate (%)", fontsize=11)
    plt.ylim(0, 105)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 2, f"{yval:.1f}%", ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "optimizer_feasibility_rate.png", dpi=200)
    plt.close()

    # 2. Mean feasible CAI by optimizer
    plt.figure(figsize=(8, 4.5))
    bars = plt.bar(methods, feasible_cais, color="#2e7d32", edgecolor="black")
    plt.title("Mean CAI Among Feasible Candidates", fontsize=13, fontweight="bold")
    plt.ylabel("Mean CAI", fontsize=11)
    plt.ylim(0, 1.1)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    for bar in bars:
        yval = bar.get_height()
        text = f"{yval:.3f}" if yval > 0 else "Infeasible"
        plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.02, text, ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "optimizer_mean_feasible_cai.png", dpi=200)
    plt.close()

    # 3. Runtime by optimizer
    plt.figure(figsize=(8, 4.5))
    bars = plt.bar(methods, runtimes, color="#b71c1c", edgecolor="black")
    plt.title("Computational Runtime by Optimizer", fontsize=13, fontweight="bold")
    plt.ylabel("Runtime (seconds)", fontsize=11)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2.0, yval, f"{yval:.3f}s", ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "optimizer_runtime.png", dpi=200)
    plt.close()

    # 4. CAI vs Edit Rate scatter & 5. CAI vs GC deviation scatter
    if detailed:
        valid = [
            r for r in detailed
            if r.get("optimized_cai") and r.get("edit_rate") and r.get("global_gc")
        ]
        if valid:
            cais = [float(r["optimized_cai"]) for r in valid]
            edits = [float(r["edit_rate"]) * 100 for r in valid]
            gc_devs = [abs(float(r["global_gc"]) - 0.50) for r in valid]
            feasible_mask = [r.get("feasible") == "True" for r in valid]

            # 4. CAI vs Edit Rate
            plt.figure(figsize=(7, 5))
            for feas, color, label in [(True, "#2e7d32", "Feasible"), (False, "#c62828", "Infeasible")]:
                x = [edits[i] for i in range(len(valid)) if feasible_mask[i] == feas]
                y = [cais[i] for i in range(len(valid)) if feasible_mask[i] == feas]
                if x:
                    plt.scatter(x, y, c=color, label=label, alpha=0.75, edgecolors="k", s=50)
            plt.title("Codon Adaptation Index vs. Sequence Edit Rate", fontsize=12, fontweight="bold")
            plt.xlabel("Edit Rate (% altered codons)", fontsize=10)
            plt.ylabel("Optimized CAI", fontsize=10)
            plt.legend()
            plt.grid(True, linestyle="--", alpha=0.5)
            plt.tight_layout()
            plt.savefig(PLOTS_DIR / "cai_vs_edit_rate.png", dpi=200)
            plt.close()

            # 5. CAI vs GC deviation
            plt.figure(figsize=(7, 5))
            for feas, color, label in [(True, "#2e7d32", "Feasible"), (False, "#c62828", "Infeasible")]:
                x = [gc_devs[i] for i in range(len(valid)) if feasible_mask[i] == feas]
                y = [cais[i] for i in range(len(valid)) if feasible_mask[i] == feas]
                if x:
                    plt.scatter(x, y, c=color, label=label, alpha=0.75, edgecolors="k", s=50)
            plt.title("Codon Adaptation Index vs. GC Deviation (|GC - 0.50|)", fontsize=12, fontweight="bold")
            plt.xlabel("Absolute GC Deviation", fontsize=10)
            plt.ylabel("Optimized CAI", fontsize=10)
            plt.legend()
            plt.grid(True, linestyle="--", alpha=0.5)
            plt.tight_layout()
            plt.savefig(PLOTS_DIR / "cai_vs_gc_deviation.png", dpi=200)
            plt.close()


def plot_constraint_sensitivity():
    summary_path = CSV_DIR / "constraint_summary.csv"
    summary = read_csv(summary_path)
    if not summary:
        return

    # 6. Feasibility rate versus constraint threshold
    cfgs = [r["config_name"] for r in summary]
    feas_rates = [float(r["feasibility_rate_percent"]) for r in summary]

    plt.figure(figsize=(9, 4.5))
    bars = plt.bar(cfgs, feas_rates, color="#4a148c", edgecolor="black")
    plt.title("Candidate Feasibility Under Varying Sequence Constraints", fontsize=13, fontweight="bold")
    plt.ylabel("Feasible Candidate Rate (%)", fontsize=11)
    plt.xticks(rotation=25, ha="right", fontsize=9)
    plt.ylim(0, 105)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 2, f"{yval:.1f}%", ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "feasibility_vs_constraint_threshold.png", dpi=200)
    plt.close()


def plot_weight_sensitivity():
    summary_path = CSV_DIR / "fitness_weight_summary.csv"
    summary = read_csv(summary_path)
    if not summary:
        return

    # 7. CAI versus fitness weight
    cfgs = [r["config_name"] for r in summary]
    cais = []
    for r in summary:
        val = r.get("mean_feasible_cai", "None")
        cais.append(float(val) if val not in ("None", "") else 0.0)

    plt.figure(figsize=(9, 4.5))
    bars = plt.bar(cfgs, cais, color="#00695c", edgecolor="black")
    plt.title("Mean Feasible CAI Across Fitness Weight Configurations", fontsize=13, fontweight="bold")
    plt.ylabel("Mean Feasible CAI", fontsize=11)
    plt.xticks(rotation=25, ha="right", fontsize=9)
    plt.ylim(0, 1.0)
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    for bar in bars:
        yval = bar.get_height()
        txt = f"{yval:.3f}" if yval > 0 else "N/A"
        plt.text(bar.get_x() + bar.get_width() / 2.0, yval + 0.02, txt, ha="center", va="bottom", fontsize=8)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "cai_vs_fitness_weight.png", dpi=200)
    plt.close()


def plot_pareto_scatter():
    pareto_path = CSV_DIR / "pareto_candidates.csv"
    candidates = read_csv(pareto_path)
    if not candidates:
        return

    # 8. Pareto scatter plot showing non-dominated candidates
    plt.figure(figsize=(8, 5.5))
    methods = list(dict.fromkeys(r["method"] for r in candidates))
    colors = {"NSGA-II": "#1565c0", "Weighted_GA": "#e65100"}

    for m in methods:
        sub = [r for r in candidates if r["method"] == m]
        cais = [float(r["cai"]) for r in sub]
        edits = [float(r["edit_rate"]) * 100 for r in sub]
        feasible = [r.get("feasible") == "True" for r in sub]

        markers = ["o" if f else "x" for f in feasible]
        for i in range(len(sub)):
            plt.scatter(
                edits[i],
                cais[i],
                color=colors.get(m, "#333333"),
                marker=markers[i],
                alpha=0.8,
                s=70,
                label=f"{m} ({'Feasible' if feasible[i] else 'Infeasible'})"
                if i == 0 or (i == 1 and feasible[0] != feasible[1])
                else "",
            )

    plt.title("Pareto Non-Dominated Front: CAI vs. Edit Rate Trade-Off", fontsize=13, fontweight="bold")
    plt.xlabel("Edit Rate (% altered codons) [Minimize]", fontsize=11)
    plt.ylabel("Codon Adaptation Index (CAI) [Maximize]", fontsize=11)
    plt.grid(True, linestyle="--", alpha=0.5)

    # Clean legend handles
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc="lower right")
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "pareto_tradeoffs_scatter.png", dpi=200)
    plt.close()


def main():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    plot_optimizer_metrics()
    plot_constraint_sensitivity()
    plot_weight_sensitivity()
    plot_pareto_scatter()
    print(f"Generated research plots saved to: {PLOTS_DIR}")


if __name__ == "__main__":
    main()
