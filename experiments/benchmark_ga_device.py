from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from optimizers.genetic_algorithm import optimize_ga
from core.verification import verify_translation

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False


def main():
    fasta = ROOT / "data/genes/fasta/gene08_lacz.fasta"
    sequence = "".join(
        line.strip()
        for line in fasta.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()

    cuda_available = TORCH_AVAILABLE and torch.cuda.is_available()
    print("PyTorch available:", TORCH_AVAILABLE)
    print("GPU available:", cuda_available)
    if cuda_available:
        print("GPU:", torch.cuda.get_device_name(0))

    devices = ["cpu"]
    if cuda_available:
        devices.append("cuda")
    else:
        print("\nNote: CUDA not available in this environment. Benchmarking CPU only.")

    for device in devices:
        print(f"\nDevice: {device.upper()}")
        start = time.perf_counter()

        result = optimize_ga(
            dna_sequence=sequence,
            organism="ecoli",
            gc_min=0.40,
            gc_max=0.60,
            population_size=100,
            generations=100,
            seed=42,
            device=device,
        )

        elapsed = time.perf_counter() - start
        optimized = result["sequence"] if isinstance(result, dict) else result
        check = verify_translation(sequence, optimized)

        print(f"Runtime: {elapsed:.4f} seconds")
        print("Protein preserved:", bool(check["valid"]))
        print("Original length:", len(sequence))
        print("Optimized length:", len(optimized))


if __name__ == "__main__":
    main()
