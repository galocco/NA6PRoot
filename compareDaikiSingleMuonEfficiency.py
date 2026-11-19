"""Compare MS single-muon efficiencies using ms_metric_function's track selection.

Run: python3 compareDaikiSingleMuonEfficiency.py
Requires the optimizer's Python dependencies (numpy, uproot, matplotlib) and ROOT
development files used by its MC particle reader.
"""

import argparse
import importlib.util
import re
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_BASE = Path("debug_daiki/layout.medAbsorber[33]=Concrete")
DEFAULT_METRIC = Path("/home/galocco/na6precparam-optimizer/metrics/ms_metric_function.py")
DIRECTORY_RE = re.compile(r"reco\.(ms|vt)MaxDeltaTanLCellsCA\[0\]=(\d+)$")
AXES = {"p": ("p (GeV/c)", 0.0, 30.0),
        "pt": ("pT (GeV/c)", 0.0, 3.0),
        "y": ("y", 1.0, 7.0)}


def load_metric(path):
    if not path.is_file():
        raise FileNotFoundError(f"Metric module not found: {path}")
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location("ms_metric_function", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def metric_counts(module, selection):
    """Run the exact metric selection and capture its trackable/matched MC muons."""
    captured = {}

    def save_plots(_output, denominator, numerator, filename, _ylabel):
        if filename == "muon_efficiency_vs_p_pt_y.png":
            captured["den"] = denominator
            captured["num"] = numerator

    def mean_efficiency(_output, denominator, numerator, _filename,
                        n_pt_bins=10, n_y_bins=10):
        hist_range = [(0.0, 3.0), (1.0, 7.0)]
        den, _, _ = np.histogram2d(denominator["pt"], denominator["y"],
                                    bins=[n_pt_bins, n_y_bins], range=hist_range)
        num, _, _ = np.histogram2d(numerator["pt"], numerator["y"],
                                    bins=[n_pt_bins, n_y_bins], range=hist_range)
        filled = den > 0
        return float((num[filled] / den[filled]).mean()) if filled.any() else 0.0

    module._save_binned_plots = save_plots
    module._mean_bin_efficiency = mean_efficiency
    # The metric writes metrics.txt. Run it in scratch space with links to the
    # input ROOT files, leaving each selection's original output untouched.
    with tempfile.TemporaryDirectory(prefix="na6p_ms_eff_") as scratch:
        scratch = Path(scratch)
        for name in ("TracksMuonSpec.root", "ClustersMuonSpec.root", "MCKine.root"):
            source = selection / name
            if name == "MCKine.root" and not source.is_file():
                source = selection.parent / name
            if not source.is_file():
                raise FileNotFoundError(source)
            (scratch / name).symlink_to(source.resolve())
        module.metric_function(str(scratch))
    return captured["den"], captured["num"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--metric-module", type=Path, default=DEFAULT_METRIC)
    args = parser.parse_args()
    module = load_metric(args.metric_module.resolve())
    selections = []
    for path in args.base_dir.iterdir():
        match = DIRECTORY_RE.fullmatch(path.name)
        if not match or not path.is_dir():
            continue
        if not (path / "TracksMuonSpec.root").is_file() or not (path / "ClustersMuonSpec.root").is_file():
            print(f"Skipping {path}: raw MuonSpec track or cluster file is missing", file=sys.stderr)
            continue
        if not (path / "MCKine.root").is_file() and not (path.parent / "MCKine.root").is_file():
            print(f"Skipping {path}: MCKine.root is missing", file=sys.stderr)
            continue
        selections.append((match.group(1), int(match.group(2)), path))
    selections.sort(key=lambda item: (item[0], item[1]))
    if not selections:
        parser.error(f"No MS or VT cut selections found in {args.base_dir}")

    figures = {name: plt.subplots(figsize=(10, 7)) for name in AXES}
    for kind, cut, path in selections:
        denominator, numerator = metric_counts(module, path)
        label = f"{kind.upper()} CA tanL cut = {cut}"
        for name, (xlabel, low, high) in AXES.items():
            den, edges = np.histogram(denominator[name], bins=20, range=(low, high))
            num, _ = np.histogram(numerator[name], bins=20, range=(low, high))
            efficiency = np.divide(num, den, out=np.zeros(20), where=den > 0)
            error = np.sqrt(np.divide(efficiency * (1 - efficiency), den,
                                      out=np.zeros(20), where=den > 0))
            centers = (edges[:-1] + edges[1:]) / 2
            fig, ax = figures[name]
            ax.errorbar(centers[den > 0], efficiency[den > 0],
                        yerr=error[den > 0], fmt="o-", markersize=7,
                        capsize=2, label=label)
            ax.set(xlabel=xlabel, ylabel="Reconstructed / trackable", ylim=(0, 1.05),
                   title="Single muons: μ⁺ + μ⁻")
        print(f"{label}: {len(numerator['p'])} / {len(denominator['p'])} muons")

    for name, (fig, ax) in figures.items():
        ax.grid(axis="y", alpha=0.35)
        ax.legend()
        fig.tight_layout()
        output = args.base_dir / f"compare_single_muon_metric_efficiency_vs_{name}.png"
        fig.savefig(output)
        plt.close(fig)
        print(f"Saved {output}")


if __name__ == "__main__":
    main()
