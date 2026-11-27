"""Overlay the projected tungsten plug edges on the MS energy bias plots."""

import argparse
import configparser
import csv
import math
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def plug_sections(layout, z_vertex):
    sections = []
    z_front = layout.getfloat("posZStartAbsorber")
    n_slices = layout.getint("nAbsorberSlices") - layout.getint("nAbsorberSlicesWall")
    for i in range(n_slices):
        z_back = z_front + layout.getfloat(f"thicknessAbsorber[{i}]")
        if layout[f"medAbsorberPlug[{i}]"] == "Tungsten":
            r_out = layout.getfloat(f"radPlug[{i}]")
            r_in = layout.getfloat(f"radPlugInt[{i}]")
            if r_out > 0:
                sections.append({
                    "slice": i,
                    "z_front_cm": z_front,
                    "z_back_cm": z_back,
                    "r_out_cm": r_out,
                    "r_in_cm": r_in,
                    "eta_outer_front": math.asinh((z_front - z_vertex) / r_out),
                    "eta_outer_back": math.asinh((z_back - z_vertex) / r_out),
                    "eta_inner_back": math.asinh((z_back - z_vertex) / r_in) if r_in > 0 else "",
                })
        z_front = z_back
    return sections


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--rapidity-width", type=float, default=0.1)
    parser.add_argument("--vertex-z", type=float, help="Reference vertex z in cm. Default is the mean target position.")
    panels = parser.add_mutually_exclusive_group()
    panels.add_argument("--two-panels", action="store_true", help="Draw two full-range panels with a shared legend above them.")
    panels.add_argument("--split-thermal", action="store_true", help="Draw thermal and other muons in separate full-range rows.")
    args = parser.parse_args()

    config = configparser.ConfigParser(interpolation=None)
    with (args.run / "reco" / "na6pLayout.ini").open() as source:
        config.read_file(source)
    layout = config["layout"]
    z_vertex = args.vertex_z
    if z_vertex is None:
        z_vertex = mean(layout.getfloat(f"posTargetZ[{i}]") for i in range(layout.getint("nTargets")))
    sections = plug_sections(layout, z_vertex)
    if not sections:
        raise ValueError("No tungsten plug sections in the layout")

    out = args.run / "ms_muon_energy_resolution"
    suffix = "" if args.rapidity_width == 0.2 else f"_dy{args.rapidity_width:g}"
    stem = f"bias_vs_rapidity{suffix}"
    def read_rows(directory):
        with (directory / f"{stem}.csv").open() as source:
            return [r for r in csv.DictReader(source) if int(r["count"]) >= 10]

    rows = read_rows(out)
    if not rows:
        raise ValueError("No rapidity bins with at least 10 tracks")

    edges = [section["eta_outer_back"] for section in sections]
    inner_edges = [section["eta_inner_back"] for section in sections if section["r_in_cm"] > 0]
    zoom_lo, zoom_hi = min(edges) - 0.3, max(edges) + 0.45
    edge_color = "#6b4c9a"
    if args.split_thermal:
        row_samples = [read_rows(out / sample) for sample in ("thermal", "other_muons")]
        if not all(row_samples):
            raise ValueError("Missing populated rapidity bins for one of the muon samples")
        labels_by_row = ("Thermal muons", "Other muons (including secondaries)")
        plot_rows = [r for sample in row_samples for r in sample]
        full_range = (
            min(float(r["y_low"]) for r in plot_rows),
            max(float(r["y_high"]) for r in plot_rows),
        )
        inner_edges = [edge for edge in inner_edges if full_range[0] <= edge <= full_range[1]]
        fig, axes = plt.subplots(2, 2, figsize=(12, 8.2), sharex=True)
        fig.subplots_adjust(left=.09, right=.98, bottom=.1, top=.84, hspace=.33, wspace=.26)
        zoom_flags = (False, False)
    elif args.two_panels:
        data_lo = min(float(r["y_low"]) for r in rows)
        data_hi = max(float(r["y_high"]) for r in rows)
        full_range = (data_lo - .02, data_hi + .02)
        inner_edges = [edge for edge in inner_edges if full_range[0] <= edge <= full_range[1]]
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), squeeze=False)
        fig.subplots_adjust(left=.09, right=.98, bottom=.16, top=.80, wspace=.26)
        zoom_flags = (False,)
    else:
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
        zoom_flags = (False, True)
    for row_index, (row, zoom) in enumerate(zip(axes, zoom_flags)):
        row_data = row_samples[row_index] if args.split_thermal else rows
        for ax, unit, ylabel in zip(row, ("MeV", "pct"), (
            r"Mean $E_{\mathrm{gen}}-E_{\mathrm{rec}}$ (MeV)",
            r"Mean $(E_{\mathrm{gen}}-E_{\mathrm{rec}})/E_{\mathrm{gen}}$ (%)",
        )):
            for edge in edges:
                ax.axvline(edge, color=edge_color, alpha=0.42, linewidth=0.7, zorder=1)
            for edge in inner_edges:
                ax.axvline(edge, color=edge_color, linestyle="--", alpha=0.75, linewidth=0.9, zorder=1)
            ax.axhline(0, color="0.45", linestyle="--", linewidth=0.8, zorder=2)
            for model, color, marker in (
                ("Ionization", "#0072b2", "o"),
                ("Ionization + radiative", "#b2182b", "s"),
            ):
                sample = [r for r in row_data if r["model"] == model]
                if zoom:
                    sample = [r for r in sample if zoom_lo <= float(r["y_mean"]) <= zoom_hi]
                ax.errorbar(
                    [float(r["y_mean"]) for r in sample],
                    [float(r[f"mean_{unit}"]) for r in sample],
                    yerr=[float(r[f"mean_sem_{unit}"]) for r in sample],
                    marker=marker, markersize=4, linestyle="none", capsize=2,
                    color=color, label=model, zorder=3,
                )
            if zoom:
                ax.set_xlim(zoom_lo, zoom_hi)
            elif args.two_panels or args.split_thermal:
                ax.set_xlim(*full_range)
            ax.set(xlabel=r"Generated muon rapidity $y_{\mathrm{gen}}$", ylabel=ylabel)
            if args.split_thermal:
                ax.set_title(labels_by_row[row_index])
                if row_index == 0:
                    ax.set_xlabel("")
            ax.grid(axis="y", alpha=0.26)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    handles.append(Line2D([0], [0], color=edge_color, linewidth=1.5))
    labels.append(r"W slice outer edge at exit (projected $\eta$)")
    if inner_edges:
        handles.append(Line2D([0], [0], color=edge_color, linestyle="--", linewidth=1.2))
        labels.append(r"Annular slice inner edge at exit (projected $\eta$)")
    fig.suptitle(f"MS energy bias at the primary vertex  |  plug projection from z = {z_vertex:g} cm")
    plot_path = out / f"{stem}_with_W_plug"
    if args.two_panels or args.split_thermal:
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, .91), ncol=len(labels), frameon=False)
        fig.savefig(f"{plot_path}.pdf")
    else:
        axes[0, 1].legend(handles, labels, loc="lower left", fontsize=8)
    fig.savefig(f"{plot_path}.png", dpi=300 if args.two_panels or args.split_thermal else 180)

    with (out / "W_plug_rapidity_edges.csv").open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=["vertex_z_cm", *sections[0]])
        writer.writeheader()
        for section in sections:
            writer.writerow({"vertex_z_cm": z_vertex, **section})


if __name__ == "__main__":
    main()
