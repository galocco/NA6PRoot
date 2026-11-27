#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath

parser = argparse.ArgumentParser()
parser.add_argument("run_directory", type=Path)
args = parser.parse_args()
folder = args.run_directory / "transitions"

with (folder / "populations.csv").open(newline="") as stream:
    populations = {(int(r["threshold"]), int(r["state"])): int(r["count"]) for r in csv.DictReader(stream)}
with (folder / "flows.csv").open(newline="") as stream:
    flows = [(int(r["from_threshold"]), int(r["to_threshold"]), int(r["from_state"]),
              int(r["to_state"]), int(r["count"])) for r in csv.DictReader(stream)]

labels = ["Unmatched", "Correct · 3 hits", "Correct · 4 hits", "Correct · 5 hits",
          "Wrong · 3 hits", "Wrong · 4 hits", "Wrong · 5 hits"]
colors = ["#8c9ba8", "#a2cda6", "#54a76e", "#167044", "#f2b8a5", "#dd795c", "#af4638"]
thresholds = (3, 4, 5)
total = sum(populations[3, state] for state in range(7))
if any(sum(populations[t, state] for state in range(7)) != total for t in thresholds):
    raise ValueError("Population totals differ")
for src, dst in zip(thresholds, thresholds[1:]):
    for state in range(7):
        if sum(n for x, y, a, b, n in flows if x == src and y == dst and a == state) != populations[src, state]:
            raise ValueError("Inconsistent outgoing flows")
        if sum(n for x, y, a, b, n in flows if x == src and y == dst and b == state) != populations[dst, state]:
            raise ValueError("Inconsistent incoming flows")

fig, ax = plt.subplots(figsize=(15, 8.8))
fig.subplots_adjust(left=.19, right=.81, top=.80, bottom=.11)
positions = {}
gap = total * .022
height = total + 6 * gap
for column, threshold in enumerate(thresholds):
    y = height
    for state in range(7):
        count = populations[threshold, state]
        if count == 0:
            continue
        y -= count
        positions[threshold, state] = (y, y + count)
        ax.add_patch(Rectangle((column - .025, y), .05, count, facecolor=colors[state], edgecolor="white", lw=.8, zorder=3))
        position = column - .045 if column == 0 else column + .045
        alignment = "right" if column == 0 else "left"
        ax.text(position, y + count/2, f"{labels[state]}  {count}", ha=alignment, va="center",
                fontsize=9.5, fontweight="medium", color="#263343", zorder=4)
        y -= gap

for left, right in zip(thresholds, thresholds[1:]):
    x0, x1 = thresholds.index(left) + .025, thresholds.index(right) - .025
    outgoing = {state: positions[left, state][0] for state in range(7) if (left, state) in positions}
    incoming = {state: positions[right, state][0] for state in range(7) if (right, state) in positions}
    selected = sorted((entry for entry in flows if entry[:2] == (left, right)), key=lambda e: (e[2], e[3]))
    for _, _, source, target, count in selected:
        y0, y1 = outgoing[source], incoming[target]
        outgoing[source] += count
        incoming[target] += count
        vertices = [(x0,y0),(x0+.42,y0),(x1-.42,y1),(x1,y1),
                    (x1,y1+count),(x1-.42,y1+count),(x0+.42,y0+count),(x0,y0+count),(x0,y0)]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
                 MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MplPath(vertices,codes), facecolor=colors[source], edgecolor="none", alpha=.34, zorder=1))

for column, threshold in enumerate(thresholds):
    ax.text(column, height+total*.07, f"VT threshold ≥ {threshold}", ha="center", va="bottom", fontsize=13, fontweight="bold")
ax.set_xlim(-.13, 2.13)
ax.set_ylim(-total*.015, height + total*.12)
ax.axis("off")
fig.text(.06,.94,"Thermal muon matching across VT thresholds", fontsize=20, fontweight="bold", color="#1f3041")
fig.text(.06,.902,f"{total:,} genuine thermal MS tracks with ≥6 hits · 1,000 simulated events", fontsize=11, color="#526474")
fig.text(.06,.045,"Each band follows the same MS track across the matching passes. Correct and wrong use MC truth.", fontsize=10, color="#526474")
fig.savefig(folder / "matching_transitions.png", dpi=190, facecolor="white")
