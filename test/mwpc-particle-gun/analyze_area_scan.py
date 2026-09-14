#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


def parse_vec(line):
    return tuple(float(x) for x in line.split("=", 1)[1].strip().split(","))


def read_geometry(path):
    values = {}
    for raw in path.read_text().splitlines():
        if raw.startswith("gas_cm="):
            values["gas_size"] = parse_vec(raw)
        elif raw.startswith("gas_centre_cm="):
            values["gas_centre"] = parse_vec(raw)
    if "gas_size" not in values or "gas_centre" not in values:
        raise RuntimeError("geometry_summary.txt is missing gas size/centre")
    return values["gas_size"], values["gas_centre"]


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: analyze_area_scan.py RUN_DIRECTORY")

    run = Path(sys.argv[1]).resolve()
    positions = run / "gun_xy.csv"
    events = run / "events.csv"
    geometry = run / "geometry_summary.txt"
    if not positions.exists() or not events.exists() or not geometry.exists():
        raise RuntimeError("Run directory must contain gun_xy.csv, events.csv and geometry_summary.txt")

    gas_size, gas_centre = read_geometry(geometry)
    half_x = gas_size[0] / 2.0
    half_y = gas_size[1] / 2.0
    gas_x, gas_y = gas_centre[0], gas_centre[1]

    with positions.open(newline="") as f:
        pos = {int(r["event"]): r for r in csv.DictReader(f)}
    with events.open(newline="") as f:
        ev = {int(r["event"]): r for r in csv.DictReader(f)}

    if pos.keys() != ev.keys():
        raise RuntimeError("gun_xy.csv and events.csv contain different event sets")

    rows = []
    false_positive = 0
    false_negative = 0
    expected_inside = 0
    observed_inside = 0

    for event in sorted(pos):
        x = float(pos[event]["x_cm"])
        y = float(pos[event]["y_cm"])
        primary_entries = int(ev[event]["primary_entries"])
        observed = primary_entries > 0

        # The 71x71 default grid uses integer-centimetre coordinates. The gas
        # boundaries are at non-integer coordinates, so no test point sits
        # exactly on a boundary where navigator conventions could be ambiguous.
        expected = abs(x - gas_x) < half_x and abs(y - gas_y) < half_y

        expected_inside += int(expected)
        observed_inside += int(observed)
        if observed and not expected:
            false_positive += 1
        if expected and not observed:
            false_negative += 1

        rows.append({
            "event": event,
            "x_cm": x,
            "y_cm": y,
            "expected_in_gas": int(expected),
            "observed_primary_entry": int(observed),
            "primary_entries": primary_entries,
            "gas_steps": int(ev[event]["gas_steps"]),
            "energy_deposit_keV": float(ev[event]["energy_deposit_keV"]),
        })

    with (run / "area_scan.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    passed = false_positive == 0 and false_negative == 0
    summary = (
        f"gas_size_cm={gas_size[0]:.6g},{gas_size[1]:.6g},{gas_size[2]:.6g}\n"
        f"gas_centre_cm={gas_centre[0]:.6g},{gas_centre[1]:.6g},{gas_centre[2]:.6g}\n"
        f"grid_points={len(rows)}\n"
        f"expected_inside={expected_inside}\n"
        f"observed_primary_entries={observed_inside}\n"
        f"false_positive={false_positive}\n"
        f"false_negative={false_negative}\n"
        f"result={'PASS' if passed else 'FAIL'}\n"
    )
    (run / "area_scan_summary.txt").write_text(summary)
    print(summary, end="")

    if not passed:
        print("\nMismatched points (first 20):")
        shown = 0
        for row in rows:
            if row["expected_in_gas"] != row["observed_primary_entry"]:
                print(f"  event {row['event']}: x={row['x_cm']}, y={row['y_cm']}, "
                      f"expected={row['expected_in_gas']}, observed={row['observed_primary_entry']}")
                shown += 1
                if shown == 20:
                    break
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
