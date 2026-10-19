#!/usr/bin/env python3
"""
SARS-CoV-2 Deletion Variants: Genomic Coordinate 2D Scatter Plot
Generates the exact Start Position vs Stop Position deletion map as presented
in Figure 1 and Figure 7B of Jiang, Dewey, & Yin (2025):
"Deletion detection in SARS-CoV-2 genomes using multiplex-PCR sequencing..."
(https://www.medrxiv.org/content/10.1101/2025.04.15.25325794v2)

Usage:
    python3 plot_deletions.py <input_file_or_dir> [output_image]
"""

import sys
import os
import glob
import math

import csv

def generate_svg(records, output_svg, max_coord=31500):
    width, height = 700, 700
    margin = 85
    plot_w = width - 2 * margin
    plot_h = height - 2 * margin

    def map_x(coord):
        return margin + (coord / max_coord) * plot_w

    def map_y(coord):
        # Y=0 is bottom
        return margin + plot_h - (coord / max_coord) * plot_h

    svg_elements = []
    # Background and frame
    svg_elements.append(f'<rect width="{width}" height="{height}" fill="#ffffff"/>')
    svg_elements.append(f'<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="#fcfcfc" stroke="#222222" stroke-width="1.5"/>')

    # Grid / diagonal line y = x (boundary where start = end)
    diag_x1, diag_y1 = map_x(0), map_y(0)
    diag_x2, diag_y2 = map_x(max_coord), map_y(max_coord)
    svg_elements.append(f'<line x1="{diag_x1:.1f}" y1="{diag_y1:.1f}" x2="{diag_x2:.1f}" y2="{diag_y2:.1f}" stroke="#aaaaaa" stroke-dasharray="4,4" stroke-width="1.2"/>')

    # Axes ticks and numbers (every 5000 nt)
    for tick in range(0, max_coord + 1, 5000):
        tx = map_x(tick)
        ty = map_y(tick)
        # X axis ticks
        svg_elements.append(f'<line x1="{tx:.1f}" y1="{margin + plot_h}" x2="{tx:.1f}" y2="{margin + plot_h + 6}" stroke="#222222" stroke-width="1.2"/>')
        svg_elements.append(f'<text x="{tx:.1f}" y="{margin + plot_h + 22}" font-family="sans-serif" font-size="12" text-anchor="middle" fill="#333333">{tick}</text>')
        # Y axis ticks
        svg_elements.append(f'<line x1="{margin - 6}" y1="{ty:.1f}" x2="{margin}" y2="{ty:.1f}" stroke="#222222" stroke-width="1.2"/>')
        svg_elements.append(f'<text x="{margin - 10}" y="{ty + 4:.1f}" font-family="sans-serif" font-size="12" text-anchor="end" fill="#333333">{tick}</text>')

    # Axis Labels
    svg_elements.append(f'<text x="{margin + plot_w/2:.1f}" y="{margin + plot_h + 55}" font-family="sans-serif" font-size="15" font-weight="bold" text-anchor="middle" fill="#111111">Stop Position (nt)</text>')
    svg_elements.append(f'<text x="{margin - 55}" y="{margin + plot_h/2:.1f}" font-family="sans-serif" font-size="15" font-weight="bold" text-anchor="middle" fill="#111111" transform="rotate(-90, {margin - 55}, {margin + plot_h/2:.1f})">Start Position (nt)</text>')

    # Title
    svg_elements.append(f'<text x="{margin + plot_w/2:.1f}" y="{margin - 40}" font-family="sans-serif" font-size="16" font-weight="bold" text-anchor="middle" fill="#111111">SARS-CoV-2 Deletion Variants (Figure 1 Reproduction)</text>')
    svg_elements.append(f'<text x="{margin + plot_w/2:.1f}" y="{margin - 18}" font-family="sans-serif" font-size="12" text-anchor="middle" fill="#666666">Mapped Junctions (n = {len(records)})</text>')

    # Draw deletion variant points
    for r in records:
        cx = map_x(r["end"])
        cy = map_y(r["start"])
        alpha = r["alpha"]
        svg_elements.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="3" fill="#0033cc" fill-opacity="{alpha:.2f}"/>')

    # Legend (log10 Frequency opacity guide)
    leg_x = margin + 30
    leg_y = margin + 35
    svg_elements.append(f'<rect x="{leg_x - 10}" y="{leg_y - 15}" width="165" height="135" fill="#ffffff" fill-opacity="0.9" stroke="#cccccc" stroke-width="1"/>')
    svg_elements.append(f'<text x="{leg_x}" y="{leg_y}" font-family="sans-serif" font-size="12" font-weight="bold" fill="#222222">log10(Frequency)</text>')
    legend_breaks = [(-2.0, 0.05, "-2.0 (0.01)"), (-1.5, 0.28, "-1.5"), (-1.0, 0.52, "-1.0 (0.10)"), (-0.5, 0.76, "-0.5"), (0.0, 1.0, " 0.0 (1.00)")]
    for i, (lf, a, lbl) in enumerate(legend_breaks):
        ly = leg_y + 20 + i * 20
        svg_elements.append(f'<circle cx="{leg_x + 10}" cy="{ly}" r="4" fill="#0033cc" fill-opacity="{a:.2f}"/>')
        svg_elements.append(f'<text x="{leg_x + 28}" y="{ly + 4}" font-family="sans-serif" font-size="11" fill="#444444">{lbl}</text>')

    svg_header = f'<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
    svg_footer = '</svg>'
    with open(output_svg, "w") as f:
        f.write(svg_header + "\n" + "\n".join(svg_elements) + "\n" + svg_footer)
    print(f"Saved vector map to {output_svg}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_deletions.py <input_file_or_dir> [output_image]")
        sys.exit(1)

    target = sys.argv[1]
    output_image = sys.argv[2] if len(sys.argv) > 2 else "deletion_variant_map.svg"

    # 1. Discover files
    files = []
    if os.path.isdir(target):
        files = glob.glob(os.path.join(target, "*_0.01_5_20_2_34.txt"))
        if not files:
            files = glob.glob(os.path.join(target, "*count.txt"))
    elif os.path.isfile(target):
        files = [target]
    else:
        files = glob.glob(target)

    if not files:
        print(f"Error: No matching files found for: {target}")
        sys.exit(1)

    records = []
    for f in files:
        try:
            with open(f, "r") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    try:
                        start = float(row["start"])
                        end = float(row["end"])
                    except (KeyError, ValueError):
                        continue
                    if start >= end:
                        continue

                    freq = None
                    if "Frequency" in row and row["Frequency"]:
                        try:
                            freq = float(row["Frequency"])
                        except ValueError:
                            pass
                    if freq is None and "depth" in row and "MinCov" in row:
                        try:
                            freq = float(row["depth"]) / max(float(row["MinCov"]), 1.0)
                        except ValueError:
                            pass
                    if freq is None:
                        freq = 0.01

                    log_freq = math.log10(max(freq, 1e-4))
                    log_freq = max(-2.0, min(0.0, log_freq))
                    # Linear alpha from -2.0 -> 0.0
                    alpha = 0.05 + ((log_freq - (-2.0)) / 2.0) * 0.95

                    records.append({
                        "start": start,
                        "end": end,
                        "freq": freq,
                        "logFreq": log_freq,
                        "alpha": alpha,
                        "sample": os.path.basename(f)
                    })
        except Exception as e:
            print(f"Warning reading {f}: {e}")

    if not records:
        print("Error: No valid deletion records found.")
        sys.exit(1)

    print(f"Loaded {len(records)} deletion variant junction(s) across {len(files)} file(s).")

    # If matplotlib and pandas are available and user didn't explicitly request .svg
    has_mpl = "matplotlib.pyplot" in sys.modules and not output_image.lower().endswith(".svg")
    if has_mpl:
        try:
            import pandas as pd
            import matplotlib.pyplot as plt
            import numpy as np

            fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
            ax.plot([0, 31500], [0, 31500], linestyle="--", color="gray", alpha=0.6, linewidth=1)
            ends = [r["end"] for r in records]
            starts = [r["start"] for r in records]
            alphas = np.array([r["alpha"] for r in records])

            colors = np.zeros((len(records), 4))
            colors[:, 0] = 0.0
            colors[:, 1] = 0.2
            colors[:, 2] = 0.8
            colors[:, 3] = alphas

            ax.scatter(ends, starts, c=colors, s=15, edgecolors="none")
            ax.set_xlim(0, 31500)
            ax.set_ylim(0, 31500)
            ax.set_aspect("equal")
            ax.set_xlabel("Stop Position (nt)", fontsize=14)
            ax.set_ylabel("Start Position (nt)", fontsize=14)
            ax.set_title("SARS-CoV-2 Deletion Variants", fontsize=13)
            ax.grid(False)
            plt.tight_layout()
            plt.savefig(output_image, dpi=300)
            print(f"Saved deletion variant map to {output_image}")
            return
        except Exception as e:
            print(f"Matplotlib rendering failed ({e}), falling back to vector SVG...")

    # Fallback or SVG direct render
    if not output_image.lower().endswith(".svg"):
        output_image = os.path.splitext(output_image)[0] + ".svg"
    generate_svg(records, output_image)

if __name__ == "__main__":
    main()

