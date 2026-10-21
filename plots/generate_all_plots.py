#!/usr/bin/env python3
import os
import csv
import numpy as np

# Ensure matplotlib uses a writable temp cache directory
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib_cache"
import matplotlib.pyplot as plt

# Output directories
out_dirs = ["/home/sun/yin/plots", "/home/sun/yin/SARS2_DVG/plots"]
for d in out_dirs:
    os.makedirs(d, exist_ok=True)

# Load data
csv_file = "/home/sun/yin/SARS2_DVG/transferred_outs/benchmark_100_samples.csv"
rows = []
with open(csv_file, "r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append(r)

# Filter successful runs for performance benchmarks
succ = [r for r in rows if "SUCCESS" in r["status"]]
print(f"Loaded {len(rows)} samples ({len(succ)} successful).")

# Extract numeric arrays
s1 = np.array([float(r["step1_fasterq_sec"]) for r in succ])
s2 = np.array([float(r["step2_trimgalore_sec"]) for r in succ])
s3_sec = np.array([float(r["step3_star_alignment_sec"]) for r in succ])
s3_min = s3_sec / 60.0
s4 = np.array([float(r["step4_standardize_filter_sec"]) for r in succ])
s5 = np.array([float(r["step5_sort_depth_r_sec"]) for r in succ])
s6 = np.array([float(r["step6_cleanup_sec"]) for r in succ])

tot_min = np.array([float(r["total_pipeline_time_min"]) for r in succ])
spots = np.array([float(r["total_spots"]) for r in succ])
mbp = np.array([float(r["total_sequenced_mbp"]) for r in succ])
sra_mb = np.array([float(r["sra_size_mb"]) for r in succ])
dels = np.array([float(r["deletions_detected"]) for r in succ])
depth = np.array([float(r["mean_coverage_depth"]) for r in succ])

# Disk arrays
# Peak scratch in GB: htcondor_disk_kb / 1024 / 1024
peak_disk_gb = np.array([float(r["htcondor_disk_kb"]) / (1024.0 * 1024.0) for r in succ])

# Final scratch in MB: parse "87M" -> 87.0
final_disk_mb = []
for r in succ:
    val_str = r["final_scratch_disk"].replace("M", "").replace("G", "")
    try:
        final_disk_mb.append(float(val_str))
    except:
        final_disk_mb.append(87.0)
final_disk_mb = np.array(final_disk_mb)

# RAM arrays
cgroup_ram_mb = np.array([float(r["cgroup_peak_ram_mb"]) for r in succ])
condor_rss_mb = np.array([float(r["htcondor_max_rss_mb"]) for r in succ])

# Clean styling helper
def save_plot(fig, filename):
    for d in out_dirs:
        path = os.path.join(d, filename)
        fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved: {filename}")
    plt.close(fig)

# -------------------------------------------------------------------------
# Plot 1: Average Time Spent on Each Step
# -------------------------------------------------------------------------
step_names = [
    "Step 1:\nFASTQ Extr.",
    "Step 2:\nTrim Galore",
    "Step 3:\nSTAR Align",
    "Step 4:\nDVG Filter",
    "Step 5:\nSort & Depth",
    "Step 6:\nCleanup"
]
step_means_sec = [np.mean(s1), np.mean(s2), np.mean(s3_sec), np.mean(s4), np.mean(s5), np.mean(s6)]
step_means_min = [s / 60.0 for s in step_means_sec]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Subplot 1: All steps in minutes
bars1 = ax1.bar(step_names, step_means_min, color="#2563eb", width=0.55)
ax1.set_ylabel("Average Time (Minutes)")
ax1.set_title("Average Time per Step (All Steps)")
ax1.grid(axis="y", linestyle="--", alpha=0.3)
for bar in bars1:
    y = bar.get_height()
    label = f"{y:.1f} m" if y >= 1.0 else f"{y*60:.0f} s"
    ax1.annotate(label, xy=(bar.get_x() + bar.get_width() / 2, y),
                 xytext=(0, 4), textcoords="offset points",
                 ha="center", va="bottom", fontsize=9, fontweight="bold")

# Subplot 2: Non-STAR steps in seconds (detailed zoom)
non_star_names = ["Step 1\n(fasterq)", "Step 2\n(TrimGalore)", "Step 4\n(DVG Filter)", "Step 5\n(Sort/Depth)", "Step 6\n(Cleanup)"]
non_star_means = [np.mean(s1), np.mean(s2), np.mean(s4), np.mean(s5), np.mean(s6)]
bars2 = ax2.bar(non_star_names, non_star_means, color="#059669", width=0.55)
ax2.set_ylabel("Average Time (Seconds)")
ax2.set_title("Pre/Post-Processing Steps (Excluding STAR)")
ax2.grid(axis="y", linestyle="--", alpha=0.3)
for bar in bars2:
    y = bar.get_height()
    ax2.annotate(f"{y:.1f} s", xy=(bar.get_x() + bar.get_width() / 2, y),
                 xytext=(0, 4), textcoords="offset points",
                 ha="center", va="bottom", fontsize=9, fontweight="bold")

plt.tight_layout()
save_plot(fig, "01_average_time_per_step.png")

# -------------------------------------------------------------------------
# Plot 2: Distribution of Total Pipeline Time per Sample
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(tot_min, bins=25, color="#2563eb", edgecolor="white", alpha=0.85)
med_time = np.median(tot_min)
ax.axvline(med_time, color="#dc2626", linestyle="--", linewidth=1.8, label=f"Median: {med_time:.1f} min")
ax.set_xlabel("Total Pipeline Time (Minutes)")
ax.set_ylabel("Number of Samples")
ax.set_title("Distribution of Total Pipeline Time (n=99)")
ax.legend(frameon=True)
ax.grid(axis="y", linestyle="--", alpha=0.3)
plt.tight_layout()
save_plot(fig, "02_pipeline_total_time_distribution.png")

# -------------------------------------------------------------------------
# Plot 3: Distribution of Time for Each Individual Pipeline Step
# -------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(12, 7))

# Step 1
axes[0, 0].hist(s1, bins=15, color="#0284c7", edgecolor="white", alpha=0.85)
axes[0, 0].set_title("Step 1: FASTQ Extraction (fasterq-dump)")
axes[0, 0].set_xlabel("Time (Seconds)")
axes[0, 0].set_ylabel("Sample Count")
axes[0, 0].grid(axis="y", linestyle="--", alpha=0.3)

# Step 2
axes[0, 1].hist(s2, bins=15, color="#0284c7", edgecolor="white", alpha=0.85)
axes[0, 1].set_title("Step 2: Quality Trimming (Trim Galore)")
axes[0, 1].set_xlabel("Time (Seconds)")
axes[0, 1].set_ylabel("Sample Count")
axes[0, 1].grid(axis="y", linestyle="--", alpha=0.3)

# Step 3
axes[0, 2].hist(s3_min, bins=20, color="#2563eb", edgecolor="white", alpha=0.85)
axes[0, 2].set_title("Step 3: STAR Chimeric Alignment")
axes[0, 2].set_xlabel("Time (Minutes)")
axes[0, 2].set_ylabel("Sample Count")
axes[0, 2].grid(axis="y", linestyle="--", alpha=0.3)

# Step 4
axes[1, 0].hist(s4, bins=15, color="#0284c7", edgecolor="white", alpha=0.85)
axes[1, 0].set_title("Step 4: DVG Filtering & Extraction")
axes[1, 0].set_xlabel("Time (Seconds)")
axes[1, 0].set_ylabel("Sample Count")
axes[1, 0].grid(axis="y", linestyle="--", alpha=0.3)

# Step 5
axes[1, 1].hist(s5, bins=15, color="#0284c7", edgecolor="white", alpha=0.85)
axes[1, 1].set_title("Step 5: Sorting, Depth & R Summary")
axes[1, 1].set_xlabel("Time (Seconds)")
axes[1, 1].set_ylabel("Sample Count")
axes[1, 1].grid(axis="y", linestyle="--", alpha=0.3)

# Step 6
axes[1, 2].hist(s6, bins=10, color="#0284c7", edgecolor="white", alpha=0.85)
axes[1, 2].set_title("Step 6: Intermediate File Cleanup")
axes[1, 2].set_xlabel("Time (Seconds)")
axes[1, 2].set_ylabel("Sample Count")
axes[1, 2].grid(axis="y", linestyle="--", alpha=0.3)

plt.tight_layout()
save_plot(fig, "03_step_time_distributions.png")

# -------------------------------------------------------------------------
# Plot 4: Storage Usage
# -------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Subplot 1: Peak scratch vs Final scratch distribution
ax1.hist(peak_disk_gb, bins=15, color="#0284c7", edgecolor="white", alpha=0.75, label="Peak Scratch (Pre-Cleanup)")
ax1.set_xlabel("Scratch Storage (GB)")
ax1.set_ylabel("Sample Count")
ax1.set_title("Scratch Storage Distribution")
ax1.axvline(10.0, color="#dc2626", linestyle="--", linewidth=1.5, label="Request Limit (10 GB)")
ax1.legend(frameon=True)
ax1.grid(axis="y", linestyle="--", alpha=0.3)

# Subplot 2: SRA Input Size vs Peak Scratch
ax2.scatter(sra_mb, peak_disk_gb, color="#2563eb", edgecolors="white", s=45, alpha=0.85)
# Add linear fit
m, b = np.polyfit(sra_mb, peak_disk_gb, 1)
x_line = np.linspace(min(sra_mb), max(sra_mb), 50)
ax2.plot(x_line, m*x_line + b, color="#dc2626", linestyle="--", linewidth=1.5)
ax2.set_xlabel("Input SRA File Size (MB)")
ax2.set_ylabel("Peak Scratch Storage (GB)")
ax2.set_title("SRA File Size vs. Peak Scratch Storage")
ax2.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
save_plot(fig, "04_storage_usage.png")

# -------------------------------------------------------------------------
# Plot 5: RAM Usage
# -------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Subplot 1: Distribution of CGroup Peak RAM vs HTCondor RSS
ax1.hist(cgroup_ram_mb, bins=18, color="#d97706", edgecolor="white", alpha=0.75, label="Container Peak (cgroups)")
ax1.hist(condor_rss_mb, bins=15, color="#059669", edgecolor="white", alpha=0.6, label="Process Tree RSS (HTCondor)")
ax1.axvline(6144, color="#dc2626", linestyle="--", linewidth=1.5, label="Request Limit (6,144 MB)")
ax1.set_xlabel("Memory Usage (MB)")
ax1.set_ylabel("Sample Count")
ax1.set_title("RAM Usage Distribution")
ax1.legend(frameon=True)
ax1.grid(axis="y", linestyle="--", alpha=0.3)

# Subplot 2: Total Sequenced Mbp vs Peak Memory
ax2.scatter(mbp, cgroup_ram_mb, color="#d97706", edgecolors="white", s=45, alpha=0.85)
ax2.axhline(6144, color="#dc2626", linestyle="--", linewidth=1.5, label="Request Limit (6,144 MB)")
ax2.set_xlabel("Total Sequenced Length (Mbp)")
ax2.set_ylabel("Container Peak RAM (MB)")
ax2.set_title("Sequenced Volume vs. Peak Memory")
ax2.legend(frameon=True)
ax2.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
save_plot(fig, "05_ram_usage.png")

# -------------------------------------------------------------------------
# Plot 6: Complexity vs Runtime - Sequencing Volume
# -------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Subplot 1: Spots vs Pipeline Time
ax1.scatter(spots / 1e6, tot_min, color="#2563eb", edgecolors="white", s=45, alpha=0.8)
m1, b1 = np.polyfit(spots / 1e6, tot_min, 1)
x_fit1 = np.linspace(min(spots / 1e6), max(spots / 1e6), 50)
ax1.plot(x_fit1, m1 * x_fit1 + b1, color="#dc2626", linestyle="--", linewidth=1.5)
ax1.set_xlabel("Sample Read Volume (Million Spots)")
ax1.set_ylabel("Total Pipeline Time (Minutes)")
ax1.set_title("Read Volume vs. Total Pipeline Time")
ax1.grid(True, linestyle="--", alpha=0.3)

# Subplot 2: Total Sequenced Mbp vs STAR Alignment Time
ax2.scatter(mbp, s3_min, color="#059669", edgecolors="white", s=45, alpha=0.8)
m2, b2 = np.polyfit(mbp, s3_min, 1)
x_fit2 = np.linspace(min(mbp), max(mbp), 50)
ax2.plot(x_fit2, m2 * x_fit2 + b2, color="#dc2626", linestyle="--", linewidth=1.5)
ax2.set_xlabel("Total Sequenced Length (Mbp)")
ax2.set_ylabel("STAR Alignment Time (Minutes)")
ax2.set_title("Sequenced Length vs. STAR Alignment Time")
ax2.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
save_plot(fig, "06_complexity_vs_runtime_volume.png")

# -------------------------------------------------------------------------
# Plot 7: Complexity vs Runtime - Recombination / DVG Complexity
# -------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

# Subplot 1: Deletions detected vs STAR Alignment Time
ax1.scatter(dels, s3_min, color="#7c3aed", edgecolors="white", s=45, alpha=0.85)
m_del, b_del = np.polyfit(dels, s3_min, 1)
x_del = np.linspace(min(dels), max(dels), 50)
ax1.plot(x_del, m_del * x_del + b_del, color="#dc2626", linestyle="--", linewidth=1.5)
ax1.set_xlabel("Deletions Detected (DVG Junctions)")
ax1.set_ylabel("STAR Alignment Time (Minutes)")
ax1.set_title("DVG Complexity vs. STAR Alignment Time")
ax1.grid(True, linestyle="--", alpha=0.3)

# Subplot 2: Deletions detected vs Total Pipeline Time
ax2.scatter(dels, tot_min, color="#7c3aed", edgecolors="white", s=45, alpha=0.85)
m_tot, b_tot = np.polyfit(dels, tot_min, 1)
ax2.plot(x_del, m_tot * x_del + b_tot, color="#dc2626", linestyle="--", linewidth=1.5)
ax2.set_xlabel("Deletions Detected (DVG Junctions)")
ax2.set_ylabel("Total Pipeline Time (Minutes)")
ax2.set_title("DVG Complexity vs. Total Pipeline Time")
ax2.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
save_plot(fig, "07_complexity_vs_runtime_deletions.png")

# -------------------------------------------------------------------------
# Plot 8: Multi-dimensional Complexity vs Runtime (Volume + DVGs)
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 5))
sc = ax.scatter(mbp, tot_min, c=dels, cmap="viridis", s=50, edgecolors="white", alpha=0.9)
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label("Deletions Detected (DVGs)")
ax.set_xlabel("Total Sequenced Length (Mbp)")
ax.set_ylabel("Total Pipeline Time (Minutes)")
ax.set_title("Runtime by Sequencing Volume and DVG Complexity")
ax.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
save_plot(fig, "08_multidimensional_complexity_runtime.png")

print("All plots generated successfully.")
