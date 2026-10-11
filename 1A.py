# -*- coding: utf-8 -*-
"""
Kaplan-Meier plot of TCGA patients with TP53 mutations at codon 175 or 273.

Mutation calls come from the MC3 v0.2.8 PUBLIC MAF file of the TCGA
Pan-Cancer Atlas. A patient is assigned to the R175 group when the protein
change starts at residue 175 (p.R175...) and to the R273 group when it
starts at residue 273 (p.R273...). Patients carrying both changes are
assigned by their first row in the MAF. Overall survival comes from the
TCGA-CDR SupplementalTableS1 (OS event and OS.time in days).

Inputs (downloaded from https://api.gdc.cancer.gov/data/...):
  Data/mc3.v0.2.8.PUBLIC.maf.gz
  Data/TCGA-CDR-SupplementalTableS1.xlsx

Outputs (written next to this script):
  Outputs/Fig1A_KM_R175_vs_R273.png
  Outputs/Fig1A_KM_R175_vs_R273.pdf
  Outputs/survival_table.csv
  Outputs/stats.txt
"""

import gzip
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

HERE = os.path.dirname(os.path.abspath(__file__))
MAF = os.path.join(HERE, "Data", "mc3.v0.2.8.PUBLIC.maf.gz")
CDR = os.path.join(HERE, "Data", "TCGA-CDR-SupplementalTableS1.xlsx")
OUTDIR = os.path.join(HERE, "Outputs")
os.makedirs(OUTDIR, exist_ok=True)

COLORS = {"R175X": "#E15759", "R273X": "#F28E2B"}
LABELS = {"R175X": "R175", "R273X": "R273"}


def read_mutations(path):
    keep = []
    header = None
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            if not (line.startswith("TP53\t") or "\tTP53\t" in line):
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row.get("Hugo_Symbol") != "TP53":
                continue
            aa = row.get("HGVSp_Short") or ""
            if aa.startswith("p.R175"):
                group = "R175X"
            elif aa.startswith("p.R273"):
                group = "R273X"
            else:
                continue
            keep.append({
                "patient": row["Tumor_Sample_Barcode"][:12],
                "group": group,
                "HGVSp_Short": aa,
            })
    mut = pd.DataFrame(keep)
    mut = mut.drop_duplicates(subset=["patient"], keep="first")
    return mut


def read_clinical(path):
    cdr = pd.read_excel(path)
    cdr = cdr[["bcr_patient_barcode", "OS", "OS.time"]].copy()
    cdr.columns = ["patient", "OS_event", "OS_time_days"]
    cdr["patient"] = cdr["patient"].astype(str).str.slice(0, 12)
    cdr["OS_time_days"] = pd.to_numeric(cdr["OS_time_days"], errors="coerce")
    cdr = cdr.dropna(subset=["OS_time_days"])
    return cdr


def draw_km(ax, df):
    fits = {}
    for group in ("R175X", "R273X"):
        d = df[df["group"] == group]
        kmf = KaplanMeierFitter(label=LABELS[group])
        kmf.fit(d["OS_time_days"], event_observed=d["OS_event"])
        fits[group] = kmf
        kmf.plot_survival_function(
            ax=ax,
            ci_show=False,
            show_censors=True,
            censor_styles={"marker": "|", "ms": 5, "mew": 1.0},
            color=COLORS[group],
            linewidth=1.4,
            legend=False,
        )
    g1 = df[df["group"] == "R175X"]
    g2 = df[df["group"] == "R273X"]
    lr = logrank_test(g1["OS_time_days"], g2["OS_time_days"],
                      g1["OS_event"], g2["OS_event"])
    p = lr.p_value
    p_text = "p < 0.0001" if p < 0.0001 else ("p = %.4f" % p if p < 0.001
                                              else "p = %.2f" % p)
    ax.text(0.05, 0.35, p_text, transform=ax.transAxes, fontsize=9)

    handles = [Line2D([0], [0], color=COLORS[g], lw=1.6) for g in ("R175X", "R273X")]
    labels = ["%s (n=%d)" % (LABELS[g], (df["group"] == g).sum())
              for g in ("R175X", "R273X")]
    ax.legend(handles, labels, title="Strata", loc="upper right",
              frameon=False, fontsize=8, title_fontsize=8, handlelength=1.2)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0.0, 0.25, 0.50, 0.75, 1.00])
    ax.set_yticklabels(["0.00", "0.25", "0.50", "0.75", "1.00"])
    ax.set_xlim(0, df["OS_time_days"].max() * 1.03)
    ax.set_xticks([0, 2000, 4000, 6000])
    ax.tick_params(labelsize=8)
    ax.set_xlabel("Time (days)", fontsize=9)
    ax.set_ylabel("Overall survival probability", fontsize=9)
    return fits, p


def main():
    mut = read_mutations(MAF)
    cdr = read_clinical(CDR)
    df = mut.merge(cdr, on="patient", how="inner")
    df = df.sort_values(["group", "patient"]).reset_index(drop=True)
    df.to_csv(os.path.join(OUTDIR, "survival_table.csv"), index=False)

    fig, ax = plt.subplots(figsize=(4.6, 4.4), dpi=300)
    fits, p = draw_km(ax, df)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "Fig1A_KM_R175_vs_R273.png"))
    fig.savefig(os.path.join(OUTDIR, "Fig1A_KM_R175_vs_R273.pdf"))
    plt.close(fig)

    lines = []
    lines.append("n R175X = %d" % (df["group"] == "R175X").sum())
    lines.append("n R273X = %d" % (df["group"] == "R273X").sum())
    lines.append("log-rank p = %.4g" % p)
    for group in ("R175X", "R273X"):
        d = df[df["group"] == group]
        med = fits[group].median_survival_time_
        med_txt = "not reached" if np.isinf(med) else "%.0f days" % med
        lines.append("%s: n=%d, events=%d, median OS = %s"
                     % (group, len(d), int(d["OS_event"].sum()), med_txt))
    with open(os.path.join(OUTDIR, "stats.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
