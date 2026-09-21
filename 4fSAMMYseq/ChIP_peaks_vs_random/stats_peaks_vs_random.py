#!/usr/bin/env python3
import sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

peaks_tsv, rand_tsv, mark, out_txt, out_pdf = sys.argv[1:]

# -------------------------------------------------
# ALWAYS read as no header and force column names
# -------------------------------------------------
def load_table(path):
    df = pd.read_csv(path, sep="\t", header=None, comment="#")

    # keep only first 6 columns if extra garbage exists
    df = df.iloc[:, :6]

    df.columns = ["chr","start","end","WT1","WT2","WT3"]
    return df

df_peaks = load_table(peaks_tsv)
df_rand  = load_table(rand_tsv)

# convert numeric safely
for c in ["WT1","WT2","WT3"]:
    df_peaks[c] = pd.to_numeric(df_peaks[c], errors="coerce")
    df_rand[c]  = pd.to_numeric(df_rand[c], errors="coerce")

# mean across replicates
x = df_peaks[["WT1","WT2","WT3"]].mean(axis=1).dropna().to_numpy()
y = df_rand [["WT1","WT2","WT3"]].mean(axis=1).dropna().to_numpy()

if len(x)==0 or len(y)==0:
    print("ERROR: empty vectors after filtering")
    sys.exit(1)

# Mann–Whitney
u,p = mannwhitneyu(x,y,alternative="two-sided")

median_peaks=float(np.median(x))
median_rand=float(np.median(y))
mean_peaks=float(np.mean(x))
mean_rand=float(np.mean(y))

# Cliff's delta
def cliffs_delta(a,b):
    gt=sum(np.sum(ai>b) for ai in a)
    lt=sum(np.sum(ai<b) for ai in a)
    return (gt-lt)/(len(a)*len(b))

delta=cliffs_delta(x,y)

# save stats
with open(out_txt,"w") as f:
    f.write(f"MARK\t{mark}\n")
    f.write(f"N_peaks\t{len(x)}\n")
    f.write(f"N_random\t{len(y)}\n")
    f.write(f"median_peaks\t{median_peaks:.6f}\n")
    f.write(f"median_random\t{median_rand:.6f}\n")
    f.write(f"mean_peaks\t{mean_peaks:.6f}\n")
    f.write(f"mean_random\t{mean_rand:.6f}\n")
    f.write(f"MannWhitney_U\t{u}\n")
    f.write(f"p_value\t{p:.6e}\n")
    f.write(f"cliffs_delta\t{delta:.6f}\n")

print(f"[STATS DONE] {mark}  p={p:.3e}")

# plot
plt.figure(figsize=(5,6))
plt.boxplot([x,y],labels=["Peaks","Random"],showfliers=False)
plt.ylabel("SAMMY signal")
plt.title(f"{mark}\np={p:.2e}")
plt.tight_layout()
plt.savefig(out_pdf)