import numpy as np
import pandas as pd

# ============================================================
# Load DCCM matrix and residue IDs
# ============================================================

cor_mat = np.loadtxt("dccm_matrix.csv", delimiter=",")
res_ids = np.loadtxt("residue_ids.txt", dtype=int)

print(f"Loaded DCCM matrix: {cor_mat.shape}")
print(f"Loaded residue list: {len(res_ids)} residues")

# ============================================================
# AChE (7E3H)
# Active-site residues obtained from docking interactions
# ============================================================

active_site = [
    72, 74, 75, 76,
    96,
    124,
    235,
    243,
    247,
    287,
    289,
    291,
    292,
    293,
    294,
    295,
    296,
    297,
    313,
    337,
    338,
    341,
    342,
    368,
    369,
    370,
    405,
    409,
    410,
    413,
    447,
    532,
    533,
    536,
    537,
    540
]

# Representative residues for major structural regions
domain_reps = {
    "N_terminal": [4, 8, 12, 16],
    "PAS_loop": [72, 76, 80, 84],
    "Mid_gorge": [199, 203, 231, 246],
    "Gorge_core": [286, 295, 337, 341],
    "Catalytic": [203, 334, 447],
    "C_terminal": [510, 520, 530, 540]
}

# ============================================================
# Check residues
# ============================================================

print("\nChecking residue availability...")

missing = [r for r in active_site if r not in res_ids]

if len(missing) == 0:
    print("All active-site residues found.")
else:
    print("Missing residues:", missing)

for domain, residues in domain_reps.items():
    miss = [r for r in residues if r not in res_ids]
    print(f"{domain}: Missing {miss}")

# ============================================================
# Helper functions
# ============================================================

def get_idx(resid):
    idx = np.where(res_ids == resid)[0]
    if len(idx) == 0:
        return None
    return idx[0]

def corr(r1, r2):
    i = get_idx(r1)
    j = get_idx(r2)

    if i is None or j is None:
        return None

    return float(cor_mat[i, j])

# ============================================================
# Active-site correlations
# ============================================================

print("\nCalculating active-site correlations...")

rows = []

for i, r1 in enumerate(active_site):

    for r2 in active_site[i + 1:]:

        c = corr(r1, r2)

        if c is not None:

            rows.append({
                "Residue_1": r1,
                "Residue_2": r2,
                "Correlation": round(c, 6)
            })

active_df = pd.DataFrame(rows)

active_df.to_csv(
    "active_site_correlations.csv",
    index=False
)

print(f"Saved active_site_correlations.csv ({len(active_df)} pairs)")

# ============================================================
# Domain correlations
# ============================================================

print("\nCalculating domain correlations...")

domain_rows = []

names = list(domain_reps.keys())

for i, d1 in enumerate(names):

    for d2 in names[i + 1:]:

        vals = []

        for r1 in domain_reps[d1]:

            for r2 in domain_reps[d2]:

                c = corr(r1, r2)

                if c is not None:
                    vals.append(c)

        vals = np.array(vals)

        if len(vals) > 0:

            domain_rows.append({
                "Domain_1": d1,
                "Domain_2": d2,
                "Mean": round(vals.mean(), 6),
                "SD": round(vals.std(), 6),
                "Minimum": round(vals.min(), 6),
                "Maximum": round(vals.max(), 6),
                "N_pairs": len(vals)
            })

domain_df = pd.DataFrame(domain_rows)

domain_df.to_csv(
    "domain_correlations.csv",
    index=False
)

print("Saved domain_correlations.csv")

# ============================================================
# Overall statistics
# ============================================================

print("\nCalculating overall statistics...")

mask = ~np.eye(len(cor_mat), dtype=bool)

off = cor_mat[mask]

stats = {
    "Residues": len(res_ids),
    "Total_pairs": len(off),
    "Mean": off.mean(),
    "Standard_deviation": off.std(),
    "Minimum": off.min(),
    "Maximum": off.max(),
    "Strong_positive_(>0.50)": np.sum(off > 0.50),
    "Moderate_positive_(>0.25)": np.sum(off > 0.25),
    "Moderate_negative_(<-0.25)": np.sum(off < -0.25),
    "Strong_negative_(<-0.50)": np.sum(off < -0.50)
}

with open("dccm_statistics.txt", "w") as f:

    f.write("DCCM SUMMARY STATISTICS\n")
    f.write("=" * 45 + "\n\n")

    for k, v in stats.items():

        if isinstance(v, float):
            f.write(f"{k:30s}: {v:.6f}\n")
        else:
            f.write(f"{k:30s}: {v}\n")

print("Saved dccm_statistics.txt")

# ============================================================
# Terminal summary
# ============================================================

print("\n==========================================")
print("DCCM ANALYSIS COMPLETED")
print("==========================================")

print(f"Residues           : {len(res_ids)}")
print(f"Active-site pairs  : {len(active_df)}")
print(f"Domain comparisons : {len(domain_df)}")

print("\nGenerated files:")
print("  active_site_correlations.csv")
print("  domain_correlations.csv")
print("  dccm_statistics.txt")

print("\nDone!")