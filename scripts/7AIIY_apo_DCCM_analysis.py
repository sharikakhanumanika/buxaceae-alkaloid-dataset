import numpy as np
import pandas as pd

# ============================================================
# Load files
# ============================================================

cor_mat = np.loadtxt("dccm_matrix.csv", delimiter=",")
res_ids = np.loadtxt("residue_ids.txt", dtype=int)

print(f"Loaded DCCM matrix: {cor_mat.shape}")
print(f"Loaded residue list: {len(res_ids)} residues")

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
# BChE (7AIY)
# ============================================================

# ============================================================
# BChE (7AIY)
# Residues selected from docking interactions
# ============================================================

active_site = [
    70,   # Asp70
    79,   # Ser79
    82,   # Trp82
    116,  # Gly116
    117,  # Gly117
    119,  # Gln119
    120,  # Thr120
    197,  # Glu197
    198,  # Ser198
    199,  # Ala199
    231,  # Trp231
    277,  # Ala277
    285,  # Pro285
    286,  # Leu286
    287,  # Ser287
    288,  # Val288
    328,  # Ala328
    329,  # Phe329
    332,  # Tyr332
    398,  # Phe398
    430,  # Trp430
    437,  # Met437
    438,  # His438
    439,  # Gly439
    440   # Tyr440
]

# Representative residues for different structural regions
domain_reps = {
    "N_terminal":  [5, 10, 20, 30],
    "Omega_loop":  [68, 72, 78, 82],
    "Gorge_entry": [280, 285, 290, 295],
    "Gorge_core":  [325, 330, 340, 350],
    "Acyl_pocket": [438, 440, 445, 450],
    "C_terminal":  [510, 515, 520, 525]
}

# ============================================================
# ACTIVE SITE CORRELATIONS
# ============================================================

print("\nCalculating active-site correlations...")

rows=[]

for i,r1 in enumerate(active_site):
    for r2 in active_site[i+1:]:

        c=corr(r1,r2)

        if c is not None:

            rows.append({
                "Residue_1":r1,
                "Residue_2":r2,
                "Correlation":round(c,6)
            })

active_df=pd.DataFrame(rows)

active_df.to_csv(
    "active_site_correlations.csv",
    index=False
)

print(f"Saved active_site_correlations.csv ({len(active_df)} pairs)")

# ============================================================
# DOMAIN CORRELATIONS
# ============================================================

print("Calculating domain correlations...")

domain_rows=[]

names=list(domain_reps.keys())

for i,d1 in enumerate(names):

    for d2 in names[i+1:]:

        vals=[]

        for r1 in domain_reps[d1]:
            for r2 in domain_reps[d2]:

                c=corr(r1,r2)

                if c is not None:
                    vals.append(c)

        vals=np.array(vals)

        domain_rows.append({
            "Domain_1":d1,
            "Domain_2":d2,
            "Mean":round(vals.mean(),6),
            "SD":round(vals.std(),6),
            "Minimum":round(vals.min(),6),
            "Maximum":round(vals.max(),6),
            "N_pairs":len(vals)
        })

domain_df=pd.DataFrame(domain_rows)

domain_df.to_csv(
    "domain_correlations.csv",
    index=False
)

print(f"Saved domain_correlations.csv")

# ============================================================
# OVERALL STATISTICS
# ============================================================

print("Calculating overall statistics...")

mask=~np.eye(len(cor_mat),dtype=bool)

off=cor_mat[mask]

stats={

    "Residues":len(res_ids),

    "Total_pairs":len(off),

    "Mean":off.mean(),

    "Standard_deviation":off.std(),

    "Minimum":off.min(),

    "Maximum":off.max(),

    "Strong_positive_(>0.50)":np.sum(off>0.50),

    "Moderate_positive_(>0.25)":np.sum(off>0.25),

    "Moderate_negative_(<-0.25)":np.sum(off<-0.25),

    "Strong_negative_(<-0.50)":np.sum(off<-0.50)

}

with open("dccm_statistics.txt","w") as f:

    f.write("DCCM SUMMARY STATISTICS\n")
    f.write("="*40+"\n\n")

    for k,v in stats.items():

        if isinstance(v,float):
            f.write(f"{k:30s}: {v:.6f}\n")
        else:
            f.write(f"{k:30s}: {v}\n")

print("Saved dccm_statistics.txt")

# ============================================================
# TERMINAL SUMMARY
# ============================================================

print("\n========================================")
print("DCCM ANALYSIS COMPLETED")
print("========================================")

print(f"Residues           : {len(res_ids)}")
print(f"Active-site pairs  : {len(active_df)}")
print(f"Domain comparisons : {len(domain_df)}")

print("\nGenerated files:")
print("  active_site_correlations.csv")
print("  domain_correlations.csv")
print("  dccm_statistics.txt")

print("\nDone!")