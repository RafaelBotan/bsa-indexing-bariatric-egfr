#!/usr/bin/env python3
"""
bootstrap_shapley.py — Bootstrap CIs for Shapley decomposition components.
Bootstrap by patient (preserving PRE-12M pair), stratified by cohort in pooled.
Output: tabelas_v4/Tab17_shapley_bootstrap.csv
"""
import pandas as pd
import numpy as np
from scipy import stats
import os, warnings
warnings.filterwarnings('ignore')

OUT = 'Y:/Base Rafael do Sergio/tabelas_v4'
os.makedirs(OUT, exist_ok=True)

N_BOOT = 5000
SEED = 42
np.random.seed(SEED)

# Load dataset
df = pd.read_csv('Y:/Base Rafael do Sergio/ESTUDO_C_DATASET_v4.csv', low_memory=False)

# Filter to 12M paired with NI available
paired = df[
    df['EGFR_CKD_EPI_PRE'].notna() &
    df['EGFR_CKD_EPI_12M'].notna() &
    df['EGFR_NI_PRE'].notna() &
    df['EGFR_NI_12M'].notna() &
    df['BSA_DUBOIS_PRE'].notna() &
    df['BSA_DUBOIS_12M'].notna()
].copy()

print(f"Paired PRE+12M with NI: {len(paired)} ({(paired['COORTE']=='Arruda').sum()} Arruda, {(paired['COORTE']=='Galvao').sum()} Galvao)")

# Compute Shapley components per patient
paired['delta_egfr_idx'] = paired['EGFR_CKD_EPI_12M'] - paired['EGFR_CKD_EPI_PRE']
paired['delta_bsa'] = paired['BSA_DUBOIS_12M'] - paired['BSA_DUBOIS_PRE']
paired['delta_ni'] = paired['EGFR_NI_12M'] - paired['EGFR_NI_PRE']
paired['egfr_mean'] = (paired['EGFR_CKD_EPI_PRE'] + paired['EGFR_CKD_EPI_12M']) / 2
paired['bsa_mean'] = (paired['BSA_DUBOIS_PRE'] + paired['BSA_DUBOIS_12M']) / 2

# Shapley decomposition
paired['comp_egfr'] = paired['delta_egfr_idx'] * paired['bsa_mean'] / 1.73
paired['comp_bsa'] = paired['egfr_mean'] * paired['delta_bsa'] / 1.73
paired['pct_bsa'] = np.abs(paired['comp_bsa']) / (np.abs(paired['comp_egfr']) + np.abs(paired['comp_bsa'])) * 100

# Verify decomposition
check = (paired['comp_egfr'] + paired['comp_bsa'] - paired['delta_ni']).abs()
print(f"Decomposition check: max residual = {check.max():.6f}")


def bootstrap_shapley(data, n_boot=N_BOOT):
    """Bootstrap by patient, return CIs for each Shapley component."""
    n = len(data)
    boot_delta_ni = np.zeros(n_boot)
    boot_comp_egfr = np.zeros(n_boot)
    boot_comp_bsa = np.zeros(n_boot)
    boot_pct_bsa = np.zeros(n_boot)

    for b in range(n_boot):
        idx = np.random.randint(0, n, size=n)
        sample = data.iloc[idx]
        boot_delta_ni[b] = sample['delta_ni'].mean()
        boot_comp_egfr[b] = sample['comp_egfr'].mean()
        boot_comp_bsa[b] = sample['comp_bsa'].mean()
        # %BSA from means, not mean of %
        total_abs = np.abs(sample['comp_egfr'].mean()) + np.abs(sample['comp_bsa'].mean())
        if total_abs > 0:
            boot_pct_bsa[b] = np.abs(sample['comp_bsa'].mean()) / total_abs * 100
        else:
            boot_pct_bsa[b] = np.nan

    return {
        'delta_ni': boot_delta_ni,
        'comp_egfr': boot_comp_egfr,
        'comp_bsa': boot_comp_bsa,
        'pct_bsa': boot_pct_bsa
    }


def ci_percentile(arr, alpha=0.05):
    """Percentile CI."""
    lo = np.nanpercentile(arr, 100 * alpha / 2)
    hi = np.nanpercentile(arr, 100 * (1 - alpha / 2))
    return lo, hi


results = []

# By cohort
for label in ['Arruda', 'Galvao']:
    sub = paired[paired['COORTE'] == label]
    n = len(sub)
    if n < 10:
        print(f"\n{label}: N={n} -- too few, skipping")
        continue

    print(f"\n{'='*60}")
    print(f"{label} (N={n}) — {N_BOOT} bootstrap replications")
    print(f"{'='*60}")

    boot = bootstrap_shapley(sub)

    for metric, col in [('Delta_NI', 'delta_ni'), ('Comp_eGFR', 'comp_egfr'),
                          ('Comp_BSA', 'comp_bsa'), ('Pct_BSA', 'pct_bsa')]:
        point = sub[col].mean() if col != 'pct_bsa' else (
            np.abs(sub['comp_bsa'].mean()) / (np.abs(sub['comp_egfr'].mean()) + np.abs(sub['comp_bsa'].mean())) * 100
        )
        lo, hi = ci_percentile(boot[col])
        print(f"  {metric}: {point:.2f} [{lo:.2f}, {hi:.2f}]")
        results.append({
            'Cohort': label, 'N': n, 'Metric': metric,
            'Point_estimate': round(point, 2),
            'CI95_lo': round(lo, 2), 'CI95_hi': round(hi, 2),
            'Boot_SE': round(np.nanstd(boot[col]), 2)
        })

# Pooled — stratified bootstrap (resample within each cohort, then combine)
print(f"\n{'='*60}")
print(f"Pooled stratified (N={len(paired)}) — {N_BOOT} bootstrap replications")
print(f"{'='*60}")

arr_data = paired[paired['COORTE'] == 'Arruda']
gal_data = paired[paired['COORTE'] == 'Galvao']
n_arr = len(arr_data)
n_gal = len(gal_data)
n_total = n_arr + n_gal

boot_pooled = {k: np.zeros(N_BOOT) for k in ['delta_ni', 'comp_egfr', 'comp_bsa', 'pct_bsa']}

for b in range(N_BOOT):
    # Resample within each cohort
    idx_a = np.random.randint(0, n_arr, size=n_arr)
    idx_g = np.random.randint(0, n_gal, size=n_gal)
    sample = pd.concat([arr_data.iloc[idx_a], gal_data.iloc[idx_g]])

    boot_pooled['delta_ni'][b] = sample['delta_ni'].mean()
    boot_pooled['comp_egfr'][b] = sample['comp_egfr'].mean()
    boot_pooled['comp_bsa'][b] = sample['comp_bsa'].mean()
    total_abs = np.abs(sample['comp_egfr'].mean()) + np.abs(sample['comp_bsa'].mean())
    boot_pooled['pct_bsa'][b] = np.abs(sample['comp_bsa'].mean()) / total_abs * 100 if total_abs > 0 else np.nan

for metric, col in [('Delta_NI', 'delta_ni'), ('Comp_eGFR', 'comp_egfr'),
                      ('Comp_BSA', 'comp_bsa'), ('Pct_BSA', 'pct_bsa')]:
    point = paired[col].mean() if col != 'pct_bsa' else (
        np.abs(paired['comp_bsa'].mean()) / (np.abs(paired['comp_egfr'].mean()) + np.abs(paired['comp_bsa'].mean())) * 100
    )
    lo, hi = ci_percentile(boot_pooled[col])
    print(f"  {metric}: {point:.2f} [{lo:.2f}, {hi:.2f}]")
    results.append({
        'Cohort': 'Pooled', 'N': n_total, 'Metric': metric,
        'Point_estimate': round(point, 2),
        'CI95_lo': round(lo, 2), 'CI95_hi': round(hi, 2),
        'Boot_SE': round(np.nanstd(boot_pooled[col]), 2)
    })

# Save
res_df = pd.DataFrame(results)
res_df.to_csv(f'{OUT}/Tab17_shapley_bootstrap.csv', index=False)
print(f"\nSaved to {OUT}/Tab17_shapley_bootstrap.csv")
print("\n" + res_df.to_string(index=False))
print("\nDone.")
