#!/usr/bin/env python3
"""
ESTUDO C v4 — Complete Analysis Script
Study of eGFR (indexed vs nonindexed) after bariatric surgery
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import sys
import warnings
warnings.filterwarnings('ignore')

# Fix encoding for Windows console
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ============================================================
# LOAD DATA
# ============================================================
BASE = 'Y:/Base Rafael do Sergio'
OUTDIR = os.path.join(BASE, 'tabelas_v4')
os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(os.path.join(BASE, 'ESTUDO_C_DATASET_v4.csv'))

# Coerce numeric columns
num_cols = [c for c in df.columns if any(k in c for k in ['PESO_','IMC_','ALTURA','CREATININA_','ACIDO_URICO_',
            'EGFR_CKD_EPI_','EGFR_NI_','BSA_DUBOIS_','IDADE','AGE_AT_EXAM_','TWL_','DELTA_','N_COMORBIDADES'])]
for c in num_cols:
    df[c] = pd.to_numeric(df[c], errors='coerce')

WINDOWS = ['3M','6M','12M','24M','36M','60M']
WINDOWS_LONG = ['12M','24M','36M','60M']  # main paired windows
COHORTS = ['Arruda', 'Galvao']

# ============================================================
# HELPER FUNCTIONS
# ============================================================
def mean_sd(s):
    s = s.dropna()
    if len(s) == 0:
        return '—'
    return f"{s.mean():.2f} ± {s.std():.2f}"

def pct_sim(s):
    s = s.dropna()
    if len(s) == 0:
        return '—'
    n_sim = (s == 'SIM').sum()
    return f"{n_sim} ({100*n_sim/len(s):.1f}%)"

def ci95(arr):
    arr = arr.dropna()
    n = len(arr)
    if n == 0:
        return np.nan, np.nan, np.nan
    m = arr.mean()
    se = arr.std() / np.sqrt(n)
    return m, m - 1.96*se, m + 1.96*se

def paired_analysis(pre, post, label=''):
    """Return dict with paired stats."""
    mask = pre.notna() & post.notna()
    pre_p, post_p = pre[mask], post[mask]
    n = mask.sum()
    if n < 3:
        return {'N': n, 'PRE': '—', 'POST': '—', 'Delta': '—', 'CI95': '—', 'p_ttest': np.nan, 'p_wilcox': np.nan}
    delta = post_p - pre_p
    m, lo, hi = ci95(delta)
    try:
        _, p_t = stats.ttest_rel(pre_p, post_p)
    except:
        p_t = np.nan
    try:
        _, p_w = stats.wilcoxon(pre_p, post_p, alternative='two-sided')
    except:
        p_w = np.nan
    return {
        'N': n,
        'PRE': f"{pre_p.mean():.2f} ± {pre_p.std():.2f}",
        'POST': f"{post_p.mean():.2f} ± {post_p.std():.2f}",
        'Delta': f"{m:.2f}",
        'CI95': f"[{lo:.2f}, {hi:.2f}]",
        'p_ttest': p_t,
        'p_wilcox': p_w
    }

def fmt_p(p):
    if pd.isna(p):
        return '—'
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"

def section_header(title):
    print('\n' + '='*80)
    print(f'  {title}')
    print('='*80)

# ============================================================
# 1. TABLE 1 — BASELINE CHARACTERISTICS
# ============================================================
section_header('1. TABLE 1 — BASELINE CHARACTERISTICS')

tab1_rows = []
for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    row = {'Cohort': cohort_label, 'N': len(subset)}
    row['Age (mean±SD)'] = mean_sd(subset['IDADE_CIRURGIA'])
    fem = subset['SEXO'].dropna()
    n_fem = (fem == 'Feminino').sum()
    row['Sex (% Female)'] = f"{n_fem} ({100*n_fem/len(fem):.1f}%)" if len(fem) > 0 else '—'
    row['IMC_PRE'] = mean_sd(subset['IMC_PRE'])
    row['PESO_PRE'] = mean_sd(subset['PESO_PRE'])
    row['ALTURA_CM'] = mean_sd(subset['ALTURA_CM'])
    row['DM2 (% SIM)'] = pct_sim(subset['DM2'])
    row['HAS (% SIM)'] = pct_sim(subset['HAS'])
    row['DRC_PREVIA (% SIM)'] = pct_sim(subset['DRC_PREVIA'])
    row['Creatinina PRE'] = mean_sd(subset['CREATININA_PRE'])
    row['eGFR indexed PRE'] = mean_sd(subset['EGFR_CKD_EPI_PRE'])
    row['eGFR NI PRE'] = mean_sd(subset['EGFR_NI_PRE'])
    row['BSA PRE'] = mean_sd(subset['BSA_DUBOIS_PRE'])
    tab1_rows.append(row)

# Compare cohorts
arr = df[df.COORTE=='Arruda']
gal = df[df.COORTE=='Galvao']
comp = {}
for var in ['IDADE_CIRURGIA','IMC_PRE','PESO_PRE','ALTURA_CM','CREATININA_PRE','EGFR_CKD_EPI_PRE','EGFR_NI_PRE','BSA_DUBOIS_PRE']:
    a = arr[var].dropna()
    g = gal[var].dropna()
    if len(a) > 1 and len(g) > 1:
        _, p = stats.ttest_ind(a, g)
        comp[var] = p
for var in ['SEXO']:
    ct = pd.crosstab(df['COORTE'], df[var])
    if ct.shape[0] >= 2 and ct.shape[1] >= 2:
        _, p, _, _ = stats.chi2_contingency(ct)
        comp[var] = p
for var in ['DM2','HAS','DRC_PREVIA']:
    ct = pd.crosstab(df['COORTE'], df[var])
    if ct.shape[0] >= 2 and ct.shape[1] >= 2:
        _, p, _, _ = stats.chi2_contingency(ct)
        comp[var] = p

tab1 = pd.DataFrame(tab1_rows)
print(tab1.to_string(index=False))
print('\nCohort comparison p-values:')
for k, v in comp.items():
    print(f"  {k}: p = {fmt_p(v)}")

tab1.to_csv(os.path.join(OUTDIR, 'Tab1_baseline.csv'), index=False)

# ============================================================
# 2. PAIRED ANALYSIS — ALL WINDOWS
# ============================================================
section_header('2. PAIRED ANALYSIS — ALL WINDOWS')

paired_rows = []
for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    for w in WINDOWS:
        # Creatinina
        cr_pre = subset['CREATININA_PRE']
        cr_post = subset.get(f'CREATININA_{w}')
        if cr_post is not None:
            res = paired_analysis(cr_pre, cr_post)
            res.update({'Cohort': cohort_label, 'Window': w, 'Variable': 'Creatinina'})
            paired_rows.append(res)

        # eGFR indexed
        eg_pre = subset['EGFR_CKD_EPI_PRE']
        eg_post = subset.get(f'EGFR_CKD_EPI_{w}')
        if eg_post is not None:
            res = paired_analysis(eg_pre, eg_post)
            res.update({'Cohort': cohort_label, 'Window': w, 'Variable': 'eGFR_indexed'})
            paired_rows.append(res)

        # eGFR nonindexed (observed weight only)
        ni_pre = subset['EGFR_NI_PRE']
        ni_post = subset.get(f'EGFR_NI_{w}')
        if ni_post is not None:
            res = paired_analysis(ni_pre, ni_post)
            res.update({'Cohort': cohort_label, 'Window': w, 'Variable': 'eGFR_NI'})
            paired_rows.append(res)

        # BSA
        bsa_pre = subset['BSA_DUBOIS_PRE']
        bsa_post = subset.get(f'BSA_DUBOIS_{w}')
        if bsa_post is not None:
            res = paired_analysis(bsa_pre, bsa_post)
            res.update({'Cohort': cohort_label, 'Window': w, 'Variable': 'BSA'})
            paired_rows.append(res)

        # Acido urico
        au_pre = subset['ACIDO_URICO_PRE']
        au_post = subset.get(f'ACIDO_URICO_{w}')
        if au_post is not None:
            res = paired_analysis(au_pre, au_post)
            res.update({'Cohort': cohort_label, 'Window': w, 'Variable': 'Acido_Urico'})
            paired_rows.append(res)

paired_df = pd.DataFrame(paired_rows)
paired_df['p_ttest'] = paired_df['p_ttest'].apply(fmt_p)
paired_df['p_wilcox'] = paired_df['p_wilcox'].apply(fmt_p)
display_cols = ['Cohort','Window','Variable','N','PRE','POST','Delta','CI95','p_ttest','p_wilcox']
for coh in ['Arruda','Galvao','Combined']:
    print(f"\n--- {coh} ---")
    sub = paired_df[paired_df.Cohort == coh][display_cols]
    print(sub.to_string(index=False))

paired_df[display_cols].to_csv(os.path.join(OUTDIR, 'Tab2_paired_overall.csv'), index=False)

# ============================================================
# 3. DIVERGENCE TEST
# ============================================================
section_header('3. DIVERGENCE TEST (Δ_indexed − Δ_nonindexed)')

div_rows = []
for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    for w in WINDOWS:
        idx_pre = subset['EGFR_CKD_EPI_PRE']
        idx_post = subset.get(f'EGFR_CKD_EPI_{w}')
        ni_pre = subset['EGFR_NI_PRE']
        ni_post = subset.get(f'EGFR_NI_{w}')
        if idx_post is None or ni_post is None:
            continue
        mask = idx_pre.notna() & idx_post.notna() & ni_pre.notna() & ni_post.notna()
        if mask.sum() < 3:
            div_rows.append({'Cohort': cohort_label, 'Window': w, 'N': mask.sum(),
                             'Mean_Div': '—', 'SD': '—', 'p': '—'})
            continue
        d_idx = idx_post[mask].values - idx_pre[mask].values
        d_ni = ni_post[mask].values - ni_pre[mask].values
        divergence = d_idx - d_ni
        m = np.nanmean(divergence)
        sd = np.nanstd(divergence, ddof=1)
        try:
            _, p = stats.ttest_1samp(divergence, 0)
        except:
            p = np.nan
        div_rows.append({'Cohort': cohort_label, 'Window': w, 'N': mask.sum(),
                         'Mean_Div': f"{m:.2f}", 'SD': f"{sd:.2f}", 'p': fmt_p(p)})

div_df = pd.DataFrame(div_rows)
print(div_df.to_string(index=False))
div_df.to_csv(os.path.join(OUTDIR, 'Tab3_divergence.csv'), index=False)

# ============================================================
# 4. DISCORDANCE ANALYSIS
# ============================================================
section_header('4. DISCORDANCE ANALYSIS')

disc_rows = []
for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    for w in WINDOWS:
        idx_pre = subset['EGFR_CKD_EPI_PRE']
        idx_post = subset.get(f'EGFR_CKD_EPI_{w}')
        ni_pre = subset['EGFR_NI_PRE']
        ni_post = subset.get(f'EGFR_NI_{w}')
        if idx_post is None or ni_post is None:
            continue
        mask = idx_pre.notna() & idx_post.notna() & ni_pre.notna() & ni_post.notna()
        n = mask.sum()
        if n < 1:
            continue
        d_idx = idx_post[mask].values - idx_pre[mask].values
        d_ni = ni_post[mask].values - ni_pre[mask].values
        disc_pos = ((d_idx > 0) & (d_ni < 0)).sum()  # idx up, NI down
        disc_neg = ((d_idx < 0) & (d_ni > 0)).sum()  # idx down, NI up
        concordant = n - disc_pos - disc_neg
        disc_rows.append({
            'Cohort': cohort_label, 'Window': w, 'N': n,
            'Concordant': f"{concordant} ({100*concordant/n:.1f}%)",
            'Disc_Pos (idx↑ NI↓)': f"{disc_pos} ({100*disc_pos/n:.1f}%)",
            'Disc_Neg (idx↓ NI↑)': f"{disc_neg} ({100*disc_neg/n:.1f}%)"
        })

disc_df = pd.DataFrame(disc_rows)
print(disc_df.to_string(index=False))

# ============================================================
# 5. SHAPLEY DECOMPOSITION
# ============================================================
section_header('5. SHAPLEY DECOMPOSITION')

shap_rows = []
for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    for w in WINDOWS:
        egfr_pre = subset['EGFR_CKD_EPI_PRE']
        egfr_post = subset.get(f'EGFR_CKD_EPI_{w}')
        bsa_pre = subset['BSA_DUBOIS_PRE']
        bsa_post = subset.get(f'BSA_DUBOIS_{w}')
        ni_pre = subset['EGFR_NI_PRE']
        ni_post = subset.get(f'EGFR_NI_{w}')
        if egfr_post is None or bsa_post is None or ni_post is None:
            continue
        mask = egfr_pre.notna() & egfr_post.notna() & bsa_pre.notna() & bsa_post.notna() & ni_pre.notna() & ni_post.notna()
        n = mask.sum()
        if n < 3:
            continue
        ep = egfr_pre[mask].values
        eo = egfr_post[mask].values
        bp = bsa_pre[mask].values
        bo = bsa_post[mask].values
        np_ = ni_pre[mask].values
        no_ = ni_post[mask].values

        egfr_mean = (ep + eo) / 2
        bsa_mean = (bp + bo) / 2
        d_egfr = eo - ep
        d_bsa = bo - bp
        d_ni = no_ - np_

        comp_egfr = d_egfr * bsa_mean / 1.73
        comp_bsa = egfr_mean * d_bsa / 1.73
        # pct BSA
        with np.errstate(divide='ignore', invalid='ignore'):
            pct_bsa_arr = np.where(d_ni != 0, comp_bsa / d_ni * 100, np.nan)

        shap_rows.append({
            'Cohort': cohort_label, 'Window': w, 'N': n,
            'Delta_NI': f"{np.nanmean(d_ni):.2f} ± {np.nanstd(d_ni, ddof=1):.2f}",
            'Comp_eGFR': f"{np.nanmean(comp_egfr):.2f} ± {np.nanstd(comp_egfr, ddof=1):.2f}",
            'Comp_BSA': f"{np.nanmean(comp_bsa):.2f} ± {np.nanstd(comp_bsa, ddof=1):.2f}",
            '%BSA (mean)': f"{np.nanmean(pct_bsa_arr):.1f}%",
            '%BSA (median)': f"{np.nanmedian(pct_bsa_arr):.1f}%"
        })

shap_df = pd.DataFrame(shap_rows)
print(shap_df.to_string(index=False))
shap_df.to_csv(os.path.join(OUTDIR, 'Tab4_shapley.csv'), index=False)

# ============================================================
# 6. KDIGO G-CATEGORY TRANSITIONS
# ============================================================
section_header('6. KDIGO G-CATEGORY TRANSITIONS')

kdigo_order = ['G1','G2','G3a','G3b','G4','G5']
kdigo_rank = {g: i for i, g in enumerate(kdigo_order)}

kdigo_rows = []
for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    for w in WINDOWS:
        kcol = f'KDIGO_{w}'
        if kcol not in subset.columns:
            continue
        mask = subset['KDIGO_PRE'].notna() & subset[kcol].notna()
        n = mask.sum()
        if n < 1:
            continue
        pre_k = subset.loc[mask, 'KDIGO_PRE']
        post_k = subset.loc[mask, kcol]
        pre_rank = pre_k.map(kdigo_rank)
        post_rank = post_k.map(kdigo_rank)
        improved = (post_rank < pre_rank).sum()  # lower rank = better
        stable = (post_rank == pre_rank).sum()
        worsened = (post_rank > pre_rank).sum()
        kdigo_rows.append({
            'Cohort': cohort_label, 'Window': w, 'N': n,
            'Improved': f"{improved} ({100*improved/n:.1f}%)",
            'Stable': f"{stable} ({100*stable/n:.1f}%)",
            'Worsened': f"{worsened} ({100*worsened/n:.1f}%)"
        })
        # Cross-tab
        ct = pd.crosstab(pre_k, post_k, margins=True)
        print(f"\n{cohort_label} — KDIGO PRE vs {w}:")
        print(ct.to_string())

kdigo_df = pd.DataFrame(kdigo_rows)
print('\nSummary:')
print(kdigo_df.to_string(index=False))
kdigo_df.to_csv(os.path.join(OUTDIR, 'Tab5_kdigo_transitions.csv'), index=False)

# ============================================================
# 7. SUBGROUP BY BASELINE RENAL GROUP (eGFR indexed PRE+12M)
# ============================================================
section_header('7. SUBGROUP BY BASELINE RENAL GROUP')

for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    print(f"\n--- {cohort_label} ---")
    mask = subset['EGFR_CKD_EPI_PRE'].notna() & subset['EGFR_CKD_EPI_12M'].notna() & subset['GRUPO_RENAL'].notna()
    sub = subset[mask]
    for grp in ['<60','60-89','90-120','>120']:
        g = sub[sub.GRUPO_RENAL == grp]
        n = len(g)
        if n < 3:
            print(f"  {grp}: N={n} (too few)")
            continue
        pre = g['EGFR_CKD_EPI_PRE']
        post = g['EGFR_CKD_EPI_12M']
        delta = post - pre
        m, lo, hi = ci95(delta)
        try:
            _, p = stats.ttest_rel(pre, post)
        except:
            p = np.nan
        print(f"  {grp}: N={n}, PRE={pre.mean():.1f}±{pre.std():.1f}, 12M={post.mean():.1f}±{post.std():.1f}, "
              f"Δ={m:.2f} [{lo:.2f},{hi:.2f}], p={fmt_p(p)}")

# ============================================================
# 8. HYPERFILTERED (>120) TRAJECTORY
# ============================================================
section_header('8. HYPERFILTERED (>120) TRAJECTORY')

for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    print(f"\n--- {cohort_label} ---")
    hyper = subset[subset['EGFR_CKD_EPI_PRE'] > 120]
    print(f"  N with eGFR PRE > 120: {len(hyper)}")
    for w in ['PRE'] + WINDOWS:
        col = f'EGFR_CKD_EPI_{w}'
        if col in hyper.columns:
            vals = hyper[col].dropna()
            if len(vals) > 0:
                print(f"  {w}: N={len(vals)}, mean={vals.mean():.2f} ± {vals.std():.2f}")

# ============================================================
# 9. URIC ACID — SUSTAINED REDUCTION
# ============================================================
section_header('9. URIC ACID — SUSTAINED REDUCTION')

for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    print(f"\n--- {cohort_label} ---")
    for w in WINDOWS:
        au_col = f'ACIDO_URICO_{w}'
        if au_col not in subset.columns:
            continue
        res = paired_analysis(subset['ACIDO_URICO_PRE'], subset[au_col])
        print(f"  {w}: N={res['N']}, PRE={res['PRE']}, POST={res['POST']}, Δ={res['Delta']} {res['CI95']}, "
              f"p_t={fmt_p(res['p_ttest'])}, p_w={fmt_p(res['p_wilcox'])}")

# ============================================================
# 10. CORRELATION WITH WEIGHT LOSS (12M)
# ============================================================
section_header('10. CORRELATION WITH WEIGHT LOSS (12M)')

for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    print(f"\n--- {cohort_label} ---")
    for x_var, y_var in [('TWL_12M','DELTA_EGFR_IDX_12M'), ('TWL_12M','DELTA_EGFR_NI_12M'), ('DELTA_BSA_12M','DELTA_EGFR_NI_12M')]:
        mask = subset[x_var].notna() & subset[y_var].notna()
        n = mask.sum()
        if n < 5:
            print(f"  {x_var} vs {y_var}: N={n} (too few)")
            continue
        r, p = stats.pearsonr(subset.loc[mask, x_var], subset.loc[mask, y_var])
        print(f"  {x_var} vs {y_var}: N={n}, r={r:.4f}, p={fmt_p(p)}")

# ============================================================
# 11. FOREST PLOT DATA
# ============================================================
section_header('11. FOREST PLOT DATA (DISCOVERY vs REPLICATION)')

forest_rows = []
for endpoint_name, pre_col, post_col in [
    ('Delta_eGFR_idx_12M', 'EGFR_CKD_EPI_PRE', 'EGFR_CKD_EPI_12M'),
    ('Delta_eGFR_NI_12M', 'EGFR_NI_PRE', 'EGFR_NI_12M'),
    ('Delta_BSA_12M', 'BSA_DUBOIS_PRE', 'BSA_DUBOIS_12M'),
]:
    cohort_effects = {}
    for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao'])]:
        mask = subset[pre_col].notna() & subset[post_col].notna()
        delta = subset.loc[mask, post_col].values - subset.loc[mask, pre_col].values
        n = len(delta)
        if n < 3:
            continue
        m, lo, hi = ci95(pd.Series(delta))
        cohort_effects[cohort_label] = {'n': n, 'mean': m, 'lo': lo, 'hi': hi, 'delta': delta}
        forest_rows.append({
            'Endpoint': endpoint_name, 'Cohort': cohort_label, 'N': n,
            'Mean': f"{m:.2f}", 'CI95': f"[{lo:.2f}, {hi:.2f}]"
        })
    # Interaction test
    if len(cohort_effects) == 2:
        _, p_int = stats.ttest_ind(cohort_effects['Arruda']['delta'], cohort_effects['Galvao']['delta'])
        forest_rows[-1]['p_interaction'] = fmt_p(p_int)
        forest_rows[-2]['p_interaction'] = fmt_p(p_int)

# Divergence (Δ_idx - Δ_NI) per cohort
cohort_div = {}
for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao'])]:
    mask = (subset['EGFR_CKD_EPI_PRE'].notna() & subset['EGFR_CKD_EPI_12M'].notna() &
            subset['EGFR_NI_PRE'].notna() & subset['EGFR_NI_12M'].notna())
    d_idx = subset.loc[mask, 'EGFR_CKD_EPI_12M'].values - subset.loc[mask, 'EGFR_CKD_EPI_PRE'].values
    d_ni = subset.loc[mask, 'EGFR_NI_12M'].values - subset.loc[mask, 'EGFR_NI_PRE'].values
    div = d_idx - d_ni
    n = len(div)
    if n >= 3:
        m, lo, hi = ci95(pd.Series(div))
        cohort_div[cohort_label] = div
        forest_rows.append({
            'Endpoint': 'Divergence_12M', 'Cohort': cohort_label, 'N': n,
            'Mean': f"{m:.2f}", 'CI95': f"[{lo:.2f}, {hi:.2f}]"
        })
if len(cohort_div) == 2:
    _, p_int = stats.ttest_ind(cohort_div['Arruda'], cohort_div['Galvao'])
    forest_rows[-1]['p_interaction'] = fmt_p(p_int)
    forest_rows[-2]['p_interaction'] = fmt_p(p_int)

forest_df = pd.DataFrame(forest_rows)
print(forest_df.to_string(index=False))
forest_df.to_csv(os.path.join(OUTDIR, 'Tab6_forest_plot_data.csv'), index=False)

# ============================================================
# 12. SENSITIVITY ANALYSES
# ============================================================
section_header('12. SENSITIVITY ANALYSES')

sens_rows = []

def sensitivity_paired_egfr(label, subset, window='12M'):
    """Run paired eGFR indexed and NI for a subset at given window."""
    results = []
    for var_name, pre_col, post_col in [
        ('eGFR_idx', 'EGFR_CKD_EPI_PRE', f'EGFR_CKD_EPI_{window}'),
        ('eGFR_NI', 'EGFR_NI_PRE', f'EGFR_NI_{window}'),
    ]:
        if post_col not in subset.columns:
            continue
        res = paired_analysis(subset[pre_col], subset[post_col])
        res.update({'Sensitivity': label, 'Variable': var_name, 'Window': window})
        results.append(res)
    return results

# a) Already built into all analyses (Arruda vs Galvao)

# b) By sex
print('\n--- By Sex (Combined, 12M) ---')
for sex in ['Feminino', 'Masculino']:
    sub = df[df.SEXO == sex]
    res_list = sensitivity_paired_egfr(f'Sex={sex}', sub)
    for r in res_list:
        sens_rows.append(r)
        print(f"  {sex} | {r['Variable']}: N={r['N']}, Δ={r['Delta']} {r['CI95']}, p_t={fmt_p(r['p_ttest'])}")

# c) By IMC group
print('\n--- By IMC Group (Combined, 12M) ---')
for imc_grp in [1.0, 2.0, 3.0]:
    sub = df[df.GRUPO_IMC == imc_grp]
    res_list = sensitivity_paired_egfr(f'IMC_grp={int(imc_grp)}', sub)
    for r in res_list:
        sens_rows.append(r)
        print(f"  IMC={int(imc_grp)} | {r['Variable']}: N={r['N']}, Δ={r['Delta']} {r['CI95']}, p_t={fmt_p(r['p_ttest'])}")

# d) Exclude eGFR < 60 at baseline
print('\n--- Exclude eGFR < 60 at baseline (Combined, 12M) ---')
sub = df[~(df.EGFR_CKD_EPI_PRE < 60)]
res_list = sensitivity_paired_egfr('Excl_eGFR<60', sub)
for r in res_list:
    sens_rows.append(r)
    print(f"  {r['Variable']}: N={r['N']}, Δ={r['Delta']} {r['CI95']}, p_t={fmt_p(r['p_ttest'])}")

# e) Restrict creatinine 0.4-3.0
print('\n--- Creatinine 0.4-3.0 restriction (Combined, 12M) ---')
sub = df[(df.CREATININA_PRE >= 0.4) & (df.CREATININA_PRE <= 3.0)]
# Also filter post if available
if 'CREATININA_12M' in sub.columns:
    sub = sub[(sub.CREATININA_12M.isna()) | ((sub.CREATININA_12M >= 0.4) & (sub.CREATININA_12M <= 3.0))]
res_list = sensitivity_paired_egfr('Cr_0.4-3.0', sub)
for r in res_list:
    sens_rows.append(r)
    print(f"  {r['Variable']}: N={r['N']}, Δ={r['Delta']} {r['CI95']}, p_t={fmt_p(r['p_ttest'])}")

sens_df = pd.DataFrame(sens_rows)
sens_df['p_ttest'] = sens_df['p_ttest'].apply(fmt_p)
sens_df['p_wilcox'] = sens_df['p_wilcox'].apply(fmt_p)
sens_df.to_csv(os.path.join(OUTDIR, 'Tab7_sensitivity.csv'), index=False)

# ============================================================
# 13. POWER ANALYSIS
# ============================================================
section_header('13. POWER ANALYSIS (Achieved Power)')

from scipy.stats import norm

def achieved_power(n, effect_size, alpha=0.05):
    """Compute achieved power for paired t-test."""
    if n < 2 or effect_size == 0 or np.isnan(effect_size):
        return np.nan
    z_alpha = norm.ppf(1 - alpha/2)
    z_power = abs(effect_size) * np.sqrt(n) - z_alpha
    return norm.cdf(z_power)

print(f"{'Cohort':<12} {'Window':<8} {'Variable':<15} {'N':>5} {'Effect_d':>10} {'Power':>8}")
print('-' * 65)

for cohort_label, subset in [('Arruda', df[df.COORTE=='Arruda']), ('Galvao', df[df.COORTE=='Galvao']), ('Combined', df)]:
    for w in WINDOWS_LONG:
        for var_name, pre_col, post_col in [
            ('eGFR_idx', 'EGFR_CKD_EPI_PRE', f'EGFR_CKD_EPI_{w}'),
            ('eGFR_NI', 'EGFR_NI_PRE', f'EGFR_NI_{w}'),
        ]:
            if post_col not in subset.columns:
                continue
            mask = subset[pre_col].notna() & subset[post_col].notna()
            n = mask.sum()
            if n < 3:
                continue
            delta = subset.loc[mask, post_col].values - subset.loc[mask, pre_col].values
            d = np.mean(delta) / np.std(delta, ddof=1) if np.std(delta, ddof=1) > 0 else 0
            pwr = achieved_power(n, d)
            print(f"{cohort_label:<12} {w:<8} {var_name:<15} {n:>5} {d:>10.4f} {pwr:>8.4f}")

# ============================================================
# FINAL SUMMARY
# ============================================================
section_header('ANALYSIS COMPLETE')
print(f"Dataset: {len(df)} patients ({df.COORTE.value_counts().to_dict()})")
print(f"Output tables saved to: {OUTDIR}/")
for f in sorted(os.listdir(OUTDIR)):
    fpath = os.path.join(OUTDIR, f)
    print(f"  {f} ({os.path.getsize(fpath)} bytes)")
print("\nDone.")
