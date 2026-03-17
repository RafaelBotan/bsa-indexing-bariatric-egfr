"""
Análises Complementares v4 — Estudo C (eGFR indexada vs não-indexada)
Tab8:  Dose-response TWL quartiles
Tab8b: Dose-response ΔBSA quartiles
Tab9:  Mosteller sensitivity analysis
Tab10: Representativeness of observed-weight subset
Tab11: Linked vs Unlinked (Galvão Sabin)
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import spearmanr, pearsonr, ttest_ind, ttest_rel, chi2_contingency
from sklearn.linear_model import LinearRegression
import warnings, os, re, sys, io
from datetime import datetime

# Fix encoding for Windows console
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

warnings.filterwarnings('ignore')

OUT = r'Y:\Base Rafael do Sergio\tabelas_v4'
os.makedirs(OUT, exist_ok=True)

# ─── Load main dataset ───
df = pd.read_csv(r'Y:\Base Rafael do Sergio\ESTUDO_C_DATASET_v4.csv')
print(f"Dataset: {len(df)} rows, {len(df.columns)} cols")
print(f"Coortes: {df['COORTE'].value_counts().to_dict()}")

# Primary sample
am = df[df['AMOSTRA_PRIMARIA_12M'] == True].copy()
print(f"AMOSTRA_PRIMARIA_12M: {len(am)} patients")
print(f"  Arruda: {(am['COORTE']=='Arruda').sum()}, Galvao: {(am['COORTE']=='Galvao').sum()}")

# Divergence
am['DIVERGENCE_12M'] = am['DELTA_EGFR_IDX_12M'] - am['DELTA_EGFR_NI_12M']


# ═══════════════════════════════════════════════════════════════
# ANALYSIS 1: DOSE-RESPONSE
# ═══════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 1a: DOSE-RESPONSE — TWL QUARTILES vs DIVERGENCE")
print("="*80)

def dose_response_analysis(data, split_var, label, cohort_label):
    """Compute quartile-based dose-response table."""
    sub = data.dropna(subset=[split_var, 'DIVERGENCE_12M']).copy()
    if len(sub) < 20:
        return None, None
    sub['Q'] = pd.qcut(sub[split_var], 4, labels=['Q1','Q2','Q3','Q4'], duplicates='drop')

    rows = []
    for q in ['Q1','Q2','Q3','Q4']:
        g = sub[sub['Q'] == q]
        if len(g) == 0:
            continue
        rows.append({
            'Cohort': cohort_label,
            'Quartile': q,
            'N': len(g),
            f'Mean_{split_var}': g[split_var].mean(),
            'Mean_DELTA_EGFR_IDX_12M': g['DELTA_EGFR_IDX_12M'].mean(),
            'Mean_DELTA_EGFR_NI_12M': g['DELTA_EGFR_NI_12M'].mean(),
            'Mean_Divergence': g['DIVERGENCE_12M'].mean(),
            'Mean_DELTA_BSA_12M': g['DELTA_BSA_12M'].mean(),
        })

    # Trend tests
    valid = sub.dropna(subset=[split_var, 'DIVERGENCE_12M'])
    rho, p_spear = spearmanr(valid[split_var], valid['DIVERGENCE_12M'])

    X = valid[split_var].values.reshape(-1,1)
    y = valid['DIVERGENCE_12M'].values
    reg = LinearRegression().fit(X, y)
    slope = reg.coef_[0]
    r2 = reg.score(X, y)
    # p-value for slope via pearsonr
    r_pear, p_pear = pearsonr(valid[split_var], valid['DIVERGENCE_12M'])

    trend = {
        'Cohort': cohort_label,
        'Variable': split_var,
        'N': len(valid),
        'Spearman_rho': rho,
        'Spearman_p': p_spear,
        'LR_slope': slope,
        'LR_R2': r2,
        'Pearson_r': r_pear,
        'Pearson_p': p_pear,
    }

    return pd.DataFrame(rows), trend

# TWL quartiles
all_tab8 = []
all_trend8 = []
for cohort_label, cohort_data in [('Arruda', am[am['COORTE']=='Arruda']),
                                    ('Galvao', am[am['COORTE']=='Galvao']),
                                    ('Combined', am)]:
    tbl, trend = dose_response_analysis(cohort_data, 'TWL_12M', 'TWL', cohort_label)
    if tbl is not None:
        all_tab8.append(tbl)
        all_trend8.append(trend)

tab8 = pd.concat(all_tab8, ignore_index=True)
tab8_trend = pd.DataFrame(all_trend8)
tab8_full = pd.concat([tab8, pd.DataFrame([{}]), tab8_trend], ignore_index=True)
tab8_full.to_csv(os.path.join(OUT, 'Tab8_dose_response_TWL.csv'), index=False)

print("\nQuartile Table:")
print(tab8.to_string(index=False, float_format='%.3f'))
print("\nTrend Tests:")
print(tab8_trend.to_string(index=False, float_format='%.4f'))

# BSA quartiles
print("\n" + "="*80)
print("ANALYSIS 1b: DOSE-RESPONSE — ΔBSA QUARTILES vs DIVERGENCE")
print("="*80)

all_tab8b = []
all_trend8b = []
for cohort_label, cohort_data in [('Arruda', am[am['COORTE']=='Arruda']),
                                    ('Galvao', am[am['COORTE']=='Galvao']),
                                    ('Combined', am)]:
    tbl, trend = dose_response_analysis(cohort_data, 'DELTA_BSA_12M', 'ΔBSA', cohort_label)
    if tbl is not None:
        all_tab8b.append(tbl)
        all_trend8b.append(trend)

tab8b = pd.concat(all_tab8b, ignore_index=True)
tab8b_trend = pd.DataFrame(all_trend8b)
tab8b_full = pd.concat([tab8b, pd.DataFrame([{}]), tab8b_trend], ignore_index=True)
tab8b_full.to_csv(os.path.join(OUT, 'Tab8b_dose_response_BSA.csv'), index=False)

print("\nQuartile Table:")
print(tab8b.to_string(index=False, float_format='%.3f'))
print("\nTrend Tests:")
print(tab8b_trend.to_string(index=False, float_format='%.4f'))


# ═══════════════════════════════════════════════════════════════
# ANALYSIS 2: MOSTELLER SENSITIVITY
# ═══════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 2: BSA MOSTELLER SENSITIVITY")
print("="*80)

am2 = am.dropna(subset=['ALTURA_CM','PESO_PRE','PESO_12M','EGFR_CKD_EPI_PRE','EGFR_CKD_EPI_12M']).copy()

# Mosteller BSA = sqrt(height_cm * weight_kg / 3600)
am2['BSA_MOST_PRE'] = np.sqrt(am2['ALTURA_CM'] * am2['PESO_PRE'] / 3600)
am2['BSA_MOST_12M'] = np.sqrt(am2['ALTURA_CM'] * am2['PESO_12M'] / 3600)

# Non-indexed eGFR via Mosteller
am2['EGFR_NI_MOST_PRE'] = am2['EGFR_CKD_EPI_PRE'] * am2['BSA_MOST_PRE'] / 1.73
am2['EGFR_NI_MOST_12M'] = am2['EGFR_CKD_EPI_12M'] * am2['BSA_MOST_12M'] / 1.73
am2['DELTA_NI_MOST'] = am2['EGFR_NI_MOST_12M'] - am2['EGFR_NI_MOST_PRE']

# DuBois delta (already in dataset)
am2['DELTA_NI_DUBOIS'] = am2['DELTA_EGFR_NI_12M']

rows9 = []
for cohort_label, cohort_data in [('Arruda', am2[am2['COORTE']=='Arruda']),
                                    ('Galvao', am2[am2['COORTE']=='Galvao']),
                                    ('Combined', am2)]:
    c = cohort_data.dropna(subset=['DELTA_NI_MOST','DELTA_NI_DUBOIS'])
    if len(c) < 5:
        continue

    t_stat, p_paired = ttest_rel(c['DELTA_NI_DUBOIS'], c['DELTA_NI_MOST'])

    # Shapley-like decomposition for Mosteller:
    # Total Δ indexed = EGFR_CKD_EPI_12M - EGFR_CKD_EPI_PRE
    # Δ NI Mosteller = EGFR_NI_MOST_12M - EGFR_NI_MOST_PRE
    # "BSA component" ≈ Δ indexed - Δ NI / (BSA_factor)
    # Simple Shapley: Δindexed = Δ_true_GFR + Δ_BSA_artifact
    # Where Δ_true_GFR ≈ Δ NI (corrected), Δ_BSA_artifact ≈ Δidx - Δ NI

    delta_idx = c['DELTA_EGFR_IDX_12M']
    delta_ni_most = c['DELTA_NI_MOST']
    delta_ni_dub = c['DELTA_NI_DUBOIS']

    # BSA changes
    delta_bsa_dub = c['DELTA_BSA_12M']
    delta_bsa_most = c['BSA_MOST_12M'] - c['BSA_MOST_PRE']

    rows9.append({
        'Cohort': cohort_label,
        'N': len(c),
        'Mean_BSA_DuBois_PRE': c['BSA_DUBOIS_PRE'].mean(),
        'Mean_BSA_Most_PRE': c['BSA_MOST_PRE'].mean(),
        'Mean_BSA_DuBois_12M': c['BSA_DUBOIS_12M'].mean(),
        'Mean_BSA_Most_12M': c['BSA_MOST_12M'].mean(),
        'Mean_ΔBSA_DuBois': delta_bsa_dub.mean(),
        'Mean_ΔBSA_Mosteller': delta_bsa_most.mean(),
        'Mean_Δ_eGFR_IDX': delta_idx.mean(),
        'Mean_Δ_NI_DuBois': delta_ni_dub.mean(),
        'Mean_Δ_NI_Mosteller': delta_ni_most.mean(),
        'Diff_DuBois_vs_Most': (delta_ni_dub - delta_ni_most).mean(),
        'Paired_t': t_stat,
        'Paired_p': p_paired,
        'Shapley_TrueGFR_DuBois': delta_ni_dub.mean(),
        'Shapley_BSA_artifact_DuBois': delta_idx.mean() - delta_ni_dub.mean(),
        'Shapley_TrueGFR_Mosteller': delta_ni_most.mean(),
        'Shapley_BSA_artifact_Mosteller': delta_idx.mean() - delta_ni_most.mean(),
    })

tab9 = pd.DataFrame(rows9)
tab9.to_csv(os.path.join(OUT, 'Tab9_mosteller_sensitivity.csv'), index=False)

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)
print(tab9.to_string(index=False, float_format='%.4f'))


# ═══════════════════════════════════════════════════════════════
# ANALYSIS 3: REPRESENTATIVENESS
# ═══════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 3: REPRESENTATIVENESS OF OBSERVED-WEIGHT SUBSET")
print("="*80)

# Group A: AMOSTRA_PRIMARIA_12M = True
# Group B: Has paired indexed 12M eGFR but NOT in primary sample
df['HAS_PAIRED_IDX_12M'] = df['EGFR_CKD_EPI_PRE'].notna() & df['EGFR_CKD_EPI_12M'].notna()

groupA = df[(df['AMOSTRA_PRIMARIA_12M'] == True)].copy()
groupB = df[(df['HAS_PAIRED_IDX_12M'] == True) & (df['AMOSTRA_PRIMARIA_12M'] != True)].copy()

print(f"\nGroup A (primary sample): {len(groupA)}")
print(f"Group B (indexed only): {len(groupB)}")

def compute_smd(a, b):
    """Standardized Mean Difference (Cohen's d with pooled SD)."""
    na, nb = a.dropna(), b.dropna()
    if len(na) < 2 or len(nb) < 2:
        return np.nan
    pooled_sd = np.sqrt(((len(na)-1)*na.std()**2 + (len(nb)-1)*nb.std()**2) / (len(na)+len(nb)-2))
    if pooled_sd == 0:
        return 0
    return (na.mean() - nb.mean()) / pooled_sd

def representativeness_compare(gA, gB, cohort_label):
    """Compare two groups on key variables."""
    continuous_vars = [
        ('Age', 'IDADE_CIRURGIA'),
        ('IMC_PRE', 'IMC_PRE'),
        ('PESO_PRE', 'PESO_PRE'),
        ('Creatinina_PRE', 'CREATININA_PRE'),
        ('eGFR_IDX_PRE', 'EGFR_CKD_EPI_PRE'),
    ]
    categorical_vars = [
        ('Sex_Fem_%', 'SEXO', 'Feminino'),
        ('DM2_%', 'DM2', 'SIM'),
        ('HAS_%', 'HAS', 'SIM'),
    ]

    rows = []
    for label, col in continuous_vars:
        a_vals = gA[col].dropna()
        b_vals = gB[col].dropna()

        smd = compute_smd(a_vals, b_vals)

        if len(a_vals) > 1 and len(b_vals) > 1:
            t, p = ttest_ind(a_vals, b_vals, equal_var=False)
        else:
            t, p = np.nan, np.nan

        rows.append({
            'Cohort': cohort_label,
            'Variable': label,
            'GroupA_N': len(a_vals),
            'GroupA_Mean': a_vals.mean(),
            'GroupA_SD': a_vals.std(),
            'GroupB_N': len(b_vals),
            'GroupB_Mean': b_vals.mean(),
            'GroupB_SD': b_vals.std(),
            'SMD': smd,
            'Test': 't-test',
            'Statistic': t,
            'P_value': p,
            'Imbalance_flag': 'YES' if abs(smd) > 0.20 else '',
        })

    for label, col, val in categorical_vars:
        a_n = len(gA[col].dropna())
        b_n = len(gB[col].dropna())
        a_count = (gA[col] == val).sum()
        b_count = (gB[col] == val).sum()
        a_pct = a_count / a_n * 100 if a_n > 0 else np.nan
        b_pct = b_count / b_n * 100 if b_n > 0 else np.nan

        # SMD for proportions
        p_a = a_count / a_n if a_n > 0 else 0
        p_b = b_count / b_n if b_n > 0 else 0
        pooled_p = (a_count + b_count) / (a_n + b_n) if (a_n + b_n) > 0 else 0
        denom = np.sqrt(pooled_p * (1 - pooled_p)) if pooled_p > 0 and pooled_p < 1 else 1
        smd_cat = (p_a - p_b) / denom if denom > 0 else 0

        # Chi-square
        if a_n > 0 and b_n > 0 and a_count > 0:
            table = [[a_count, a_n - a_count], [b_count, b_n - b_count]]
            try:
                chi2, p_chi, _, _ = chi2_contingency(table)
            except:
                chi2, p_chi = np.nan, np.nan
        else:
            chi2, p_chi = np.nan, np.nan

        rows.append({
            'Cohort': cohort_label,
            'Variable': label,
            'GroupA_N': a_n,
            'GroupA_Mean': a_pct,
            'GroupA_SD': np.nan,
            'GroupB_N': b_n,
            'GroupB_Mean': b_pct,
            'GroupB_SD': np.nan,
            'SMD': smd_cat,
            'Test': 'chi2',
            'Statistic': chi2,
            'P_value': p_chi,
            'Imbalance_flag': 'YES' if abs(smd_cat) > 0.20 else '',
        })

    return pd.DataFrame(rows)

all_tab10 = []
for cohort_label in ['Arruda', 'Galvao', 'Combined']:
    if cohort_label == 'Combined':
        gA = groupA
        gB = groupB
    else:
        gA = groupA[groupA['COORTE'] == cohort_label]
        gB = groupB[groupB['COORTE'] == cohort_label]

    print(f"\n{cohort_label}: Group A={len(gA)}, Group B={len(gB)}")
    tbl = representativeness_compare(gA, gB, cohort_label)
    all_tab10.append(tbl)

tab10 = pd.concat(all_tab10, ignore_index=True)
tab10.to_csv(os.path.join(OUT, 'Tab10_representativeness.csv'), index=False)

print("\n" + tab10.to_string(index=False, float_format='%.3f'))

# Flag imbalances
flagged = tab10[tab10['Imbalance_flag'] == 'YES']
if len(flagged) > 0:
    print(f"\n⚠ IMBALANCES DETECTED (SMD > 0.20):")
    for _, r in flagged.iterrows():
        print(f"  {r['Cohort']} — {r['Variable']}: SMD = {r['SMD']:.3f}")
else:
    print("\n✓ No meaningful imbalances detected (all SMD ≤ 0.20)")


# ═══════════════════════════════════════════════════════════════
# ANALYSIS 4: LINKED vs UNLINKED (Galvão Sabin)
# ═══════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("ANALYSIS 4: LINKED vs UNLINKED IN GALVÃO SABIN")
print("="*80)

mega = pd.read_csv(r'Y:\Base Rafael do Sergio\MEGA_BASE_GALVAO.csv', encoding='utf-8-sig')
link = pd.read_csv(r'Y:\Base Rafael do Sergio\LINKAGEM_GALVAO_CLINICA.csv')
mapa_sexo = pd.read_csv(r'Y:\Base Rafael do Sergio\MAPA_NOME_SEXO.csv')

print(f"MEGA_BASE: {len(mega)} rows")
print(f"LINKAGEM: {len(link)} rows, unique nome_sabin: {link['nome_sabin'].nunique()}")

# Clean NOME in MEGA_BASE — use NOME_NORM if available
if 'NOME_NORM' in mega.columns:
    mega['nome_clean'] = mega['NOME_NORM'].str.strip().str.upper()
else:
    mega['nome_clean'] = mega['NOME'].str.strip().str.upper()
    # Remove "Id. Paciente : XXXX" suffix
    mega['nome_clean'] = mega['nome_clean'].str.replace(r'\s*ID\.?\s*PACIENTE\s*:?\s*\d+', '', regex=True).str.strip()

# Linked names (uppercase)
linked_names = set(link['nome_sabin'].str.strip().str.upper().unique())
print(f"Linked names: {len(linked_names)}")

# Parse DN_NORM for age
mega['DN_PARSED'] = pd.to_datetime(mega['DN_NORM'], errors='coerce', dayfirst=True)
# Also try DATA_NASCIMENTO
if mega['DN_PARSED'].isna().all():
    mega['DN_PARSED'] = pd.to_datetime(mega['DATA_NASCIMENTO'], errors='coerce', dayfirst=True)

# Parse DATA_ATENDIMENTO for visit dates
mega['DATA_ATEND_PARSED'] = pd.to_datetime(mega['DATA_ATENDIMENTO'], errors='coerce', dayfirst=True)
# Also try DATA_COLETA
mega['DATA_COLETA_PARSED'] = pd.to_datetime(mega['DATA_COLETA'].str[:10], errors='coerce', format='%Y-%m-%d')

# Use whichever date is available
mega['VISIT_DATE'] = mega['DATA_ATEND_PARSED'].fillna(mega['DATA_COLETA_PARSED'])

# Creatinine
mega['CREAT_NUM'] = pd.to_numeric(mega['CREATININA'], errors='coerce')

# Per-patient summary
patient_summary = mega.groupby('nome_clean').agg(
    DN=('DN_PARSED', 'first'),
    n_visits=('VISIT_DATE', 'count'),
    first_visit=('VISIT_DATE', 'min'),
    last_visit=('VISIT_DATE', 'max'),
    first_creat=('CREAT_NUM', 'first'),
).reset_index()

# Age at reference date (2024-01-01)
ref_date = pd.Timestamp('2024-01-01')
patient_summary['age'] = (ref_date - patient_summary['DN']).dt.days / 365.25
# Time span in days
patient_summary['time_span_days'] = (patient_summary['last_visit'] - patient_summary['first_visit']).dt.days

# Linkage status
patient_summary['linked'] = patient_summary['nome_clean'].isin(linked_names)

# Sex inference from first name
mapa_dict = dict(zip(mapa_sexo['first_name'].str.upper(), mapa_sexo['sexo']))
patient_summary['first_name'] = patient_summary['nome_clean'].str.split().str[0]
patient_summary['sexo_inferred'] = patient_summary['first_name'].map(mapa_dict)

print(f"\nUnique patients in MEGA_BASE: {len(patient_summary)}")
print(f"  Linked: {patient_summary['linked'].sum()}")
print(f"  Unlinked: {(~patient_summary['linked']).sum()}")

gL = patient_summary[patient_summary['linked'] == True]
gU = patient_summary[patient_summary['linked'] == False]

rows11 = []

# Continuous comparisons
for label, col in [('Age', 'age'), ('First_Creat', 'first_creat'),
                    ('N_visits', 'n_visits'), ('Time_span_days', 'time_span_days')]:
    lv = gL[col].dropna()
    uv = gU[col].dropna()
    smd = compute_smd(lv, uv)
    if len(lv) > 1 and len(uv) > 1:
        t, p = ttest_ind(lv, uv, equal_var=False)
    else:
        t, p = np.nan, np.nan

    rows11.append({
        'Variable': label,
        'Linked_N': len(lv),
        'Linked_Mean': lv.mean(),
        'Linked_SD': lv.std(),
        'Unlinked_N': len(uv),
        'Unlinked_Mean': uv.mean(),
        'Unlinked_SD': uv.std(),
        'SMD': smd,
        'T_stat': t,
        'P_value': p,
        'Imbalance_flag': 'YES' if abs(smd) > 0.20 else '',
    })

# Sex (categorical)
l_fem = (gL['sexo_inferred'] == 'Feminino').sum()
l_n = gL['sexo_inferred'].notna().sum()
u_fem = (gU['sexo_inferred'] == 'Feminino').sum()
u_n = gU['sexo_inferred'].notna().sum()

l_pct = l_fem / l_n * 100 if l_n > 0 else np.nan
u_pct = u_fem / u_n * 100 if u_n > 0 else np.nan

p_l = l_fem / l_n if l_n > 0 else 0
p_u = u_fem / u_n if u_n > 0 else 0
pp = (l_fem + u_fem) / (l_n + u_n) if (l_n + u_n) > 0 else 0
denom_sex = np.sqrt(pp * (1 - pp)) if 0 < pp < 1 else 1
smd_sex = (p_l - p_u) / denom_sex

try:
    chi2, p_chi, _, _ = chi2_contingency([[l_fem, l_n - l_fem], [u_fem, u_n - u_fem]])
except:
    chi2, p_chi = np.nan, np.nan

rows11.append({
    'Variable': 'Sex_Fem_%',
    'Linked_N': l_n,
    'Linked_Mean': l_pct,
    'Linked_SD': np.nan,
    'Unlinked_N': u_n,
    'Unlinked_Mean': u_pct,
    'Unlinked_SD': np.nan,
    'SMD': smd_sex,
    'T_stat': chi2,
    'P_value': p_chi,
    'Imbalance_flag': 'YES' if abs(smd_sex) > 0.20 else '',
})

tab11 = pd.DataFrame(rows11)
tab11.to_csv(os.path.join(OUT, 'Tab11_linked_vs_unlinked.csv'), index=False)

print("\n" + tab11.to_string(index=False, float_format='%.3f'))

flagged11 = tab11[tab11['Imbalance_flag'] == 'YES']
if len(flagged11) > 0:
    print(f"\n⚠ IMBALANCES DETECTED (SMD > 0.20):")
    for _, r in flagged11.iterrows():
        print(f"  {r['Variable']}: SMD = {r['SMD']:.3f}")
else:
    print("\n✓ No meaningful imbalances detected (all SMD ≤ 0.20)")


# ═══════════════════════════════════════════════════════════════
# INTERPRETATION SUMMARY
# ═══════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("INTERPRETATION NOTES")
print("="*80)

print("""
ANALYSIS 1 (Dose-Response):
- If Spearman rho is significantly positive for TWL vs Divergence, greater weight loss
  causes more indexed-to-nonindexed divergence (BSA confounding is dose-dependent).
- ΔBSA quartiles should show even clearer relationship since BSA change is the
  direct mediator.

ANALYSIS 2 (Mosteller):
- If DuBois and Mosteller non-indexed deltas are nearly identical (small paired-t),
  the BSA formula choice is robust — findings don't depend on which BSA formula.
- Shapley decomposition quantifies how much of the indexed Δ is "true GFR change"
  vs "BSA artifact".

ANALYSIS 3 (Representativeness):
- SMD < 0.10 = negligible difference, 0.10-0.20 = small, > 0.20 = potentially meaningful.
- If no imbalances detected, the observed-weight subset generalizes well to all
  patients with paired indexed eGFR.

ANALYSIS 4 (Linked vs Unlinked):
- Assesses selection bias in the Galvão linkage process.
- If linked patients differ substantially from unlinked, linkage bias should be
  acknowledged as a limitation.
""")

print("\n✓ ALL ANALYSES COMPLETE")
print(f"Output files saved to: {OUT}")
for f in sorted(os.listdir(OUT)):
    if f.startswith('Tab'):
        fp = os.path.join(OUT, f)
        print(f"  {f} ({os.path.getsize(fp)} bytes)")
