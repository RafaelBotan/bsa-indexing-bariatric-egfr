# BSA Indexing Artifact After Bariatric Surgery

**Body surface area indexing can reverse the apparent direction of creatinine-based estimated glomerular filtration rate change 12 months after bariatric surgery**

> Companion repository for the article submitted to *Nephrology Dialysis Transplantation* (NDT).

## Summary

This repository contains the analysis code and aggregated results for a dual-cohort study demonstrating that BSA indexing can reverse the apparent direction of creatinine-based eGFR change after Roux-en-Y gastric bypass (RYGB).

### Key Findings (12 months post-surgery)

| Metric | Discovery (Arruda) | Replication (Galvao) | Pooled |
|--------|:--:|:--:|:--:|
| N indexed | 268 | 283 | 551 |
| N nonindexed | 154 | 96 | 250 |
| Delta eGFR indexed | +3.60 | +2.80 | +3.19 |
| Delta eGFR nonindexed | -16.19 | -18.70 | -17.15 |
| Divergence (mL/min) | 19.22 | 21.21 | 19.99 |
| Shapley %BSA | 84.9% [74.8, 97.2] | 88.1% [82.1, 94.5] | 86.1% [79.6, 93.8] |
| Directional discordance | 43.5% | 57.3% | 48.8% [42.7, 55.0] |
| Reverse pattern | 0% | 0% | 0% |

## Study Design

- **Discovery cohort (Arruda):** 1,869 RYGB patients, CPF-linked to Sabin laboratory
- **Replication cohort (Galvao):** 10,872 patients, name-linked to Sabin laboratory (tiers 1+2: N=2,300)
- **Shared laboratory:** Sabin Laboratorios, Brasilia, Brazil
- **eGFR equation:** CKD-EPI 2021 race-free
- **BSA formula:** DuBois (primary), Mosteller (sensitivity)
- **Decomposition:** Shapley bilinear identity with 5,000-replicate bootstrap CIs

## Repository Structure

```
scripts/          Analysis scripts (Python 3.13)
  analise_estudo_c_v4.py         Main analysis pipeline
  analise_complementar_v4.py     Sensitivity analyses
  analise_blindagem_v4.py        Blinding and temporal validation
  bootstrap_shapley.py           Bootstrap CIs for Shapley components
  gerar_docx_completo.py         Manuscript DOCX generator

tables/           Aggregated results (no individual patient data)
  Tab1_baseline.csv              Baseline characteristics
  Tab2_paired_overall.csv        12M paired changes
  Tab3_divergence.csv            Indexed-nonindexed divergence
  Tab4_shapley.csv               Shapley decomposition
  Tab5_kdigo_transitions.csv     G-category transitions
  Tab6_forest_plot_data.csv      Forest plot data
  Tab7_sensitivity.csv           Subgroup sensitivities
  Tab8_dose_response_TWL.csv     Dose-response by %TWL
  Tab8b_dose_response_BSA.csv    Dose-response by delta BSA
  Tab9_mosteller_sensitivity.csv Mosteller vs DuBois
  Tab10_representativeness.csv   NI subset representativeness
  Tab11_linked_vs_unlinked.csv   Linked vs unlinked comparison
  Tab12_age_coverage.csv         Age/date coverage
  Tab12b_exam_dates.csv          Exam date distribution
  Tab13_directional_discordance.csv  Directional discordance
  Tab14_pre_sensitivity.csv      PRE window sensitivity
  Tab15_lme_results.csv          LME coefficients
  Tab16_galvao_date_validated.csv    Galvao date-validated sensitivity
  Tab17_shapley_bootstrap.csv    Bootstrap CIs for Shapley

figures/          Publication-ready figures (300 DPI)
  fig1_flowchart.png             STROBE/RECORD flow diagram
  fig2_shapley.png               Shapley decomposition
  fig3_forest_discordance.png    Forest plot + discordance
  fig4_dose_response.png         Dose-response
  sfig1-sfig7                    Supplementary figures
```

## Data Availability

Individual-level patient data cannot be shared publicly due to privacy regulations (LGPD — Brazilian General Data Protection Law). The dataset contains linked clinical and laboratory records from two surgical practices with identifiable information (names, dates of birth, national IDs).

Aggregated results are provided in the `tables/` directory. Requests for de-identified data access should be directed to the corresponding author.

## Reporting Guidelines

This study follows the **STROBE** statement and **RECORD** extension for studies using routinely collected health data with linkage.

## Requirements

- Python 3.13+
- pandas, numpy, scipy, statsmodels, matplotlib
- python-docx (for manuscript generation)

## License

Code: MIT License
Data/Results: CC-BY 4.0

## Citation

[To be updated upon publication]
