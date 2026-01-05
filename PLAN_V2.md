# JANUS-S V2 Analysis Plan

## Status: In Progress
**Last updated:** 2026-01-05

---

## Results Summary (V2 Analysis)

| Dataset | N SNe | JANUS q₀ | χ²/dof | ΛCDM Ωₘ | χ²/dof | ΔAIC |
|---------|-------|----------|--------|---------|--------|------|
| Pantheon+ | 1701 | -0.0126 ± 0.0103 | 1.046 | 0.362 ± 0.019 | 1.031 | +25.2 |
| DES-SN5YR | 1820 | -0.0278 ± 0.0130 | 0.915 | 0.330 ± 0.015 | 0.898 | +31.0 |
| Combined | 3182 | -0.0958 ± 0.0064 | 1.091 | 0.256 ± 0.008 | 1.030 | +195.0 |

---

## Next Steps

### 1. Analyser les résultats [DONE]
- [x] Comprendre pourquoi q₀ varie entre datasets (-0.013 → -0.096)
- [x] Investiguer la tension dans le dataset Combined
- [x] Comparer q₀ combiné (-0.096) avec référence 2018 (-0.087)
- [x] Analyser les corrélations entre paramètres

**FINDINGS:**
- Tension de calibration ~0.1 mag entre Pantheon+ et DES-SN5YR
- q₀ Combined (-0.096) ≈ Ref 2018 (-0.087) ✓
- Couverture redshift très différente: Pantheon+ (z̄=0.16) vs DES (z̄=0.48)

### 2. Générer des figures publication [DONE]
- [x] Diagramme de Hubble avec résidus (JANUS vs ΛCDM) → fig1_hubble_v2.pdf
- [x] Comparaison paramètres (q₀, Ωₘ, ΔAIC) → fig2_parameters_v2.pdf
- [x] Qualité d'ajustement χ²/dof → fig3_chi2_v2.pdf
- [x] Tableau récapitulatif → fig4_summary_v2.pdf

### 3. Diagnostics MCMC [DONE]
- [x] Vérifier la convergence (autocorrélation) → tau=17-56, ESS>1100
- [x] Calculer Gelman-Rubin → R-hat < 1.03 (tous OK)
- [x] Tracer les chaînes → mcmc_trace_plots.pdf
- [x] Augmenter n_steps → Non nécessaire, convergence confirmée

**RÉSULTATS CONVERGENCE:**
| Modèle/Dataset | R-hat max | ESS min | Status |
|----------------|-----------|---------|--------|
| janus/Combined | 1.018 | 1960 | ✓ OK |
| janus/DES-SN5YR | 1.014 | 1934 | ✓ OK |
| janus/Pantheon+ | 1.026 | 1242 | ✓ OK |
| lcdm/Combined | 1.026 | 1694 | ✓ OK |
| lcdm/DES-SN5YR | 1.019 | 2360 | ✓ OK |
| lcdm/Pantheon+ | 1.017 | 1142 | ✓ OK |

### 4. Analyses complémentaires [DONE]
- [x] Tester différentes valeurs de H₀ (67, 70, 73) → h0_sensitivity.pdf
- [x] Résultat: q₀ et Ωₘ sont INDÉPENDANTS de H₀ (offset absorbe la dépendance)
- [ ] Explorer contraintes BAO/CMB si pertinent (optionnel)

**RÉSULTAT CLÉ - INDÉPENDANCE DE H₀:**
Les paramètres q₀ (JANUS) et Ωₘ (ΛCDM) sont totalement indépendants de H₀
car le paramètre offset marginalise sur la constante de Hubble.
→ Robustesse des résultats V2 confirmée

### 5. Rédaction article [DONE]
- [x] Mettre à jour l'article avec résultats V2 → article_v2.tex
- [x] Tableau comparatif JANUS vs ΛCDM (Table 1)
- [x] Discussion des implications cosmologiques
- [x] Compilation PDF réussie (8 pages)

**ARTICLE V2 CRÉÉ:**
`publications/article/article_v2.tex` (compilé: article_v2.pdf)
- Abstract mis à jour avec résultats V2
- Nouveau tableau résultats (Pantheon+, DES-SN5YR, Combined)
- Diagnostics MCMC documentés
- Discussion calibration tension
- Figures V2 intégrées

---

## Files

- Results: `results/v2_analysis/v2_results.json`
- Checkpoints: `results/v2_analysis/checkpoints/`
- Figures: `results/v2_analysis/figures/`
- Analysis script: `code/run_v2_checkpointed.py`
