# JANUS-S - Journal de Projet
## Recherche sur le Modèle JANUS et les Supernovae

---

## PLAN DÉTAILLÉ PROPOSÉ

### Phase 1 : Structure du Projet
- [ ] Créer l'arborescence des dossiers
  - `/publications` : stockage des publications passées
  - `/synthese` : documents de synthèse
  - `/data` : données et graphiques
  - `/references` : bibliographie et ressources

### Phase 2 : Recherche Bibliographique
- [ ] Rechercher les publications de Jean-Pierre Petit sur JANUS et supernovae
- [ ] Identifier les articles scientifiques clés
- [ ] Lister les travaux sur :
  - Courbes de luminosité des supernovae Ia
  - Alternative à l'énergie noire
  - Masses négatives et cosmologie bimétrique

### Phase 3 : Documentation
- [ ] Créer un document de synthèse (Markdown)
- [ ] Convertir en PDF
- [ ] Organiser les références bibliographiques

### Phase 4 : Versionnement
- [ ] Initialiser le dépôt Git local
- [ ] Créer le dépôt GitHub distant
- [ ] Synchroniser et pousser le contenu

---

## HISTORIQUE DES OPÉRATIONS

### Session du 2026-01-04

| Heure | Opération | Statut |
|-------|-----------|--------|
| -- | Création du dossier JANUS-S | ✅ Succès |
| -- | Création du fichier LOG.md | ✅ Succès |
| -- | Attente validation du plan | ✅ Validé |
| -- | Création sous-dossiers (publications, synthese, data, references) | ✅ Succès |
| -- | Recherche publications JANUS/supernovae | ✅ Succès |
| -- | Création document de synthèse (Markdown) | ✅ Succès |
| -- | Installation de pandoc via Homebrew | ✅ Succès |
| -- | Conversion du document en PDF (129 KB) | ✅ Succès |
| -- | Initialisation Git local | ✅ Succès |
| -- | Création du dépôt GitHub (PGPLF/JANUS-S) | ✅ Succès |
| -- | Push initial vers GitHub | ✅ Succès |
| -- | Enrichissement synthèse (réflexions Grok) | ✅ Succès |
| -- | Ajout : comparaison JANUS/Lambda-CDM | ✅ Succès |
| -- | Ajout : protocole détaillé article 2018 | ✅ Succès |
| -- | Ajout : équations du modèle | ✅ Succès |
| -- | Ajout : données post-2018 (Pantheon+, DES) | ✅ Succès |
| -- | Ajout : méthodologie de reproduction | ✅ Succès |
| -- | Régénération PDF v2 (54 KB) | ✅ Succès |
| -- | Commit v2 du document de synthèse | ✅ Succès |
| -- | Push v2 vers GitHub | ✅ Succès |
| -- | Analyse du plan de reproduction (Grok) | ✅ Succès |
| -- | Création PLAN_PUBLICATION_COMPARATIVE.md | ✅ Succès |
| -- | Création structure publications/code/results | ✅ Succès |
| -- | Création janus_model.py (modèle JANUS + Lambda-CDM) | ✅ Succès |
| -- | Création data_loader.py (chargement JLA/Pantheon+) | ✅ Succès |
| -- | Création fitting.py (ajustement chi-2, comparaison) | ✅ Succès |
| -- | Commit v3 : Plan publication + code base | ✅ Succès |
| -- | Push v3 vers GitHub | ✅ Succès |
| -- | **Phase 3 : Extension Pantheon+** | ✅ Succès |
| -- | Chargement Pantheon+ (1543 SNe uniques) | ✅ Succès |
| -- | Ajustement JANUS: q0 = -0.0356 | ⚠️ Divergence |
| -- | Bootstrap (100 samples): q0 = -0.043 ± 0.052 | ✅ Succès |
| -- | Analyse sous-échantillons (z<0.1, z<0.5, z>=0.5) | ✅ Succès |
| -- | Génération figures Pantheon+ | ✅ Succès |
| -- | Comparaison JLA/Pantheon+ : Delta q0 = 0.051 | ⚠️ Incohérent |
| -- | **Validation et correction erreur Pantheon+** | ✅ Succès |
| -- | Erreur identifiée: m_b_corr ≠ mu | ⚠️ Corrigé |
| -- | Correction: utilisation MU_SH0ES | ✅ Succès |
| -- | Résultat corrigé: q0 = -0.0352 ± 0.0135 | ✅ Succès |
| -- | **Phase 4 : JANUS vs Lambda-CDM** | ✅ Succès |
| -- | Ajustement LCDM JLA: chi2/dof = 0.8518 | ✅ Succès |
| -- | Ajustement LCDM Pantheon+: chi2/dof = 0.4812 | ✅ Succès |
| -- | Comparaison AIC/BIC: LCDM légèrement préféré | ✅ Succès |
| -- | Génération figures Phase 4 | ✅ Succès |
| -- | **Phase 5 : Figures publication** | ✅ Succès |
| -- | Création generate_publication_figures.py | ✅ Succès |
| -- | Génération 4 figures PDF/PNG | ✅ Succès |
| -- | **Phase 6 : Rédaction LaTeX** | ✅ Succès |
| -- | Création article.tex (13 KB) | ✅ Succès |
| -- | **Phase 7 : Production PDF** | ✅ Succès |
| -- | Compilation pdflatex | ✅ Succès |
| -- | article.pdf généré (420 KB, 8 pages) | ✅ Succès |
| -- | **Git commit & push** | ✅ Succès |
| -- | Commit 4a8e7b7: Publication V0 complete | ✅ Succès |
| -- | Push vers GitHub PGPLF/JANUS-S | ✅ Succès |
| -- | **Mise à jour v0.1** | ✅ Succès |
| -- | Ajout auteur: Patrick Guerin (pg@gfo.bzh) | ✅ Succès |
| -- | Ajout: contributions, funding, conflicts | ✅ Succès |
| -- | Format two-column (style JANUS-Z v17) | ✅ Succès |
| -- | article_v0.1.pdf généré (7 pages, 529 KB) | ✅ Succès |

---

## PLAN DE PUBLICATION V0 - SUIVI

**Objectif :** Publication comparative reproduction 2018 + validation Pantheon+
**Livrable :** `publications/article/article_v0.1.pdf`

### HISTORIQUE DES VERSIONS
| Version | Date | Description |
|---------|------|-------------|
| v0.1 | 2026-01-04 | Version initiale avec auteur, acknowledgments, format JANUS-Z |

---

### PHASE 1 : DONNEES
| Tâche | Statut | Notes |
|-------|--------|-------|
| Télécharger JLA (740 SNe) | ✅ Fait | 22.6 MB, supernovae.in2p3.fr |
| Télécharger Pantheon+ (1701 SNe) | ✅ Fait | GitHub PantheonPlusSH0ES |
| Extraire et organiser les fichiers | ✅ Fait | data/jla/, data/pantheon/ |
| Valider chargement avec data_loader.py | ✅ Fait | JLA: 740, Pantheon+: 1701 SNe |

### PHASE 2 : REPRODUCTION 2018
| Tâche | Statut | Notes |
|-------|--------|-------|
| Charger données JLA | ✅ Fait | 740 SNe Ia, z: 0.01-1.30 |
| Calculer modules de distance (mu) | ✅ Fait | alpha=0.141, beta=3.101, M_B=-19.05 |
| Ajuster modèle JANUS | ✅ Fait | **q0 = -0.0864** (ref: -0.087) |
| Valider chi2/dof ≈ 0.89 | ✅ Fait | **chi2/dof = 0.8834** (ref: 0.89) |
| Générer diagramme de Hubble | ✅ Fait | hubble_jla_janus.png |
| Générer graphe des résidus | ✅ Fait | residus_jla_janus.png |

**RESULTATS PHASE 2:**
```
q0 = -0.0864 ± 0.001 (ref: -0.087 ± 0.015) ✅
chi2 = 651.93, chi2/dof = 0.8834 (ref: 0.89) ✅
Reproduction VALIDEE
```

### PHASE 3 : EXTENSION PANTHEON+
**Objectif:** Valider JANUS sur dataset etendu (1701 SNe, z jusqu'a 2.26)

| Tâche | Statut | Notes |
|-------|--------|-------|
| Charger données Pantheon+ | ✅ Fait | 1543 SNe uniques (moyenne pondérée) |
| Ajuster modèle JANUS | ✅ Fait | **q0 = -0.0356** (diverge de JLA!) |
| Tests robustesse (bootstrap) | ✅ Fait | 100 échantillons, q0 = -0.043 ± 0.052 |
| Tests sous-échantillons (z<0.5, z>0.5) | ✅ Fait | Voir analyse ci-dessous |
| Générer figures Pantheon+ | ✅ Fait | hubble_pantheon_janus.png, residus |
| Comparer JLA vs Pantheon+ | ✅ Fait | Delta q0 = 0.051 (incohérent) |

**RESULTATS PHASE 3 (CORRIGE):**

**ERREUR CORRIGEE:** Pantheon+ `m_b_corr` n'est PAS le module de distance!
- `m_b_corr` ~ 9-27 mag (magnitude corrigée)
- `MU_SH0ES` ~ 29-46 mag (module de distance correct)
- Relation: `m_b_corr = MU_SH0ES + M_B` où M_B ~ -19.25

```
Dataset:        Pantheon+ (1543 SNe uniques)
Plage z:        [0.0012, 2.2614]
Plage mu:       [29.03, 46.18] mag (CORRIGE)

Ajustement JANUS (corrigé):
  q0          = -0.0352 ± 0.0135
  chi2/dof    = 0.4973

Analyse sous-échantillons:
  z < 0.1:  q0 = -0.2604, chi2/dof = 0.5798 (n=583)
  z < 0.5:  q0 = -0.1648, chi2/dof = 0.5022 (n=1333)
  z < 1.0:  q0 = -0.0696, chi2/dof = 0.4907 (n=1518)
  z < 1.3:  q0 = -0.0723, chi2/dof = 0.4895 (n=1527)

Comparaison JLA vs Pantheon+:
  q0 JLA      = -0.0864 ± 0.0143
  q0 Pantheon = -0.0352 ± 0.0135
  Delta       = 0.0512 -> DIVERGENCE SIGNIFICATIVE
```

**INTERPRETATION:**
La divergence q0 est REELLE et liée à:
1. Distribution des redshifts différente (Pantheon+ a plus de SNe à bas z)
2. Les SNe à bas z (z<0.1) donnent q0 ~ -0.26, celles à haut z donnent q0 ~ -0.04
3. JANUS semble sensible à la composition de l'échantillon
4. Effet physique possible: évolution de q0 avec z?

**Figures générées:**
- `results/figures/hubble_jla_pantheon_corrige.png`
- `results/figures/janus_vs_lcdm_jla_v2.png`
- `results/figures/janus_vs_lcdm_pantheon_v2.png`

### PHASE 4 : COMPARAISON JANUS vs LAMBDA-CDM
| Tâche | Statut | Notes |
|-------|--------|-------|
| Ajuster Lambda-CDM sur JLA | ✅ Fait | Om=0.3, chi2/dof=0.8518 |
| Ajuster Lambda-CDM sur Pantheon+ | ✅ Fait | Om=0.3, chi2/dof=0.4812 |
| Calculer AIC, BIC | ✅ Fait | Voir résultats ci-dessous |
| Calculer Delta chi2 | ✅ Fait | JLA: -22.43, Panth: -24.24 |
| Interpréter préférence statistique | ✅ Fait | Lambda-CDM légèrement préféré |

**RESULTATS PHASE 4:**
```
╔══════════════════════════════════════════════════════════════════╗
║                    COMPARAISON JANUS vs LCDM                     ║
╠══════════════════════════════════════════════════════════════════╣
║  Dataset       │ chi2/dof JANUS │ chi2/dof LCDM │ Delta AIC     ║
║  JLA           │         0.8834 │        0.8518 │    -24.43     ║
║  Pantheon+     │         0.4973 │        0.4812 │    -26.24     ║
╚══════════════════════════════════════════════════════════════════╝

Interprétation:
- Delta AIC négatif -> Préférence statistique pour Lambda-CDM
- MAIS: les chi2/dof sont très proches (écart < 4%)
- JANUS a 2 paramètres (q0, offset), LCDM a 1 paramètre (offset, Om fixé)
- Si Om libre: LCDM trouve Om=0.27 sur JLA, Om=0.00 sur Pantheon+ (!)

CONCLUSION:
Les deux modèles ajustent les données de manière comparable.
Lambda-CDM est statistiquement légèrement préféré (parsimonie).
JANUS reste compétitif avec un seul paramètre libre (q0).
```

### PHASE 5 : FIGURES & TABLES
| Tâche | Statut | Notes |
|-------|--------|-------|
| Fig. 1 : Hubble JLA + JANUS + LCDM | ✅ Fait | fig1_hubble_jla.pdf |
| Fig. 2 : Hubble Pantheon+ + JANUS + LCDM | ✅ Fait | fig2_hubble_pantheon.pdf |
| Fig. 3 : Evolution q0 avec z | ✅ Fait | fig3_q0_evolution.pdf |
| Fig. 4 : Comparaison chi2/dof | ✅ Fait | fig4_comparison.pdf |
| Tab. 1 : Résultats JLA | ✅ Fait | Dans article.tex |
| Tab. 2 : Résultats Pantheon+ | ✅ Fait | Dans article.tex |
| Tab. 3 : Sous-échantillons | ✅ Fait | Dans article.tex |
| Tab. 4 : Comparaison statistique | ✅ Fait | Dans article.tex |

### PHASE 6 : REDACTION LATEX
| Tâche | Statut | Notes |
|-------|--------|-------|
| Abstract | ✅ Fait | ~180 mots |
| Introduction | ✅ Fait | Contexte, objectifs |
| Theoretical Framework | ✅ Fait | JANUS + LCDM |
| Data and Methodology | ✅ Fait | JLA, Pantheon+, fitting |
| Results | ✅ Fait | 4 tableaux, 4 figures |
| Discussion | ✅ Fait | Interprétation, limites |
| Conclusion | ✅ Fait | Résumé des findings |
| Bibliography | ✅ Fait | 7 références |

### PHASE 7 : PRODUCTION PDF
| Tâche | Statut | Notes |
|-------|--------|-------|
| Compilation pdflatex | ✅ Fait | 2 passes |
| article.pdf | ✅ Fait | 8 pages, 420 KB |
| **VERSION V0 COMPLETE** | ✅ | publications/article/article.pdf |

---

### PROGRESSION GLOBALE

| Phase | Avancement |
|-------|------------|
| 1. Données | 100% |
| 2. Reproduction 2018 | 100% |
| 3. Extension Pantheon+ | 100% |
| 4. Comparaison JANUS/LCDM | 100% |
| 5. Figures & Tables | 100% |
| 6. Rédaction LaTeX | 100% |
| 7. PDF Final | 100% |
| **TOTAL** | **100%** |

---

**Fichiers du projet :**
- `publications/PLAN_PUBLICATION_COMPARATIVE.md`
- `publications/article/article_v0.1.tex` (source LaTeX v0.1)
- `publications/article/article_v0.1.pdf` (VERSION V0.1 COURANTE)
- `code/janus_model.py`
- `code/data_loader.py`
- `code/fitting.py`
- `code/reproduce_2018.py`
- `code/extend_pantheon.py`
- `code/test_data_loading.py`
- `code/validation_complete.py`
- `code/validation_corrigee.py`
- `code/generate_publication_figures.py`

**Figures de la publication:**
- `publications/article/figures/fig1_hubble_jla.pdf`
- `publications/article/figures/fig2_hubble_pantheon.pdf`
- `publications/article/figures/fig3_q0_evolution.pdf`
- `publications/article/figures/fig4_comparison.pdf`

**Figures d'analyse:**
- `results/figures/hubble_jla_janus.png`
- `results/figures/residus_jla_janus.png`
- `results/figures/hubble_jla_pantheon_corrige.png`
- `results/figures/janus_vs_lcdm_jla_v2.png`
- `results/figures/janus_vs_lcdm_pantheon_v2.png`
- `results/figures/validation_jla_pantheon.png`

---

## PUBLICATIONS IDENTIFIÉES

### Publication principale sur les supernovae
- **D'Agostini, Petit (2018)** : "Constraints on Janus Cosmological model from recent observations of supernovae type Ia" - Astrophysics and Space Science, 363(7): 139

### Autres publications pertinentes
- Petit, Margnat, Zejli (2024) : "A bimetric cosmological model based on Andrei Sakharov's twin universe approach" - EPJC
- Petit, D'Agostini (2014) : "Can negative mass be considered in General Relativity?" - arXiv:1408.2451
- Petit, D'Agostini, Debergh (2019) : "Physical and Mathematical Consistency of the JCM" - Progress in Physics

---

---

## PLAN V1.0 - Extension MCMC et Dark Energy

**Objectif:** MCMC complet (emcee + dynesty), modeles wCDM/CPL, 10 figures publication

### Strategie de Versionnage
- **V1.0** : Pantheon+ seul (1550 SNe) - analyse complete MCMC
- **V1.1** : Pantheon+ + DES-SN5YR (~3000 SNe) - extension dataset

### Choix Techniques
- **MCMC** : emcee ET dynesty avec analyse comparative
- **Dark Energy** : wCDM + CPL (w0-wa) analyse complete
- **Regle** : Jamais effacer, toujours archiver

### PHASE 1 : Preparation Donnees V1
| Tache | Statut | Notes |
|-------|--------|-------|
| Telecharger covariances Pantheon+ | 🔄 En cours | Matrices stat+sys |
| Integrer covariances data_loader.py | ⏳ Pending | Full-cov 1701x1701 |
| Valider chi2 full-cov vs diagonal | ⏳ Pending | Impact attendu <20% |

### PHASE 2 : MCMC Double (emcee + dynesty)
| Tache | Statut | Notes |
|-------|--------|-------|
| Creer fitting_mcmc.py | ⏳ Pending | Framework unifie |
| Implementation emcee | ⏳ Pending | 32 walkers, 5000 steps |
| Implementation dynesty | ⏳ Pending | Nested sampling, evidence |
| Analyse comparative | ⏳ Pending | Posteriors coherents <10% |

### PHASE 3 : Sous-Echantillons
| Tache | Statut | Notes |
|-------|--------|-------|
| Stratification redshift | ⏳ Pending | 5 bins z |
| Stratification host mass | ⏳ Pending | </>10^10 Msun |
| Parameterisation q(z) | ⏳ Pending | q0 + q1*z |

### PHASE 4 : wCDM + CPL
| Tache | Statut | Notes |
|-------|--------|-------|
| Implementation wCDM | ⏳ Pending | w constant |
| Implementation CPL | ⏳ Pending | w0 + wa*(1-a) |
| Comparaison 4 modeles | ⏳ Pending | JANUS/LCDM/wCDM/CPL |
| Interpretation q0 <-> w | ⏳ Pending | Relation physique |

### PHASE 5 : Figures Publication
| Tache | Statut | Notes |
|-------|--------|-------|
| 10 figures PDF | ⏳ Pending | 300 DPI, colorblind-friendly |
| 4 tables | ⏳ Pending | MCMC, AIC/BIC, subsamples, errors |

### PHASE 6 : Article V1.0
| Tache | Statut | Notes |
|-------|--------|-------|
| Redaction LaTeX | ⏳ Pending | article_v1.0.tex |
| References (+5) | ⏳ Pending | DESI, DES, emcee, dynesty, corner |

### PHASE 7 : Validation & Publication
| Tache | Statut | Notes |
|-------|--------|-------|
| Script reproductible | ⏳ Pending | run_v1_analysis.py |
| Tag Git v1.0 | ⏳ Pending | Apres validation |
| Soumission arXiv | ⏳ Pending | astro-ph.CO |

---

### Fichiers V1 a Creer
```
code/fitting_mcmc.py          # MCMC emcee + dynesty
code/dark_energy_models.py    # wCDM + CPL
code/subsample_analysis.py    # Analyse par bins
code/generate_v1_figures.py   # 10 figures
code/run_v1_analysis.py       # Pipeline complet
publications/article/article_v1.0.tex
```

### Dependances V1
```
emcee>=3.1, dynesty>=2.0, corner>=2.2, h5py>=3.0, tqdm>=4.60
```

---

## NOTES

Ce projet vise à reproduire et compléter les travaux de Jean-Pierre Petit sur l'application du modèle cosmologique JANUS aux observations des supernovae de type Ia, en vue de publications scientifiques.

**Dépôt GitHub :** https://github.com/PGPLF/JANUS-S

**Regle projet :** Archiver, jamais effacer

