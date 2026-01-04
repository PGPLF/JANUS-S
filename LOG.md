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

---

## PLAN DE PUBLICATION V0 - SUIVI

**Objectif :** Publication comparative reproduction 2018 + validation Pantheon+
**Livrable :** `publications/article/article.pdf`

---

### PHASE 1 : DONNEES
| Tâche | Statut | Notes |
|-------|--------|-------|
| Télécharger JLA (740 SNe) | ⬜ À faire | supernovae.in2p3.fr |
| Télécharger Pantheon+ (1550 SNe) | ⬜ À faire | github.com/PantheonPlusSH0ES |
| Extraire et organiser les fichiers | ⬜ À faire | data/jla/, data/pantheon/ |
| Valider chargement avec data_loader.py | ⬜ À faire | Test unitaire |

### PHASE 2 : REPRODUCTION 2018
| Tâche | Statut | Notes |
|-------|--------|-------|
| Charger données JLA | ⬜ À faire | 740 SNe Ia |
| Calculer modules de distance (mu) | ⬜ À faire | alpha=0.141, beta=3.101 |
| Ajuster modèle JANUS | ⬜ À faire | Objectif: q0 = -0.087 |
| Valider chi2/dof ≈ 0.89 | ⬜ À faire | Comparaison article 2018 |
| Générer diagramme de Hubble | ⬜ À faire | results/figures/ |
| Générer graphe des résidus | ⬜ À faire | results/figures/ |

### PHASE 3 : EXTENSION PANTHEON+
| Tâche | Statut | Notes |
|-------|--------|-------|
| Charger données Pantheon+ | ⬜ À faire | 1550 SNe Ia |
| Ajuster modèle JANUS | ⬜ À faire | Nouveau q0 |
| Tests robustesse (bootstrap) | ⬜ À faire | 1000 échantillons |
| Tests sous-échantillons (z<0.5, z>0.5) | ⬜ À faire | Cohérence |

### PHASE 4 : COMPARAISON JANUS vs LAMBDA-CDM
| Tâche | Statut | Notes |
|-------|--------|-------|
| Ajuster Lambda-CDM sur JLA | ⬜ À faire | Om=0.3, OL=0.7 |
| Ajuster Lambda-CDM sur Pantheon+ | ⬜ À faire | |
| Calculer AIC, BIC | ⬜ À faire | Critères de sélection |
| Calculer Delta chi2 | ⬜ À faire | |
| Interpréter préférence statistique | ⬜ À faire | |

### PHASE 5 : FIGURES & TABLES
| Tâche | Statut | Notes |
|-------|--------|-------|
| Fig. 1 : Hubble JLA + JANUS | ⬜ À faire | |
| Fig. 2 : Hubble JLA + LCDM | ⬜ À faire | |
| Fig. 3 : Hubble Pantheon+ + JANUS | ⬜ À faire | |
| Fig. 4 : Résidus comparatifs | ⬜ À faire | |
| Fig. 5 : Contours MCMC (optionnel) | ⬜ À faire | |
| Tab. 1 : Résultats JLA | ⬜ À faire | |
| Tab. 2 : Résultats Pantheon+ | ⬜ À faire | |
| Tab. 3 : Comparaison statistique | ⬜ À faire | |

### PHASE 6 : REDACTION LATEX
| Tâche | Statut | Notes |
|-------|--------|-------|
| Abstract | ⬜ À faire | 150-200 mots |
| Introduction | ⬜ À faire | Contexte, objectifs |
| Méthodes | ⬜ À faire | Données, modèles, ajustement |
| Résultats | ⬜ À faire | Tableaux, figures |
| Discussion | ⬜ À faire | Interprétation, limites |
| Conclusion | ⬜ À faire | |
| Bibliographie (BibTeX) | ⬜ À faire | references.bib |

### PHASE 7 : PRODUCTION PDF
| Tâche | Statut | Notes |
|-------|--------|-------|
| Compilation LaTeX | ⬜ À faire | pdflatex + bibtex |
| Relecture | ⬜ À faire | |
| Version finale V0 | ⬜ À faire | article.pdf |

---

### PROGRESSION GLOBALE

| Phase | Avancement |
|-------|------------|
| 1. Données | 0% |
| 2. Reproduction 2018 | 0% |
| 3. Extension Pantheon+ | 0% |
| 4. Comparaison | 0% |
| 5. Figures & Tables | 0% |
| 6. Rédaction | 0% |
| 7. PDF Final | 0% |
| **TOTAL** | **0%** |

---

**Fichiers du projet :**
- `publications/PLAN_PUBLICATION_COMPARATIVE.md`
- `code/janus_model.py`
- `code/data_loader.py`
- `code/fitting.py`

---

## PUBLICATIONS IDENTIFIÉES

### Publication principale sur les supernovae
- **D'Agostini, Petit (2018)** : "Constraints on Janus Cosmological model from recent observations of supernovae type Ia" - Astrophysics and Space Science, 363(7): 139

### Autres publications pertinentes
- Petit, Margnat, Zejli (2024) : "A bimetric cosmological model based on Andrei Sakharov's twin universe approach" - EPJC
- Petit, D'Agostini (2014) : "Can negative mass be considered in General Relativity?" - arXiv:1408.2451
- Petit, D'Agostini, Debergh (2019) : "Physical and Mathematical Consistency of the JCM" - Progress in Physics

---

## NOTES

Ce projet vise à reproduire et compléter les travaux de Jean-Pierre Petit sur l'application du modèle cosmologique JANUS aux observations des supernovae de type Ia, en vue de publications scientifiques.

**Dépôt GitHub :** https://github.com/PGPLF/JANUS-S

