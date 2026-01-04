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

## PLAN DE PUBLICATION

**Objectif :** Publication comparative reproduction 2018 + validation Pantheon+

**Phases :**
1. Préparation (datasets, environnement)
2. Reproduction article 2018 (JLA, 740 SNe)
3. Extension Pantheon+ (1550 SNe)
4. Comparaison JANUS vs Lambda-CDM (AIC, BIC)
5. Analyse MCMC (erreurs bayésiennes)
6. Rédaction article LaTeX
7. Production PDF final

**Fichiers créés :**
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

