# Synthese : Modele Cosmologique JANUS et Supernovae

## Document de Reference pour le Projet JANUS-S
**Date de creation :** 4 janvier 2026
**Derniere mise a jour :** 4 janvier 2026 (v2 - enrichissement)
**Objectif :** Reproduire et completer les travaux de Jean-Pierre Petit sur les supernovae

---

## 1. Introduction au Modele JANUS

Le modele cosmologique JANUS (Janus Cosmological Model - JCM) est une theorie bimetrique developpee par Jean-Pierre Petit depuis 1977. Ce modele fusionne :
- La relativite generale d'Albert Einstein
- Les travaux d'Andrei Sakharov en physique des particules et cosmologie (univers jumeaux)
- La geometrie symplectique de Jean-Marie Souriau

### Principe fondamental
Le modele propose l'existence de deux secteurs couples :
- **Secteur positif** : masses positives (notre univers observable)
- **Secteur negatif** : masses negatives (secteur miroir, non observable directement)

Les interactions suivent des regles specifiques :
- Masses de meme signe : attraction selon la loi de Newton
- Masses de signes opposes : repulsion selon une "loi anti-Newton"

### Caracteristiques cles
- Univers vu comme une hypersurface a 4 dimensions avec topologie fermee
- Elimination des singularites
- Explication de l'asymetrie matiere-antimatiere
- Energie globale negative dominee par les masses negatives

---

## 2. Comparaison JANUS vs Lambda-CDM

### 2.1 Le modele standard Lambda-CDM

| Composante | Proportion | Role |
|------------|------------|------|
| Matiere ordinaire | ~5% | Matiere baryonique observable |
| Matiere noire (CDM) | ~27% | Explique rotation des galaxies, structures |
| Energie noire (Lambda) | ~68% | Explique l'acceleration de l'expansion |

**Interpretation des supernovae dans Lambda-CDM :**
- L'acceleration de l'expansion s'explique par l'energie noire
- Constante cosmologique Lambda exerce une pression repulsive
- Necessite des composantes hypothetiques non directement observees

### 2.2 Le modele JANUS

| Aspect | Modele JANUS |
|--------|--------------|
| Matiere noire | Remplacee par effets des masses negatives |
| Energie noire | Remplacee par repulsion masses +/- |
| Parametres libres | 1 seul (q_0) vs plusieurs pour Lambda-CDM |
| Constante cosmologique | Aucune |

**Interpretation des supernovae dans JANUS :**
- L'acceleration vient de la repulsion gravitationnelle entre masses positives et negatives
- Les masses negatives dominent energetiquement (energie globale negative)
- Meme ajustement aux donnees mais explication physique differente

### 2.3 Tableau comparatif

| Critere | Lambda-CDM | JANUS |
|---------|------------|-------|
| Composantes invisibles | Oui (matiere noire, energie noire) | Non |
| Parametres libres | Plusieurs (Omega_m, Omega_Lambda, w, ...) | 1 (q_0) |
| Ajustement SNe Ia | Excellent | Excellent (comparable ou superieur) |
| Constante cosmologique | Oui (ad hoc) | Non |
| Explication physique | Pression repulsive hypothetique | Antigravite naturelle |

---

## 3. Protocole de l'Article 2018

### 3.1 Reference complete

**Titre :** "Constraints on Janus Cosmological model from recent observations of supernovae type Ia"

| Champ | Information |
|-------|-------------|
| Auteurs | G. D'Agostini, J.-P. Petit |
| Journal | Astrophysics and Space Science |
| Volume | 363, Issue 7, Article 139 |
| Date | Juin 2018 |
| DOI | 10.1007/s10509-018-3365-3 |

### 3.2 Donnees utilisees

- **Source :** Compilation SDSS-II + SNLS (ensemble JLA - Joint Light-curve Analysis)
- **Reference :** Betoule et al. (2014)
- **Nombre :** 740 supernovae de type Ia
- **Couverture :** Redshifts de z proche a z ~ 1.3

**Parametres observes par supernova :**
- m_B* : magnitude de pic en bande B au repos
- X_1 : etirement temporel de la courbe de lumiere
- C : couleur au maximum de luminosite

### 3.3 Equations du modele

**Equations de champ bimétriques :**

Pour le secteur positif (+) :
```
R^(+)_μν - (1/2) R^(+) g^(+)_μν = χ [ T^(+)_μν + (a^(-)_3 / a^(+)_3) T^(-)_μν ]
```

Pour le secteur negatif (-) :
```
R^(-)_μν - (1/2) R^(-) g^(-)_μν = -χ [ T^(-)_μν + (a^(+)_3 / a^(-)_3) T^(+)_μν ]
```

ou R est le tenseur de Ricci, T les tenseurs energie-impulsion, χ une constante, a les facteurs d'echelle.

**Solution exacte (ere de poussiere) :**

Facteur d'echelle :
```
a_(+)(u) = (α²/2) cosh²(u)
```

Temps cosmique :
```
t_(+)(u) = (α²c/2) [1 + sinh(2u)/2 + u]
```

avec α² = -(8πG/3c²)E, ou E < 0 (densite energetique negative dominante).

**Distance comobile :**
```
r = (c / a_0 H_0) × [q_0 z + (1-q_0)(1 - √(1+2q_0 z))] / [q_0² (1+z)]
```

Analogue a la relation de Mattig, adaptee pour q_0 < 0.

**Module de distance :**
```
μ = m*_B - M_B + α X_1 - β C
```

ou M_B, α, β sont des parametres de nuisance fixes aux valeurs Lambda-CDM de Betoule et al.

### 3.4 Methode d'ajustement

1. Calcul du module de distance theorique pour chaque z
2. Ajustement aux 740 points de donnees
3. Minimisation du χ² avec q_0 et une constante comme parametres libres
4. Comparaison via diagrammes de Hubble (magnitude vs z)

### 3.5 Resultats obtenus

| Parametre | Valeur |
|-----------|--------|
| q_0 (parametre de deceleration) | -0.087 ± 0.015 |
| χ²/d.o.f. | 657/738 |
| Age de l'univers T_0 | ~15 Gyr (pour H_0 = 70 km/s/Mpc) |

**Interpretation :** q_0 < 0 indique une acceleration de l'expansion, expliquee par les masses negatives sans energie noire.

---

## 4. Publications Completes

### 4.1 Articles sur arXiv

| Annee | Titre | Reference |
|-------|-------|-----------|
| 2024 | A bimetric cosmological model based on Andrei Sakharov's twin universe approach | arXiv:2412.04644 |
| 2014 | Can negative mass be considered in General Relativity? | arXiv:1408.2451 |
| 2008 | Bigravity: a bimetric model of the Universe with variable constants, including VSL | arXiv:0803.1362 |
| 2007 | Bigravity as an interpretation of the cosmic acceleration | arXiv:0712.0067 |

### 4.2 Articles dans des Revues a Comite de Lecture

| Annee | Titre | Journal |
|-------|-------|---------|
| 2024 | Study of symmetries through the action on torsors of the Janus symplectic group | Reviews in Mathematical Physics |
| 2024 | A bimetric cosmological model based on Andrei Sakharov's twin universe approach | European Physical Journal C, 84: 1226 |
| 2019 | Physical and Mathematical Consistency of the Janus Cosmological Model (JCM) | Progress in Physics, 15(1): 38-47 |
| 2018 | On evidence for negative energies and masses in the Dirac equation | Journal of Physics Communications, 2(11): 115012 |
| 2018 | Constraints on Janus Cosmological model from recent observations of supernovae type Ia | Astrophysics and Space Science, 363(7): 139 |
| 2015 | Lagrangian derivation of the two coupled field equations in the Janus cosmological model | Astrophysics and Space Science |
| 2014 | Cosmological bimetric model with interacting positive and negative masses | Modern Physics Letters A, 29(34): 1450182 |
| 2014 | Negative mass hypothesis in cosmology and the nature of dark energy | Astrophysics and Space Science, 354(2): 611-615 |

### 4.3 Publications Fondatrices

| Annee | Titre | Journal |
|-------|-------|---------|
| 1995 | Twin Universe Cosmology | Astrophysics and Space Science, 226: 273-307 |
| 1994 | The missing mass problem | Il Nuovo Cimento, Vol.109: 697-710 |
| 1988 | Cosmological model with variable velocity of light | Modern Physics Letters A3: 1527 |

---

## 5. Donnees Disponibles Post-2018

### 5.1 Catalogues de supernovae recents

| Catalogue | Annee | Nombre SNe | Caracteristiques |
|-----------|-------|------------|------------------|
| JLA | 2014 | 740 | Utilise dans l'article 2018 |
| Pantheon | 2018 | ~1048 | Scolnic et al., z jusqu'a ~2 |
| Pantheon+ | 2022 | 1701 courbes / 1550 SNe | Brout et al., meilleure calibration |
| DES (5 ans) | 2024 | ~1500 | Resultats finaux janvier 2024 |

### 5.2 Sources de donnees

- **JLA :** Archives SNLS/SDSS, CDS Strasbourg
- **Pantheon/Pantheon+ :** GitHub du Supernova Cosmology Project
- **DES :** Dark Energy Survey Data Release
- **Autres :** NASA/IPAC, Union3 (2022), donnees JWST (haut z)

### 5.3 Interet pour le projet JANUS-S

Les nouvelles donnees permettraient de :
1. Tester la robustesse du modele sur plus de supernovae
2. Explorer des redshifts plus eleves (z > 2 avec JWST)
3. Comparer statistiquement JANUS vs Lambda-CDM (AIC/BIC)
4. Detecter d'eventuels ecarts a bas vs haut redshift

---

## 6. Methodologie pour Reproduire les Travaux

### 6.1 Etapes de reproduction

1. **Telecharger les donnees JLA** (740 SNe Ia)
   - Fichiers : m_B*, X_1, C, z, covariances
   - Sources : CDS, archives SNLS/SDSS

2. **Implementer les equations JANUS**
   - Distance comobile r(z, q_0)
   - Module de distance μ(z, q_0)

3. **Ajustement numerique**
   - Minimisation χ² avec q_0 comme parametre libre
   - Fixer M_B, α, β aux valeurs Lambda-CDM pour comparaison directe

4. **Validation**
   - Reproduire q_0 ≈ -0.087
   - Verifier χ²/d.o.f. ≈ 657/738

### 6.2 Outils recommandes

| Outil | Usage |
|-------|-------|
| Python + NumPy/SciPy | Calculs et minimisation χ² |
| emcee | Analyse MCMC si necessaire |
| Astropy | Manipulation des donnees astronomiques |
| sncosmo | Package specialise supernovae |
| Matplotlib | Diagrammes de Hubble, residus |

### 6.3 Extension avec donnees recentes

1. Telecharger Pantheon+ (disponible sur GitHub)
2. Adapter le code pour le nouveau format
3. Tenir compte des covariances systematiques ameliorees
4. Comparer les ajustements JANUS vs Lambda-CDM
5. Calculer criteres AIC/BIC pour comparaison de modeles

---

## 7. Pistes pour Completer les Travaux

### 7.1 Extensions scientifiques

- **Nouvelles donnees :** Pantheon+, DES, futures LSST/Euclid
- **Haut redshift :** Supernovae z > 2 (JWST)
- **Comparaison statistique :** Tests AIC/BIC, validation croisee
- **Autres observables :** CMB (Planck), BAO (DESI)

### 7.2 Contributions possibles

1. Premiere analyse JANUS sur Pantheon+ (non publiee a ce jour)
2. Comparaison rigoureuse des ajustements avec Lambda-CDM
3. Etude des tensions cosmologiques (H_0) dans le cadre JANUS
4. Predictions testables pour futures observations

### 7.3 Points a approfondir

- Derivation mathematique complete des equations de champ
- Implementation numerique optimisee
- Analyse des incertitudes systematiques
- Coherence avec CMB et structures a grande echelle

---

## 8. Ressources et Liens

### Sites Web Officiels
- Site JANUS : https://januscosmologicalmodel.com/
- Page Jean-Pierre Petit : http://www.jp-petit.org/

### Bases de Donnees Scientifiques
- NASA ADS : https://ui.adsabs.harvard.edu/
- arXiv : https://arxiv.org/
- HAL : https://hal.science/
- CDS Strasbourg : https://cds.u-strasbg.fr/

### Donnees Supernovae
- Pantheon+ GitHub : https://github.com/PantheonPlusSH0ES
- JLA : via CDS ou archives SNLS
- DES : https://www.darkenergysurvey.org/

### Profils Chercheurs
- ResearchGate : https://www.researchgate.net/profile/Jean-Pierre-Petit
- Academia : https://www.academia.edu/

---

## 9. Bibliographie

1. D'Agostini, G.; Petit, J.-P. (2018). "Constraints on Janus Cosmological model from recent observations of supernovae type Ia". *Astrophysics and Space Science*, 363(7): 139.

2. Petit, J.-P.; Margnat, F.; Zejli, H. (2024). "A bimetric cosmological model based on Andrei Sakharov's twin universe approach". *European Physical Journal C*, 84: 1226.

3. Petit, J.-P.; D'Agostini, G. (2014). "Can negative mass be considered in General Relativity?". arXiv:1408.2451.

4. Petit, J.-P.; D'Agostini, G.; Debergh, N. (2019). "Physical and Mathematical Consistency of the Janus Cosmological Model (JCM)". *Progress in Physics*, 15(1): 38-47.

5. Petit, J.-P. (1995). "Twin Universe Cosmology". *Astrophysics and Space Science*, 226: 273-307.

6. Debergh, N.; Petit, J.-P.; D'Agostini, G. (2018). "On evidence for negative energies and masses in the Dirac equation through a unitary time-reversal operator". *Journal of Physics Communications*, 2(11): 115012.

7. Betoule, M. et al. (2014). "Improved cosmological constraints from a joint analysis of the SDSS-II and SNLS supernova samples". *Astronomy & Astrophysics*, 568: A22.

8. Scolnic, D. et al. (2018). "The Complete Light-curve Sample of Spectroscopically Confirmed SNe Ia from Pan-STARRS1 and Cosmological Constraints from the Combined Pantheon Sample". *The Astrophysical Journal*, 859(2): 101.

9. Brout, D. et al. (2022). "The Pantheon+ Analysis: Cosmological Constraints". *The Astrophysical Journal*, 938(2): 110.

---

## Historique des versions

| Version | Date | Modifications |
|---------|------|---------------|
| v1 | 04/01/2026 | Creation initiale |
| v2 | 04/01/2026 | Ajout comparaison JANUS/Lambda-CDM, protocole 2018 detaille, equations, donnees post-2018, methodologie reproduction |

---

*Document genere dans le cadre du projet JANUS-S*
