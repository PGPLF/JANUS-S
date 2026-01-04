# Plan de Publication Comparative
## Reproduction et Validation du Modele JANUS sur les Supernovae de Type Ia

**Version :** 1.0
**Date :** 4 janvier 2026
**Auteurs :** [A completer]
**Journal cible :** Astrophysics and Space Science / European Physical Journal C

---

## ANALYSE DU PLAN INITIAL

### Points forts
- Methodologie technique detaillee et actionable
- Code Python fonctionnel et reproductible
- References aux sources de donnees correctes
- Parametres de nuisance documentes

### Points a ameliorer pour une publication academique

| Aspect | Plan initial | Amelioration requise |
|--------|--------------|---------------------|
| Structure | Tutoriel technique | Format article scientifique |
| Comparaison | Reproduction simple | Analyse comparative rigoureuse |
| Statistiques | Chi-2 basique | Tests statistiques complets (AIC, BIC, likelihood ratio) |
| Incertitudes | Erreurs simples | Analyse MCMC, bootstrap, systematics |
| Donnees | JLA 2014 uniquement | Ajout Pantheon+ 2022 pour validation |
| Discussion | Validation resultats | Interpretation physique, limites, perspectives |

---

## PLAN DE PUBLICATION AMELIORE

### PHASE 1 : PREPARATION (Semaine 1-2)

#### 1.1 Acquisition des ressources

**Articles de reference :**
- [ ] D'Agostini & Petit (2018) - Article original
- [ ] Betoule et al. (2014) - Dataset JLA
- [ ] Scolnic et al. (2018) - Pantheon
- [ ] Brout et al. (2022) - Pantheon+
- [ ] Petit et al. (2024) - Article EPJC recent

**Datasets :**
| Dataset | Lien | Nombre SNe | Usage |
|---------|------|------------|-------|
| JLA | supernovae.in2p3.fr | 740 | Reproduction 2018 |
| Pantheon+ | github.com/PantheonPlusSH0ES | 1550 | Validation etendue |

**Fichiers requis JLA :**
- `jla_lcparams.txt` : Parametres courbes de lumiere
- `C_stat.dat`, `C_cal.dat`, etc. : Matrices covariance

**Fichiers requis Pantheon+ :**
- `Pantheon+SH0ES.dat` : Donnees completes
- Matrices covariance associees

#### 1.2 Environnement technique

```bash
# Environnement Python
conda create -n janus python=3.10
conda activate janus
pip install numpy scipy pandas matplotlib astropy emcee corner sncosmo
```

---

### PHASE 2 : REPRODUCTION DE 2018 (Semaine 2-3)

#### 2.1 Implementation du modele JANUS

**Equations a implementer :**

```python
# janus_model.py

import numpy as np
from scipy.optimize import minimize

class JanusCosmology:
    """
    Implementation du modele cosmologique JANUS
    Base sur D'Agostini & Petit (2018)
    """

    def __init__(self, H0=70.0):
        self.H0 = H0  # km/s/Mpc
        self.c = 299792.458  # km/s

    def distance_modulus(self, z, q0, M_B=-19.05):
        """
        Module de distance theorique mu(z) dans JCM

        Parameters:
        -----------
        z : array-like
            Redshift
        q0 : float
            Parametre de deceleration (< 0 pour acceleration)
        M_B : float
            Magnitude absolue de reference

        Returns:
        --------
        mu : array-like
            Module de distance
        """
        # Condition de validite
        if np.any(1 + 2 * q0 * z <= 0):
            return np.full_like(z, np.inf)

        # Equation (XX) de l'article
        sqrt_term = np.sqrt(1 + 2 * q0 * z)
        numerator = z + z**2 * (1 - q0) / (1 + q0 * z + sqrt_term)

        # Distance luminosite en Mpc
        d_L = (self.c / self.H0) * numerator

        # Module de distance
        mu = 5 * np.log10(d_L) + 25

        return mu

    def comoving_distance(self, z, q0):
        """Distance comobile r(z)"""
        factor = self.c / self.H0
        num = q0 * z + (1 - q0) * (1 - np.sqrt(1 + 2 * q0 * z))
        den = q0**2 * (1 + z)
        return factor * (num / den)

    def universe_age(self, q0):
        """Age de l'univers T0 en Gyr"""
        if q0 >= 0:
            return np.nan
        u0 = np.arcsinh(np.sqrt(-1 / (2 * q0)))
        # Conversion en Gyr
        H0_per_Gyr = self.H0 * 1.022e-3  # km/s/Mpc -> 1/Gyr
        T0 = (1 / H0_per_Gyr) * (1 + np.sinh(2*u0)/2 + u0) / (2 * (-q0)**1.5)
        return T0
```

#### 2.2 Chargement et preparation des donnees

```python
# data_loader.py

import pandas as pd
import numpy as np

def load_jla_data(filepath='data/jla_lcparams.txt'):
    """
    Charge le dataset JLA (740 SNe Ia)
    """
    data = pd.read_csv(filepath, delim_whitespace=True)

    return {
        'z': data['zcmb'].values,
        'mb': data['mb'].values,
        'x1': data['x1'].values,
        'color': data['color'].values,
        'dmb': data['dmb'].values,
        'dx1': data['dx1'].values,
        'dcolor': data['dcolor'].values,
        'n_sne': len(data)
    }

def compute_distance_modulus(data, alpha=0.141, beta=3.101, M_B=-19.05):
    """
    Calcule le module de distance standardise
    mu = m_B - M_B + alpha * x1 - beta * C
    """
    mu = data['mb'] - M_B + alpha * data['x1'] - beta * data['color']

    # Propagation des erreurs (simplifiee)
    sigma_mu = np.sqrt(
        data['dmb']**2 +
        (alpha * data['dx1'])**2 +
        (beta * data['dcolor'])**2
    )

    return mu, sigma_mu
```

#### 2.3 Ajustement et validation

**Objectifs de reproduction :**
| Parametre | Valeur 2018 | Tolerance |
|-----------|-------------|-----------|
| q0 | -0.087 | ± 0.015 |
| chi2/dof | 0.89 | ± 0.05 |
| T0 (H0=70) | ~15 Gyr | ± 1 Gyr |

```python
# fitting.py

from scipy.optimize import minimize, differential_evolution
import emcee

def chi2_simple(params, z, mu_obs, sigma_mu, model):
    """Chi-2 avec erreurs diagonales"""
    q0, offset = params
    mu_theory = model.distance_modulus(z, q0) + offset
    return np.sum(((mu_obs - mu_theory) / sigma_mu)**2)

def chi2_covariance(params, z, mu_obs, cov_matrix, model):
    """Chi-2 avec matrice de covariance complete"""
    q0, offset = params
    mu_theory = model.distance_modulus(z, q0) + offset
    residuals = mu_obs - mu_theory
    cov_inv = np.linalg.inv(cov_matrix)
    return residuals @ cov_inv @ residuals

def fit_janus(z, mu_obs, sigma_mu, model, method='Nelder-Mead'):
    """Ajustement du modele JANUS"""
    initial = [-0.1, 0.0]
    bounds = [(-0.5, -0.001), (-2, 2)]

    result = minimize(
        chi2_simple,
        initial,
        args=(z, mu_obs, sigma_mu, model),
        method=method,
        bounds=bounds
    )

    return {
        'q0': result.x[0],
        'offset': result.x[1],
        'chi2': result.fun,
        'dof': len(z) - 2,
        'chi2_reduced': result.fun / (len(z) - 2),
        'success': result.success
    }
```

---

### PHASE 3 : EXTENSION PANTHEON+ (Semaine 3-4)

#### 3.1 Adaptation pour Pantheon+

- Charger le dataset Pantheon+ (1550 SNe)
- Adapter le format des donnees
- Appliquer le meme protocole d'ajustement
- Comparer les resultats

#### 3.2 Tests de robustesse

| Test | Description |
|------|-------------|
| Sous-echantillons | Ajuster sur z < 0.5, z > 0.5 separement |
| Bootstrap | 1000 reechantillonnages pour erreurs |
| MCMC | Exploration complete de l'espace des parametres |
| Leave-one-out | Sensibilite aux outliers |

---

### PHASE 4 : COMPARAISON AVEC LAMBDA-CDM (Semaine 4-5)

#### 4.1 Implementation Lambda-CDM

```python
def mu_lcdm(z, Omega_m=0.3, Omega_L=0.7, H0=70):
    """Module de distance dans Lambda-CDM (integration numerique)"""
    from scipy.integrate import quad

    def E(z, Om, OL):
        return np.sqrt(Om * (1+z)**3 + OL)

    def integrand(z, Om, OL):
        return 1 / E(z, Om, OL)

    c = 299792.458  # km/s
    d_H = c / H0

    mu = np.zeros_like(z, dtype=float)
    for i, zi in enumerate(z):
        integral, _ = quad(integrand, 0, zi, args=(Omega_m, Omega_L))
        d_L = d_H * (1 + zi) * integral
        mu[i] = 5 * np.log10(d_L) + 25

    return mu
```

#### 4.2 Criteres de comparaison

| Critere | Formule | Interpretation |
|---------|---------|----------------|
| Chi-2 reduit | chi2 / dof | < 1 : bon ajustement |
| AIC | chi2 + 2k | Plus petit = meilleur (penalise params) |
| BIC | chi2 + k*ln(n) | Plus petit = meilleur (penalise plus) |
| Delta AIC | AIC_LCDM - AIC_JANUS | > 10 : forte evidence |

```python
def model_comparison(chi2_janus, chi2_lcdm, n_data, k_janus=2, k_lcdm=3):
    """Comparaison statistique des modeles"""

    # AIC
    aic_janus = chi2_janus + 2 * k_janus
    aic_lcdm = chi2_lcdm + 2 * k_lcdm
    delta_aic = aic_lcdm - aic_janus

    # BIC
    bic_janus = chi2_janus + k_janus * np.log(n_data)
    bic_lcdm = chi2_lcdm + k_lcdm * np.log(n_data)
    delta_bic = bic_lcdm - bic_janus

    # Likelihood ratio
    lr = np.exp(-0.5 * (chi2_janus - chi2_lcdm))

    return {
        'AIC_JANUS': aic_janus,
        'AIC_LCDM': aic_lcdm,
        'Delta_AIC': delta_aic,
        'BIC_JANUS': bic_janus,
        'BIC_LCDM': bic_lcdm,
        'Delta_BIC': delta_bic,
        'Likelihood_ratio': lr
    }
```

---

### PHASE 5 : ANALYSE MCMC (Semaine 5-6)

#### 5.1 Exploration bayesienne

```python
import emcee
import corner

def log_likelihood(params, z, mu_obs, sigma_mu, model):
    q0, offset = params
    if q0 >= 0 or q0 < -0.5:
        return -np.inf
    mu_theory = model.distance_modulus(z, q0) + offset
    chi2 = np.sum(((mu_obs - mu_theory) / sigma_mu)**2)
    return -0.5 * chi2

def log_prior(params):
    q0, offset = params
    if -0.5 < q0 < 0 and -5 < offset < 5:
        return 0.0
    return -np.inf

def log_probability(params, z, mu_obs, sigma_mu, model):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, z, mu_obs, sigma_mu, model)

def run_mcmc(z, mu_obs, sigma_mu, model, nwalkers=32, nsteps=5000):
    """Execute l'analyse MCMC"""
    ndim = 2
    initial = np.array([-0.1, 0.0])
    pos = initial + 1e-3 * np.random.randn(nwalkers, ndim)

    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_probability,
        args=(z, mu_obs, sigma_mu, model)
    )
    sampler.run_mcmc(pos, nsteps, progress=True)

    # Burn-in
    flat_samples = sampler.get_chain(discard=1000, thin=15, flat=True)

    # Resultats
    q0_mcmc = np.percentile(flat_samples[:, 0], [16, 50, 84])
    offset_mcmc = np.percentile(flat_samples[:, 1], [16, 50, 84])

    return {
        'samples': flat_samples,
        'q0_median': q0_mcmc[1],
        'q0_err_low': q0_mcmc[1] - q0_mcmc[0],
        'q0_err_high': q0_mcmc[2] - q0_mcmc[1],
        'sampler': sampler
    }
```

---

### PHASE 6 : REDACTION DE L'ARTICLE (Semaine 6-8)

#### 6.1 Structure de l'article

```
1. ABSTRACT (150-200 mots)
   - Contexte : test du modele JANUS sur SNe Ia
   - Methode : reproduction 2018 + extension Pantheon+
   - Resultats : valeurs q0, comparaison LCDM
   - Conclusion : validation/infirmation du modele

2. INTRODUCTION (1-2 pages)
   - Probleme de l'energie noire et matiere noire
   - Modeles alternatifs : masses negatives, bimetrique
   - Historique du modele JANUS (Petit 1977-2024)
   - Objectifs de cette etude

3. THEORETICAL FRAMEWORK (2-3 pages)
   3.1 Le modele cosmologique JANUS
       - Equations de champ bimetriques
       - Solution exacte pour l'ere de poussiere
       - Relation magnitude-redshift
   3.2 Le modele standard Lambda-CDM
       - Equations de Friedmann
       - Parametres cosmologiques
   3.3 Supernovae de type Ia comme chandelles standard
       - Standardisation des courbes de lumiere
       - Module de distance

4. DATA AND METHODS (2-3 pages)
   4.1 Datasets
       - JLA (740 SNe, Betoule 2014)
       - Pantheon+ (1550 SNe, Brout 2022)
   4.2 Preprocessing
       - Standardisation (alpha, beta, M_B)
       - Corrections (biais, selection)
   4.3 Fitting procedure
       - Minimisation chi-2
       - Analyse MCMC
   4.4 Model comparison
       - AIC, BIC, likelihood ratio

5. RESULTS (3-4 pages)
   5.1 Reproduction de l'article 2018
       - Ajustement JLA
       - Comparaison avec resultats publies
   5.2 Extension Pantheon+
       - Nouveaux ajustements
       - Tests de robustesse
   5.3 Comparaison JANUS vs Lambda-CDM
       - Criteres statistiques
       - Diagrammes de Hubble
       - Residus

6. DISCUSSION (2-3 pages)
   6.1 Interpretation des resultats
   6.2 Avantages du modele JANUS
   6.3 Limites et critiques
   6.4 Coherence avec autres observables (CMB, BAO)
   6.5 Predictions testables

7. CONCLUSION (0.5-1 page)
   - Resume des resultats
   - Implications pour la cosmologie
   - Perspectives futures

ACKNOWLEDGMENTS
REFERENCES
APPENDIX : Code et donnees supplementaires
```

#### 6.2 Figures requises

| Figure | Description |
|--------|-------------|
| Fig. 1 | Schema du modele JANUS (deux secteurs) |
| Fig. 2 | Diagramme de Hubble JLA avec fit JANUS |
| Fig. 3 | Residus JLA |
| Fig. 4 | Diagramme de Hubble Pantheon+ avec fit JANUS |
| Fig. 5 | Residus Pantheon+ |
| Fig. 6 | Comparaison JANUS vs Lambda-CDM |
| Fig. 7 | Contours MCMC (q0 vs offset) |
| Fig. 8 | Distribution posterieure de q0 |

#### 6.3 Tables requises

| Table | Description |
|-------|-------------|
| Tab. 1 | Parametres des datasets |
| Tab. 2 | Resultats ajustement JLA (reproduction 2018) |
| Tab. 3 | Resultats ajustement Pantheon+ |
| Tab. 4 | Comparaison statistique JANUS vs Lambda-CDM |
| Tab. 5 | Tests de robustesse |

---

### PHASE 7 : PRODUCTION PDF (Semaine 8)

#### 7.1 Format LaTeX

```latex
\documentclass[twocolumn]{aastex631}
% ou \documentclass{aa} pour A&A

\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{natbib}

\begin{document}

\title{Reproduction and Validation of the Janus Cosmological Model
       Constraints from Type Ia Supernovae: A Comparative Study}

\author{[Authors]}
\affiliation{[Affiliations]}

\begin{abstract}
...
\end{abstract}

\keywords{cosmology: theory -- cosmology: observations --
          supernovae: general -- dark energy}

% Corps de l'article
...

\bibliography{references}

\end{document}
```

#### 7.2 Compilation

```bash
# Compilation PDF
pdflatex article.tex
bibtex article
pdflatex article.tex
pdflatex article.tex
```

---

## CALENDRIER RECAPITULATIF

| Phase | Semaine | Livrables |
|-------|---------|-----------|
| 1. Preparation | 1-2 | Datasets, environnement, code base |
| 2. Reproduction | 2-3 | Resultats JLA, validation 2018 |
| 3. Extension | 3-4 | Resultats Pantheon+, robustesse |
| 4. Comparaison | 4-5 | JANUS vs LCDM, criteres stat |
| 5. MCMC | 5-6 | Analyse bayesienne, erreurs |
| 6. Redaction | 6-8 | Article complet, figures, tables |
| 7. Production | 8 | PDF final, soumission |

---

## CRITERES DE SUCCES

### Reproduction 2018
- [ ] q0 = -0.087 ± 0.02
- [ ] chi2/dof < 1.0
- [ ] Diagrammes similaires aux Figs. 3-7

### Extension
- [ ] Ajustement Pantheon+ convergent
- [ ] Tests robustesse coherents
- [ ] MCMC bien explore

### Comparaison
- [ ] Delta AIC/BIC calcules
- [ ] Interpretation claire
- [ ] Discussion equilibree

### Publication
- [ ] Format journal respecte
- [ ] Figures haute qualite
- [ ] References completes
- [ ] Code reproductible (GitHub)

---

## FICHIERS DU PROJET

```
JANUS-S/
├── publications/
│   ├── PLAN_PUBLICATION_COMPARATIVE.md  (ce fichier)
│   ├── article/
│   │   ├── article.tex
│   │   ├── figures/
│   │   └── references.bib
│   └── drafts/
├── code/
│   ├── janus_model.py
│   ├── data_loader.py
│   ├── fitting.py
│   ├── comparison.py
│   ├── mcmc_analysis.py
│   └── notebooks/
│       ├── 01_data_exploration.ipynb
│       ├── 02_reproduction_2018.ipynb
│       ├── 03_pantheon_extension.ipynb
│       └── 04_model_comparison.ipynb
├── data/
│   ├── jla/
│   └── pantheon/
├── results/
│   ├── figures/
│   └── tables/
├── synthese/
├── references/
└── LOG.md
```

---

*Plan de publication v1.0 - Projet JANUS-S*
