#!/usr/bin/env python3
"""
reproduce_2018.py
Reproduction de l'article D'Agostini & Petit (2018)
"Constraints on Janus Cosmological model from recent observations of supernovae type Ia"

Objectif: Reproduire q0 = -0.087 ± 0.015, chi2/dof ≈ 0.89

Auteur: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
import os
import sys

# Configuration
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12

# ==============================================================================
# MODELE JANUS
# ==============================================================================

class JanusCosmology:
    """Modele cosmologique JANUS"""

    def __init__(self, H0=70.0):
        self.H0 = H0
        self.c = 299792.458  # km/s

    def distance_modulus(self, z, q0):
        """
        Module de distance theorique mu(z) dans le modele JANUS
        Equation de D'Agostini & Petit (2018)
        """
        z = np.atleast_1d(z)

        # Condition de validite
        condition = 1 + 2 * q0 * z
        if np.any(condition <= 0):
            return np.full_like(z, np.inf, dtype=float)

        sqrt_term = np.sqrt(condition)

        # Distance luminosite adimensionnee
        # D_L = (c/H0) * z * [1 + z(1-q0)/(1 + q0*z + sqrt(1+2*q0*z))]
        term = z + z**2 * (1 - q0) / (1 + q0 * z + sqrt_term)

        # Distance luminosite en Mpc
        d_L = (self.c / self.H0) * term

        # Module de distance
        mu = 5 * np.log10(d_L) + 25

        return mu

    def universe_age(self, q0):
        """Age de l'univers en Gyr"""
        if q0 >= 0:
            return np.nan

        # u0 a l'epoque actuelle
        u0 = np.arcsinh(np.sqrt(-1 / (2 * q0)))

        # Temps de Hubble en Gyr
        t_H = 977.8 / self.H0  # Gyr

        # Age
        term = (1 + np.sinh(2 * u0) / 2 + u0)
        T0 = t_H * term / (2 * (-q0)**1.5)

        return T0


# ==============================================================================
# CHARGEMENT DES DONNEES
# ==============================================================================

def load_jla_data(filepath):
    """Charge le dataset JLA"""
    print(f"Chargement de {filepath}...")
    data = pd.read_csv(filepath, sep=r'\s+')
    print(f"  -> {len(data)} supernovae chargees")
    return data


def compute_distance_modulus(data, alpha=0.141, beta=3.101, M_B=-19.05):
    """
    Calcule le module de distance standardise
    mu = m_B - M_B + alpha * x1 - beta * color
    """
    mu = data['mb'] - M_B + alpha * data['x1'] - beta * data['color']

    # Propagation des erreurs
    sigma_mu = np.sqrt(
        data['dmb']**2 +
        (alpha * data['dx1'])**2 +
        (beta * data['dcolor'])**2
    )

    return mu.values, sigma_mu.values


# ==============================================================================
# AJUSTEMENT
# ==============================================================================

def chi2_function(params, z, mu_obs, sigma_mu, model):
    """Fonction chi2 a minimiser"""
    q0, offset = params

    # Contraintes sur q0
    if q0 >= 0 or q0 < -0.5:
        return 1e10

    mu_theory = model.distance_modulus(z, q0) + offset

    if np.any(~np.isfinite(mu_theory)):
        return 1e10

    residuals = (mu_obs - mu_theory) / sigma_mu
    return np.sum(residuals**2)


def fit_janus(z, mu_obs, sigma_mu, H0=70.0):
    """Ajuste le modele JANUS aux donnees"""
    model = JanusCosmology(H0=H0)

    # Estimation initiale
    initial = [-0.1, 0.0]

    # Optimisation
    result = minimize(
        chi2_function,
        initial,
        args=(z, mu_obs, sigma_mu, model),
        method='Nelder-Mead',
        options={'maxiter': 10000, 'xatol': 1e-8, 'fatol': 1e-8}
    )

    q0_fit = result.x[0]
    offset_fit = result.x[1]
    chi2_min = result.fun
    n_data = len(z)
    dof = n_data - 2

    # Age de l'univers
    T0 = model.universe_age(q0_fit)

    return {
        'q0': q0_fit,
        'offset': offset_fit,
        'chi2': chi2_min,
        'dof': dof,
        'chi2_reduced': chi2_min / dof,
        'age_Gyr': T0,
        'success': result.success,
        'model': model
    }


# ==============================================================================
# VISUALISATION
# ==============================================================================

def plot_hubble_diagram(z, mu_obs, sigma_mu, result, output_path):
    """Trace le diagramme de Hubble"""
    fig, ax = plt.subplots(figsize=(12, 8))

    # Donnees
    ax.errorbar(z, mu_obs, yerr=sigma_mu, fmt='o', markersize=3,
                alpha=0.5, color='blue', label=f'JLA (n={len(z)})',
                elinewidth=0.5, capsize=0)

    # Modele JANUS
    z_model = np.logspace(np.log10(z.min()), np.log10(z.max()), 200)
    mu_model = result['model'].distance_modulus(z_model, result['q0']) + result['offset']
    ax.plot(z_model, mu_model, 'r-', linewidth=2,
            label=f"JANUS: q₀ = {result['q0']:.4f}")

    ax.set_xscale('log')
    ax.set_xlabel('Redshift z', fontsize=14)
    ax.set_ylabel('Module de distance μ (mag)', fontsize=14)
    ax.set_title('Diagramme de Hubble - Reproduction D\'Agostini & Petit (2018)', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure sauvegardee: {output_path}")


def plot_residuals(z, mu_obs, sigma_mu, result, output_path):
    """Trace les residus"""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Calcul des residus
    mu_model = result['model'].distance_modulus(z, result['q0']) + result['offset']
    residuals = mu_obs - mu_model

    # Plot
    ax.errorbar(z, residuals, yerr=sigma_mu, fmt='o', markersize=3,
                alpha=0.5, color='blue', elinewidth=0.5, capsize=0)
    ax.axhline(y=0, color='r', linestyle='-', linewidth=2)

    # Statistiques
    mean_res = np.mean(residuals)
    std_res = np.std(residuals)
    ax.axhline(y=mean_res, color='green', linestyle='--', linewidth=1,
               label=f'Moyenne: {mean_res:.4f}')
    ax.fill_between([z.min(), z.max()],
                    [mean_res - std_res, mean_res - std_res],
                    [mean_res + std_res, mean_res + std_res],
                    alpha=0.2, color='green', label=f'±1σ: {std_res:.3f}')

    ax.set_xscale('log')
    ax.set_xlabel('Redshift z', fontsize=14)
    ax.set_ylabel('Résidus (mag)', fontsize=14)
    ax.set_title('Résidus μ_obs - μ_JANUS', fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure sauvegardee: {output_path}")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("\n" + "=" * 70)
    print("REPRODUCTION DE L'ARTICLE D'AGOSTINI & PETIT (2018)")
    print("Constraints on Janus Cosmological model from SNe Ia")
    print("=" * 70)

    # Repertoire du projet
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(project_dir)
    print(f"\nRepertoire: {os.getcwd()}")

    # Creer le dossier results/figures si necessaire
    os.makedirs('results/figures', exist_ok=True)

    # -------------------------------------------------------------------------
    # ETAPE 1: Chargement des donnees
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 1: Chargement des donnees JLA")
    print("-" * 70)

    jla_path = 'data/jla/jla_likelihood_v6/data/jla_lcparams.txt'
    data = load_jla_data(jla_path)

    z = data['zcmb'].values
    print(f"  Redshift: z_min = {z.min():.4f}, z_max = {z.max():.4f}")

    # -------------------------------------------------------------------------
    # ETAPE 2: Calcul des modules de distance
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 2: Calcul des modules de distance standardises")
    print("-" * 70)

    # Parametres de nuisance (Betoule 2014)
    alpha = 0.141
    beta = 3.101
    M_B = -19.05

    print(f"  Parametres: alpha = {alpha}, beta = {beta}, M_B = {M_B}")

    mu_obs, sigma_mu = compute_distance_modulus(data, alpha, beta, M_B)
    print(f"  mu: min = {mu_obs.min():.2f}, max = {mu_obs.max():.2f}")
    print(f"  sigma_mu moyen: {sigma_mu.mean():.3f}")

    # -------------------------------------------------------------------------
    # ETAPE 3: Ajustement du modele JANUS
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 3: Ajustement du modele JANUS")
    print("-" * 70)

    result = fit_janus(z, mu_obs, sigma_mu, H0=70.0)

    print(f"\n  RESULTATS:")
    print(f"  -" * 35)
    print(f"  q0          = {result['q0']:.6f}")
    print(f"  offset      = {result['offset']:.6f}")
    print(f"  chi2        = {result['chi2']:.2f}")
    print(f"  dof         = {result['dof']}")
    print(f"  chi2/dof    = {result['chi2_reduced']:.4f}")
    print(f"  Age univers = {result['age_Gyr']:.2f} Gyr")
    print(f"  Convergence = {result['success']}")

    # -------------------------------------------------------------------------
    # ETAPE 4: Validation
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 4: Validation vs Article 2018")
    print("-" * 70)

    q0_ref = -0.087
    chi2_dof_ref = 657 / 738  # ≈ 0.89

    delta_q0 = abs(result['q0'] - q0_ref)
    delta_chi2 = abs(result['chi2_reduced'] - chi2_dof_ref)

    print(f"\n  Comparaison:")
    print(f"  {'Parametre':<15} {'Obtenu':<15} {'Reference':<15} {'Ecart':<15}")
    print(f"  {'-'*60}")
    print(f"  {'q0':<15} {result['q0']:<15.4f} {q0_ref:<15.4f} {delta_q0:<15.4f}")
    print(f"  {'chi2/dof':<15} {result['chi2_reduced']:<15.4f} {chi2_dof_ref:<15.4f} {delta_chi2:<15.4f}")

    # Validation
    q0_ok = delta_q0 < 0.02
    chi2_ok = delta_chi2 < 0.1

    print(f"\n  Validation:")
    print(f"  q0 dans tolerance (±0.02):      {'✅ OUI' if q0_ok else '❌ NON'}")
    print(f"  chi2/dof dans tolerance (±0.1): {'✅ OUI' if chi2_ok else '❌ NON'}")

    # -------------------------------------------------------------------------
    # ETAPE 5: Generation des figures
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 5: Generation des figures")
    print("-" * 70)

    plot_hubble_diagram(z, mu_obs, sigma_mu, result,
                        'results/figures/hubble_jla_janus.png')

    plot_residuals(z, mu_obs, sigma_mu, result,
                   'results/figures/residus_jla_janus.png')

    # -------------------------------------------------------------------------
    # RESUME
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("RESUME - REPRODUCTION 2018")
    print("=" * 70)
    print(f"""
    Dataset:        JLA (740 SNe Ia)

    Resultats:
      q0          = {result['q0']:.4f} (ref: -0.087)
      chi2        = {result['chi2']:.2f}
      chi2/dof    = {result['chi2_reduced']:.4f} (ref: 0.89)
      Age univers = {result['age_Gyr']:.1f} Gyr (H0=70)

    Figures generees:
      - results/figures/hubble_jla_janus.png
      - results/figures/residus_jla_janus.png

    Statut: {'✅ REPRODUCTION REUSSIE' if q0_ok and chi2_ok else '⚠️ ECARTS DETECTES'}
    """)
    print("=" * 70)

    return result


if __name__ == "__main__":
    result = main()
