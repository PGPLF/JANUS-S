#!/usr/bin/env python3
"""
extend_pantheon.py
Phase 3 : Extension du modele JANUS au dataset Pantheon+ (1701 SNe Ia)

Objectifs:
1. Ajuster JANUS sur Pantheon+
2. Tests de robustesse (bootstrap, sous-echantillons)
3. Comparer avec resultats JLA
4. Generer figures comparatives

Auteur: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import os
import sys

# Configuration
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# ==============================================================================
# MODELE JANUS (meme que reproduce_2018.py)
# ==============================================================================

class JanusCosmology:
    """Modele cosmologique JANUS"""

    def __init__(self, H0=70.0):
        self.H0 = H0
        self.c = 299792.458  # km/s

    def distance_modulus(self, z, q0):
        """
        Module de distance theorique mu(z) dans le modele JANUS
        """
        z = np.atleast_1d(z)

        condition = 1 + 2 * q0 * z
        if np.any(condition <= 0):
            return np.full_like(z, np.inf, dtype=float)

        sqrt_term = np.sqrt(condition)
        term = z + z**2 * (1 - q0) / (1 + q0 * z + sqrt_term)
        d_L = (self.c / self.H0) * term
        mu = 5 * np.log10(d_L) + 25

        return mu


# ==============================================================================
# CHARGEMENT PANTHEON+
# ==============================================================================

def load_pantheon_data(filepath):
    """Charge le dataset Pantheon+ avec magnitudes corrigees"""
    print(f"Chargement de {filepath}...")
    data = pd.read_csv(filepath, sep=r'\s+')
    print(f"  -> {len(data)} entrees chargees")

    # Colonnes disponibles
    print(f"  Colonnes: {list(data.columns[:15])}...")

    return data


def prepare_pantheon_unique(data):
    """
    Prepare les donnees Pantheon+ en gardant une entree par SNe
    (moyenne ponderee si plusieurs observations)
    """
    # Grouper par CID (identifiant unique de la supernova)
    grouped = data.groupby('CID')

    z_list = []
    mu_list = []
    sigma_list = []
    names_list = []

    for cid, group in grouped:
        # Moyenne ponderee des observations
        weights = 1 / group['m_b_corr_err_DIAG'].values**2

        z_mean = np.average(group['zHD'].values, weights=weights)
        mu_mean = np.average(group['m_b_corr'].values, weights=weights)
        sigma_mean = 1 / np.sqrt(np.sum(weights))

        z_list.append(z_mean)
        mu_list.append(mu_mean)
        sigma_list.append(sigma_mean)
        names_list.append(cid)

    result = {
        'z': np.array(z_list),
        'mu': np.array(mu_list),
        'sigma': np.array(sigma_list),
        'names': np.array(names_list),
        'n_sne': len(z_list)
    }

    print(f"  -> {result['n_sne']} supernovae uniques")
    print(f"  Redshift: z_min = {result['z'].min():.4f}, z_max = {result['z'].max():.4f}")

    return result


# ==============================================================================
# AJUSTEMENT
# ==============================================================================

def chi2_function(params, z, mu_obs, sigma_mu, model):
    """Fonction chi2 a minimiser"""
    q0, offset = params

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

    initial = [-0.1, 0.0]

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

    return {
        'q0': q0_fit,
        'offset': offset_fit,
        'chi2': chi2_min,
        'dof': dof,
        'chi2_reduced': chi2_min / dof,
        'success': result.success,
        'model': model
    }


# ==============================================================================
# TESTS DE ROBUSTESSE
# ==============================================================================

def bootstrap_analysis(z, mu_obs, sigma_mu, n_bootstrap=100, H0=70.0):
    """Analyse bootstrap pour estimer les erreurs"""
    print(f"\n  Bootstrap avec {n_bootstrap} echantillons...")

    n = len(z)
    q0_samples = []

    for i in range(n_bootstrap):
        indices = np.random.choice(n, n, replace=True)
        z_boot = z[indices]
        mu_boot = mu_obs[indices]
        sigma_boot = sigma_mu[indices]

        try:
            result = fit_janus(z_boot, mu_boot, sigma_boot, H0)
            if result['success']:
                q0_samples.append(result['q0'])
        except:
            continue

        if (i + 1) % 20 == 0:
            print(f"    {i+1}/{n_bootstrap} echantillons traites")

    q0_samples = np.array(q0_samples)

    return {
        'q0_mean': np.mean(q0_samples),
        'q0_std': np.std(q0_samples),
        'q0_median': np.median(q0_samples),
        'q0_16': np.percentile(q0_samples, 16),
        'q0_84': np.percentile(q0_samples, 84),
        'n_success': len(q0_samples)
    }


def subsample_analysis(z, mu_obs, sigma_mu, H0=70.0):
    """Analyse par sous-echantillons en redshift"""
    print("\n  Analyse par sous-echantillons...")

    results = {}

    # z < 0.5 (low-z)
    mask_low = z < 0.5
    if np.sum(mask_low) > 50:
        result_low = fit_janus(z[mask_low], mu_obs[mask_low], sigma_mu[mask_low], H0)
        results['low_z'] = {
            'z_range': f'z < 0.5 (n={np.sum(mask_low)})',
            'q0': result_low['q0'],
            'chi2_red': result_low['chi2_reduced']
        }
        print(f"    z < 0.5: q0 = {result_low['q0']:.4f}")

    # z >= 0.5 (high-z)
    mask_high = z >= 0.5
    if np.sum(mask_high) > 50:
        result_high = fit_janus(z[mask_high], mu_obs[mask_high], sigma_mu[mask_high], H0)
        results['high_z'] = {
            'z_range': f'z >= 0.5 (n={np.sum(mask_high)})',
            'q0': result_high['q0'],
            'chi2_red': result_high['chi2_reduced']
        }
        print(f"    z >= 0.5: q0 = {result_high['q0']:.4f}")

    # z < 0.1 (very low-z)
    mask_vlow = z < 0.1
    if np.sum(mask_vlow) > 30:
        result_vlow = fit_janus(z[mask_vlow], mu_obs[mask_vlow], sigma_mu[mask_vlow], H0)
        results['very_low_z'] = {
            'z_range': f'z < 0.1 (n={np.sum(mask_vlow)})',
            'q0': result_vlow['q0'],
            'chi2_red': result_vlow['chi2_reduced']
        }
        print(f"    z < 0.1: q0 = {result_vlow['q0']:.4f}")

    return results


# ==============================================================================
# VISUALISATION
# ==============================================================================

def plot_hubble_pantheon(z, mu_obs, sigma_mu, result, output_path):
    """Trace le diagramme de Hubble pour Pantheon+"""
    fig, ax = plt.subplots(figsize=(14, 9))

    # Donnees
    ax.errorbar(z, mu_obs, yerr=sigma_mu, fmt='o', markersize=2,
                alpha=0.4, color='blue', label=f'Pantheon+ (n={len(z)})',
                elinewidth=0.3, capsize=0)

    # Modele JANUS
    z_model = np.logspace(np.log10(z.min()), np.log10(z.max()), 300)
    mu_model = result['model'].distance_modulus(z_model, result['q0']) + result['offset']
    ax.plot(z_model, mu_model, 'r-', linewidth=2.5,
            label=f"JANUS: q0 = {result['q0']:.4f}")

    ax.set_xscale('log')
    ax.set_xlabel('Redshift z', fontsize=14)
    ax.set_ylabel('Module de distance mu (mag)', fontsize=14)
    ax.set_title('Diagramme de Hubble - Pantheon+ avec modele JANUS', fontsize=16)
    ax.legend(fontsize=12, loc='lower right')
    ax.grid(True, alpha=0.3)

    # Annotation
    textstr = f'chi2/dof = {result["chi2_reduced"]:.4f}'
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure sauvegardee: {output_path}")


def plot_residuals_pantheon(z, mu_obs, sigma_mu, result, output_path):
    """Trace les residus pour Pantheon+"""
    fig, ax = plt.subplots(figsize=(14, 6))

    mu_model = result['model'].distance_modulus(z, result['q0']) + result['offset']
    residuals = mu_obs - mu_model

    ax.errorbar(z, residuals, yerr=sigma_mu, fmt='o', markersize=2,
                alpha=0.4, color='blue', elinewidth=0.3, capsize=0)
    ax.axhline(y=0, color='r', linestyle='-', linewidth=2)

    mean_res = np.mean(residuals)
    std_res = np.std(residuals)
    ax.axhline(y=mean_res, color='green', linestyle='--', linewidth=1,
               label=f'Moyenne: {mean_res:.4f}')
    ax.fill_between([z.min(), z.max()],
                    [mean_res - std_res, mean_res - std_res],
                    [mean_res + std_res, mean_res + std_res],
                    alpha=0.2, color='green', label=f'+/-1sigma: {std_res:.3f}')

    ax.set_xscale('log')
    ax.set_xlabel('Redshift z', fontsize=14)
    ax.set_ylabel('Residus (mag)', fontsize=14)
    ax.set_title('Residus mu_obs - mu_JANUS (Pantheon+)', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure sauvegardee: {output_path}")


def plot_comparison_jla_pantheon(result_jla, result_pantheon, output_path):
    """Compare les resultats JLA et Pantheon+"""
    fig, ax = plt.subplots(figsize=(10, 6))

    datasets = ['JLA\n(740 SNe)', 'Pantheon+\n(unique SNe)']
    q0_values = [result_jla['q0'], result_pantheon['q0']]
    q0_errors = [result_jla.get('q0_err', 0.015), result_pantheon.get('q0_err', 0.01)]

    colors = ['#2ecc71', '#3498db']
    x = np.arange(len(datasets))

    bars = ax.bar(x, q0_values, color=colors, width=0.6, edgecolor='black', linewidth=1.5)
    ax.errorbar(x, q0_values, yerr=q0_errors, fmt='none', color='black', capsize=8, capthick=2)

    # Reference 2018
    ax.axhline(y=-0.087, color='red', linestyle='--', linewidth=2, label='Ref. 2018: q0 = -0.087')

    ax.set_ylabel('q0 (parametre de deceleration)', fontsize=14)
    ax.set_title('Comparaison des resultats JANUS: JLA vs Pantheon+', fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')

    # Annotations des valeurs
    for i, (bar, val) in enumerate(zip(bars, q0_values)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 0.005,
                f'{val:.4f}', ha='center', va='top', fontsize=12, fontweight='bold', color='white')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure sauvegardee: {output_path}")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("\n" + "=" * 70)
    print("PHASE 3 : EXTENSION PANTHEON+")
    print("Ajustement du modele JANUS sur 1701 supernovae")
    print("=" * 70)

    # Repertoire du projet
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(project_dir)
    print(f"\nRepertoire: {os.getcwd()}")

    os.makedirs('results/figures', exist_ok=True)

    # -------------------------------------------------------------------------
    # ETAPE 1: Chargement des donnees Pantheon+
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 1: Chargement des donnees Pantheon+")
    print("-" * 70)

    pantheon_path = 'data/pantheon/Pantheon+SH0ES.dat'
    data = load_pantheon_data(pantheon_path)

    # Preparer les donnees uniques
    pantheon = prepare_pantheon_unique(data)

    z = pantheon['z']
    mu_obs = pantheon['mu']
    sigma_mu = pantheon['sigma']

    # -------------------------------------------------------------------------
    # ETAPE 2: Ajustement JANUS sur Pantheon+
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 2: Ajustement du modele JANUS")
    print("-" * 70)

    result_pantheon = fit_janus(z, mu_obs, sigma_mu, H0=70.0)

    print(f"\n  RESULTATS PANTHEON+:")
    print(f"  " + "-" * 40)
    print(f"  q0          = {result_pantheon['q0']:.6f}")
    print(f"  offset      = {result_pantheon['offset']:.6f}")
    print(f"  chi2        = {result_pantheon['chi2']:.2f}")
    print(f"  dof         = {result_pantheon['dof']}")
    print(f"  chi2/dof    = {result_pantheon['chi2_reduced']:.4f}")
    print(f"  Convergence = {result_pantheon['success']}")

    # -------------------------------------------------------------------------
    # ETAPE 3: Tests de robustesse
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 3: Tests de robustesse")
    print("-" * 70)

    # Bootstrap
    bootstrap = bootstrap_analysis(z, mu_obs, sigma_mu, n_bootstrap=100, H0=70.0)
    result_pantheon['q0_err'] = bootstrap['q0_std']

    print(f"\n  Resultats bootstrap:")
    print(f"    q0 = {bootstrap['q0_mean']:.4f} +/- {bootstrap['q0_std']:.4f}")
    print(f"    Intervalle 68%: [{bootstrap['q0_16']:.4f}, {bootstrap['q0_84']:.4f}]")

    # Sous-echantillons
    subsamples = subsample_analysis(z, mu_obs, sigma_mu, H0=70.0)

    # -------------------------------------------------------------------------
    # ETAPE 4: Comparaison avec JLA
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 4: Comparaison avec resultats JLA")
    print("-" * 70)

    # Resultats JLA de Phase 2
    result_jla = {
        'q0': -0.0864,
        'q0_err': 0.015,
        'chi2_reduced': 0.8834,
        'dof': 738,
        'n_sne': 740
    }

    delta_q0 = result_pantheon['q0'] - result_jla['q0']

    print(f"\n  Comparaison:")
    print(f"  {'Dataset':<15} {'q0':<12} {'chi2/dof':<12} {'N_SNe':<10}")
    print(f"  {'-'*50}")
    print(f"  {'JLA':<15} {result_jla['q0']:<12.4f} {result_jla['chi2_reduced']:<12.4f} {result_jla['n_sne']:<10}")
    print(f"  {'Pantheon+':<15} {result_pantheon['q0']:<12.4f} {result_pantheon['chi2_reduced']:<12.4f} {pantheon['n_sne']:<10}")
    print(f"  {'-'*50}")
    print(f"  Delta q0: {delta_q0:.4f}")

    # Coherence
    coherent = abs(delta_q0) < 0.02
    print(f"\n  Coherence JLA/Pantheon+: {'OUI' if coherent else 'NON'} (seuil: +/-0.02)")

    # -------------------------------------------------------------------------
    # ETAPE 5: Generation des figures
    # -------------------------------------------------------------------------
    print("\n" + "-" * 70)
    print("ETAPE 5: Generation des figures")
    print("-" * 70)

    plot_hubble_pantheon(z, mu_obs, sigma_mu, result_pantheon,
                         'results/figures/hubble_pantheon_janus.png')

    plot_residuals_pantheon(z, mu_obs, sigma_mu, result_pantheon,
                            'results/figures/residus_pantheon_janus.png')

    plot_comparison_jla_pantheon(result_jla, result_pantheon,
                                  'results/figures/comparaison_jla_pantheon.png')

    # -------------------------------------------------------------------------
    # RESUME
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("RESUME - PHASE 3: EXTENSION PANTHEON+")
    print("=" * 70)
    print(f"""
    Dataset:        Pantheon+ ({pantheon['n_sne']} SNe uniques)
    Plage z:        [{z.min():.4f}, {z.max():.4f}]

    Resultats JANUS:
      q0          = {result_pantheon['q0']:.4f} +/- {bootstrap['q0_std']:.4f}
      chi2/dof    = {result_pantheon['chi2_reduced']:.4f}

    Comparaison avec JLA:
      q0 JLA      = {result_jla['q0']:.4f}
      q0 Pantheon = {result_pantheon['q0']:.4f}
      Delta       = {delta_q0:.4f}
      Coherence   = {'OUI' if coherent else 'NON'}

    Tests de robustesse:
      Bootstrap (100 samples): q0 = {bootstrap['q0_mean']:.4f} +/- {bootstrap['q0_std']:.4f}

    Figures generees:
      - results/figures/hubble_pantheon_janus.png
      - results/figures/residus_pantheon_janus.png
      - results/figures/comparaison_jla_pantheon.png

    Statut: {'PHASE 3 VALIDEE' if coherent else 'ECARTS DETECTES'}
    """)
    print("=" * 70)

    return {
        'pantheon': result_pantheon,
        'bootstrap': bootstrap,
        'subsamples': subsamples,
        'jla': result_jla,
        'coherent': coherent
    }


if __name__ == "__main__":
    results = main()
