#!/usr/bin/env python3
"""
validation_complete.py
Validation complete des donnees et calculs + Phase 4

Objectifs:
1. Verifier les donnees brutes JLA et Pantheon+
2. Valider les formules du modele JANUS
3. Comparer sur meme plage de redshift
4. Ajuster Lambda-CDM (Phase 4)
5. Comparaison statistique JANUS vs Lambda-CDM

Auteur: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import quad
import os
import sys

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

# ==============================================================================
# CONSTANTES
# ==============================================================================

C_LIGHT = 299792.458  # km/s
H0 = 70.0  # km/s/Mpc

# ==============================================================================
# MODELES COSMOLOGIQUES
# ==============================================================================

class JanusCosmology:
    """Modele JANUS - D'Agostini & Petit 2018"""

    def __init__(self, H0=70.0):
        self.H0 = H0
        self.c = C_LIGHT

    def distance_modulus(self, z, q0):
        """Module de distance JANUS"""
        z = np.atleast_1d(z)

        condition = 1 + 2 * q0 * z
        if np.any(condition <= 0):
            return np.full_like(z, np.inf, dtype=float)

        sqrt_term = np.sqrt(condition)
        term = z + z**2 * (1 - q0) / (1 + q0 * z + sqrt_term)
        d_L = (self.c / self.H0) * term
        mu = 5 * np.log10(d_L) + 25

        return mu


class LambdaCDM:
    """Modele Lambda-CDM standard"""

    def __init__(self, H0=70.0, Omega_m=0.3, Omega_L=0.7):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_L = Omega_L
        self.c = C_LIGHT

    def E(self, z):
        """Fonction E(z) = H(z)/H0"""
        return np.sqrt(self.Omega_m * (1 + z)**3 + self.Omega_L)

    def comoving_distance(self, z):
        """Distance comobile en Mpc"""
        if np.isscalar(z):
            integral, _ = quad(lambda zp: 1/self.E(zp), 0, z)
            return (self.c / self.H0) * integral
        else:
            result = np.zeros_like(z, dtype=float)
            for i, zi in enumerate(z):
                integral, _ = quad(lambda zp: 1/self.E(zp), 0, zi)
                result[i] = (self.c / self.H0) * integral
            return result

    def luminosity_distance(self, z):
        """Distance luminosite en Mpc"""
        d_c = self.comoving_distance(z)
        return d_c * (1 + z)

    def distance_modulus(self, z):
        """Module de distance"""
        d_L = self.luminosity_distance(z)
        return 5 * np.log10(d_L) + 25


# ==============================================================================
# CHARGEMENT ET VALIDATION DES DONNEES
# ==============================================================================

def load_and_validate_jla(filepath):
    """Charge et valide JLA"""
    print("\n" + "=" * 60)
    print("VALIDATION JLA")
    print("=" * 60)

    data = pd.read_csv(filepath, sep=r'\s+')

    print(f"Fichier: {filepath}")
    print(f"Lignes: {len(data)}")
    print(f"Colonnes: {list(data.columns)}")

    # Verifications
    assert len(data) == 740, f"Attendu 740, trouve {len(data)}"
    assert 'zcmb' in data.columns
    assert 'mb' in data.columns
    assert 'x1' in data.columns
    assert 'color' in data.columns

    # Statistiques
    print(f"\nRedshift (zcmb):")
    print(f"  min: {data['zcmb'].min():.4f}")
    print(f"  max: {data['zcmb'].max():.4f}")
    print(f"  mean: {data['zcmb'].mean():.4f}")

    print(f"\nMagnitude (mb):")
    print(f"  min: {data['mb'].min():.2f}")
    print(f"  max: {data['mb'].max():.2f}")
    print(f"  mean: {data['mb'].mean():.2f}")

    # Calculer mu standardise
    alpha = 0.141
    beta = 3.101
    M_B = -19.05

    mu = data['mb'] - M_B + alpha * data['x1'] - beta * data['color']
    sigma = np.sqrt(data['dmb']**2 + (alpha * data['dx1'])**2 + (beta * data['dcolor'])**2)

    print(f"\nModule de distance (mu = mb - M_B + alpha*x1 - beta*color):")
    print(f"  alpha = {alpha}, beta = {beta}, M_B = {M_B}")
    print(f"  mu min: {mu.min():.2f}")
    print(f"  mu max: {mu.max():.2f}")
    print(f"  mu mean: {mu.mean():.2f}")
    print(f"  sigma mean: {sigma.mean():.3f}")

    return {
        'z': data['zcmb'].values,
        'mu': mu.values,
        'sigma': sigma.values,
        'n': len(data),
        'name': 'JLA'
    }


def load_and_validate_pantheon(filepath):
    """Charge et valide Pantheon+"""
    print("\n" + "=" * 60)
    print("VALIDATION PANTHEON+")
    print("=" * 60)

    data = pd.read_csv(filepath, sep=r'\s+')

    print(f"Fichier: {filepath}")
    print(f"Lignes: {len(data)}")
    print(f"Colonnes principales: {list(data.columns[:20])}")

    # Verifications
    assert 'zHD' in data.columns, "Colonne zHD manquante"
    assert 'm_b_corr' in data.columns, "Colonne m_b_corr manquante"

    print(f"\nNombre d'entrees: {len(data)}")
    print(f"SNe uniques (CID): {data['CID'].nunique()}")

    # Statistiques sur donnees brutes
    print(f"\nRedshift (zHD):")
    print(f"  min: {data['zHD'].min():.4f}")
    print(f"  max: {data['zHD'].max():.4f}")
    print(f"  mean: {data['zHD'].mean():.4f}")

    print(f"\nMagnitude corrigee (m_b_corr):")
    print(f"  min: {data['m_b_corr'].min():.2f}")
    print(f"  max: {data['m_b_corr'].max():.2f}")
    print(f"  mean: {data['m_b_corr'].mean():.2f}")

    print(f"\nMU_SH0ES (module de distance fourni):")
    print(f"  min: {data['MU_SH0ES'].min():.2f}")
    print(f"  max: {data['MU_SH0ES'].max():.2f}")
    print(f"  mean: {data['MU_SH0ES'].mean():.2f}")

    # IMPORTANT: Verifier la relation entre m_b_corr et MU_SH0ES
    diff = data['m_b_corr'] - data['MU_SH0ES']
    print(f"\nDifference m_b_corr - MU_SH0ES:")
    print(f"  mean: {diff.mean():.4f}")
    print(f"  std: {diff.std():.4f}")
    print(f"  -> Cela correspond a M_B !")

    # Moyenne ponderee par SNe unique
    grouped = data.groupby('CID')
    z_list, mu_list, sigma_list = [], [], []

    for cid, group in grouped:
        weights = 1 / group['m_b_corr_err_DIAG'].values**2
        z_mean = np.average(group['zHD'].values, weights=weights)
        mu_mean = np.average(group['m_b_corr'].values, weights=weights)
        sigma_mean = 1 / np.sqrt(np.sum(weights))

        z_list.append(z_mean)
        mu_list.append(mu_mean)
        sigma_list.append(sigma_mean)

    z_arr = np.array(z_list)
    mu_arr = np.array(mu_list)
    sigma_arr = np.array(sigma_list)

    print(f"\nApres moyenne ponderee:")
    print(f"  N SNe uniques: {len(z_arr)}")
    print(f"  z: [{z_arr.min():.4f}, {z_arr.max():.4f}]")
    print(f"  mu: [{mu_arr.min():.2f}, {mu_arr.max():.2f}]")

    return {
        'z': z_arr,
        'mu': mu_arr,
        'sigma': sigma_arr,
        'n': len(z_arr),
        'name': 'Pantheon+'
    }


# ==============================================================================
# AJUSTEMENTS
# ==============================================================================

def fit_janus(z, mu_obs, sigma_mu, name=""):
    """Ajuste JANUS"""
    model = JanusCosmology(H0=70.0)

    def chi2(params):
        q0, offset = params
        if q0 >= 0 or q0 < -0.5:
            return 1e10
        mu_th = model.distance_modulus(z, q0) + offset
        if np.any(~np.isfinite(mu_th)):
            return 1e10
        return np.sum(((mu_obs - mu_th) / sigma_mu)**2)

    result = minimize(chi2, [-0.1, 0.0], method='Nelder-Mead',
                      options={'maxiter': 10000, 'xatol': 1e-8})

    q0, offset = result.x
    chi2_val = result.fun
    dof = len(z) - 2

    print(f"\n  JANUS sur {name}:")
    print(f"    q0 = {q0:.6f}")
    print(f"    offset = {offset:.4f}")
    print(f"    chi2 = {chi2_val:.2f}, chi2/dof = {chi2_val/dof:.4f}")

    return {'q0': q0, 'offset': offset, 'chi2': chi2_val, 'dof': dof,
            'chi2_red': chi2_val/dof, 'model': model, 'name': name}


def fit_lcdm(z, mu_obs, sigma_mu, name="", Omega_m=0.3):
    """Ajuste Lambda-CDM (offset seulement, Omega_m fixe)"""
    model = LambdaCDM(H0=70.0, Omega_m=Omega_m, Omega_L=1-Omega_m)

    # Precalculer mu theorique
    mu_th_base = model.distance_modulus(z)

    def chi2(offset):
        mu_th = mu_th_base + offset
        return np.sum(((mu_obs - mu_th) / sigma_mu)**2)

    result = minimize(chi2, [0.0], method='Nelder-Mead')

    offset = result.x[0]
    chi2_val = result.fun
    dof = len(z) - 1  # 1 parametre libre (offset)

    print(f"\n  Lambda-CDM sur {name} (Omega_m={Omega_m}):")
    print(f"    offset = {offset:.4f}")
    print(f"    chi2 = {chi2_val:.2f}, chi2/dof = {chi2_val/dof:.4f}")

    return {'Omega_m': Omega_m, 'offset': offset, 'chi2': chi2_val, 'dof': dof,
            'chi2_red': chi2_val/dof, 'model': model, 'name': name}


def fit_lcdm_with_omega(z, mu_obs, sigma_mu, name=""):
    """Ajuste Lambda-CDM avec Omega_m libre"""
    def chi2(params):
        Omega_m, offset = params
        if Omega_m <= 0 or Omega_m >= 1:
            return 1e10
        model = LambdaCDM(H0=70.0, Omega_m=Omega_m, Omega_L=1-Omega_m)
        mu_th = model.distance_modulus(z) + offset
        return np.sum(((mu_obs - mu_th) / sigma_mu)**2)

    result = minimize(chi2, [0.3, 0.0], method='Nelder-Mead',
                      options={'maxiter': 10000})

    Omega_m, offset = result.x
    chi2_val = result.fun
    dof = len(z) - 2

    print(f"\n  Lambda-CDM (Omega_m libre) sur {name}:")
    print(f"    Omega_m = {Omega_m:.4f}")
    print(f"    offset = {offset:.4f}")
    print(f"    chi2 = {chi2_val:.2f}, chi2/dof = {chi2_val/dof:.4f}")

    return {'Omega_m': Omega_m, 'offset': offset, 'chi2': chi2_val, 'dof': dof,
            'chi2_red': chi2_val/dof, 'name': name}


# ==============================================================================
# COMPARAISONS
# ==============================================================================

def compare_same_z_range(jla, pantheon):
    """Compare JLA et Pantheon+ sur meme plage de z"""
    print("\n" + "=" * 60)
    print("COMPARAISON SUR MEME PLAGE DE REDSHIFT")
    print("=" * 60)

    # Plage commune: z dans [0.01, 1.30] (plage JLA)
    z_min, z_max = 0.01, 1.30

    # Filtrer JLA
    mask_jla = (jla['z'] >= z_min) & (jla['z'] <= z_max)
    jla_filt = {
        'z': jla['z'][mask_jla],
        'mu': jla['mu'][mask_jla],
        'sigma': jla['sigma'][mask_jla]
    }
    print(f"\nJLA filtre [{z_min}, {z_max}]: {len(jla_filt['z'])} SNe")

    # Filtrer Pantheon+
    mask_pan = (pantheon['z'] >= z_min) & (pantheon['z'] <= z_max)
    pan_filt = {
        'z': pantheon['z'][mask_pan],
        'mu': pantheon['mu'][mask_pan],
        'sigma': pantheon['sigma'][mask_pan]
    }
    print(f"Pantheon+ filtre [{z_min}, {z_max}]: {len(pan_filt['z'])} SNe")

    # Ajuster JANUS sur les deux
    print("\n" + "-" * 40)
    print("Ajustements JANUS sur plage commune:")
    print("-" * 40)

    janus_jla = fit_janus(jla_filt['z'], jla_filt['mu'], jla_filt['sigma'], "JLA_filtre")
    janus_pan = fit_janus(pan_filt['z'], pan_filt['mu'], pan_filt['sigma'], "Pantheon+_filtre")

    delta_q0 = abs(janus_jla['q0'] - janus_pan['q0'])
    print(f"\n  Delta q0 = {delta_q0:.4f}")

    return janus_jla, janus_pan, jla_filt, pan_filt


def model_comparison_stats(result_janus, result_lcdm, n_data):
    """Calcule AIC, BIC et autres statistiques"""
    k_janus = 2  # q0 + offset
    k_lcdm = 1   # offset seulement (Omega_m fixe)

    chi2_janus = result_janus['chi2']
    chi2_lcdm = result_lcdm['chi2']

    # AIC = chi2 + 2k
    aic_janus = chi2_janus + 2 * k_janus
    aic_lcdm = chi2_lcdm + 2 * k_lcdm
    delta_aic = aic_lcdm - aic_janus

    # BIC = chi2 + k * ln(n)
    bic_janus = chi2_janus + k_janus * np.log(n_data)
    bic_lcdm = chi2_lcdm + k_lcdm * np.log(n_data)
    delta_bic = bic_lcdm - bic_janus

    return {
        'AIC_janus': aic_janus, 'AIC_lcdm': aic_lcdm, 'delta_AIC': delta_aic,
        'BIC_janus': bic_janus, 'BIC_lcdm': bic_lcdm, 'delta_BIC': delta_bic,
        'delta_chi2': chi2_lcdm - chi2_janus
    }


# ==============================================================================
# VISUALISATION
# ==============================================================================

def plot_validation_comparison(jla, pantheon, output_path):
    """Compare visuellement JLA et Pantheon+"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Hubble diagram overlay
    ax = axes[0, 0]
    ax.scatter(jla['z'], jla['mu'], s=5, alpha=0.5, label=f"JLA (n={jla['n']})", c='blue')
    ax.scatter(pantheon['z'], pantheon['mu'], s=5, alpha=0.5, label=f"Pantheon+ (n={pantheon['n']})", c='red')
    ax.set_xscale('log')
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('mu (mag)')
    ax.set_title('Comparaison JLA vs Pantheon+')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. Distribution des redshifts
    ax = axes[0, 1]
    ax.hist(jla['z'], bins=50, alpha=0.5, label='JLA', density=True)
    ax.hist(pantheon['z'], bins=50, alpha=0.5, label='Pantheon+', density=True)
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Densite')
    ax.set_title('Distribution des redshifts')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3. Distribution des erreurs
    ax = axes[1, 0]
    ax.hist(jla['sigma'], bins=50, alpha=0.5, label='JLA', density=True)
    ax.hist(pantheon['sigma'], bins=50, alpha=0.5, label='Pantheon+', density=True)
    ax.set_xlabel('sigma_mu (mag)')
    ax.set_ylabel('Densite')
    ax.set_title('Distribution des erreurs')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. mu vs z avec bins
    ax = axes[1, 1]
    z_bins = np.logspace(-2.5, 0.5, 20)

    for data, name, color in [(jla, 'JLA', 'blue'), (pantheon, 'Pantheon+', 'red')]:
        z_centers, mu_means, mu_stds = [], [], []
        for i in range(len(z_bins)-1):
            mask = (data['z'] >= z_bins[i]) & (data['z'] < z_bins[i+1])
            if np.sum(mask) > 3:
                z_centers.append(np.mean(data['z'][mask]))
                mu_means.append(np.mean(data['mu'][mask]))
                mu_stds.append(np.std(data['mu'][mask]))

        ax.errorbar(z_centers, mu_means, yerr=mu_stds, fmt='o-', label=name,
                    color=color, capsize=3, markersize=5)

    ax.set_xscale('log')
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('mu moyen (mag)')
    ax.set_title('Moyenne par bin de redshift')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nFigure sauvegardee: {output_path}")


def plot_model_fits(z, mu_obs, sigma_mu, janus_result, lcdm_result, name, output_path):
    """Compare les ajustements JANUS et LCDM"""
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), height_ratios=[2, 1])

    # Hubble diagram
    ax = axes[0]
    ax.errorbar(z, mu_obs, yerr=sigma_mu, fmt='o', markersize=2, alpha=0.4,
                color='gray', elinewidth=0.3, capsize=0, label=f'{name} data')

    z_model = np.logspace(np.log10(z.min()), np.log10(z.max()), 200)

    # JANUS
    mu_janus = janus_result['model'].distance_modulus(z_model, janus_result['q0']) + janus_result['offset']
    ax.plot(z_model, mu_janus, 'b-', lw=2.5,
            label=f"JANUS: q0={janus_result['q0']:.4f}, chi2/dof={janus_result['chi2_red']:.3f}")

    # LCDM
    mu_lcdm = lcdm_result['model'].distance_modulus(z_model) + lcdm_result['offset']
    ax.plot(z_model, mu_lcdm, 'r--', lw=2.5,
            label=f"LCDM: Om={lcdm_result['Omega_m']:.2f}, chi2/dof={lcdm_result['chi2_red']:.3f}")

    ax.set_xscale('log')
    ax.set_ylabel('mu (mag)')
    ax.set_title(f'Comparaison JANUS vs Lambda-CDM - {name}')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)

    # Residuals
    ax = axes[1]
    mu_janus_data = janus_result['model'].distance_modulus(z, janus_result['q0']) + janus_result['offset']
    mu_lcdm_data = lcdm_result['model'].distance_modulus(z) + lcdm_result['offset']

    res_janus = mu_obs - mu_janus_data
    res_lcdm = mu_obs - mu_lcdm_data

    ax.scatter(z, res_janus, s=3, alpha=0.4, color='blue', label='Residus JANUS')
    ax.scatter(z, res_lcdm, s=3, alpha=0.4, color='red', label='Residus LCDM')
    ax.axhline(0, color='black', lw=1)
    ax.set_xscale('log')
    ax.set_xlabel('Redshift z')
    ax.set_ylabel('Residus (mag)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-1, 1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure sauvegardee: {output_path}")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("\n" + "#" * 70)
    print("# VALIDATION COMPLETE ET PHASE 4")
    print("# Verification donnees + Comparaison JANUS vs Lambda-CDM")
    print("#" * 70)

    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(project_dir)
    os.makedirs('results/figures', exist_ok=True)

    # =========================================================================
    # PARTIE 1: VALIDATION DES DONNEES
    # =========================================================================
    print("\n" + "=" * 70)
    print("PARTIE 1: VALIDATION DES DONNEES BRUTES")
    print("=" * 70)

    jla = load_and_validate_jla('data/jla/jla_likelihood_v6/data/jla_lcparams.txt')
    pantheon = load_and_validate_pantheon('data/pantheon/Pantheon+SH0ES.dat')

    # Visualisation comparative
    plot_validation_comparison(jla, pantheon, 'results/figures/validation_jla_pantheon.png')

    # =========================================================================
    # PARTIE 2: COMPARAISON SUR MEME PLAGE Z
    # =========================================================================
    janus_jla_filt, janus_pan_filt, jla_filt, pan_filt = compare_same_z_range(jla, pantheon)

    # =========================================================================
    # PARTIE 3: AJUSTEMENTS COMPLETS
    # =========================================================================
    print("\n" + "=" * 70)
    print("PARTIE 3: AJUSTEMENTS JANUS ET LAMBDA-CDM")
    print("=" * 70)

    # JLA complet
    print("\n--- JLA complet ---")
    janus_jla = fit_janus(jla['z'], jla['mu'], jla['sigma'], "JLA")
    lcdm_jla = fit_lcdm(jla['z'], jla['mu'], jla['sigma'], "JLA", Omega_m=0.3)
    lcdm_jla_free = fit_lcdm_with_omega(jla['z'], jla['mu'], jla['sigma'], "JLA")

    # Pantheon+ complet
    print("\n--- Pantheon+ complet ---")
    janus_pantheon = fit_janus(pantheon['z'], pantheon['mu'], pantheon['sigma'], "Pantheon+")
    lcdm_pantheon = fit_lcdm(pantheon['z'], pantheon['mu'], pantheon['sigma'], "Pantheon+", Omega_m=0.3)
    lcdm_pantheon_free = fit_lcdm_with_omega(pantheon['z'], pantheon['mu'], pantheon['sigma'], "Pantheon+")

    # Figures comparatives
    plot_model_fits(jla['z'], jla['mu'], jla['sigma'], janus_jla, lcdm_jla,
                    "JLA", 'results/figures/janus_vs_lcdm_jla.png')
    plot_model_fits(pantheon['z'], pantheon['mu'], pantheon['sigma'], janus_pantheon, lcdm_pantheon,
                    "Pantheon+", 'results/figures/janus_vs_lcdm_pantheon.png')

    # =========================================================================
    # PARTIE 4: COMPARAISON STATISTIQUE
    # =========================================================================
    print("\n" + "=" * 70)
    print("PARTIE 4: COMPARAISON STATISTIQUE JANUS vs LAMBDA-CDM")
    print("=" * 70)

    stats_jla = model_comparison_stats(janus_jla, lcdm_jla, jla['n'])
    stats_pantheon = model_comparison_stats(janus_pantheon, lcdm_pantheon, pantheon['n'])

    print("\n--- JLA ---")
    print(f"  Delta chi2 (LCDM - JANUS): {stats_jla['delta_chi2']:.2f}")
    print(f"  Delta AIC:  {stats_jla['delta_AIC']:.2f}")
    print(f"  Delta BIC:  {stats_jla['delta_BIC']:.2f}")

    print("\n--- Pantheon+ ---")
    print(f"  Delta chi2 (LCDM - JANUS): {stats_pantheon['delta_chi2']:.2f}")
    print(f"  Delta AIC:  {stats_pantheon['delta_AIC']:.2f}")
    print(f"  Delta BIC:  {stats_pantheon['delta_BIC']:.2f}")

    # =========================================================================
    # RESUME FINAL
    # =========================================================================
    print("\n" + "=" * 70)
    print("RESUME FINAL - VALIDATION ET PHASE 4")
    print("=" * 70)

    print("""
    ╔══════════════════════════════════════════════════════════════════╗
    ║                    RESULTATS JANUS                               ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║  Dataset       │ N SNe │ q0        │ chi2/dof │ Offset          ║""")
    print(f"    ║  JLA           │ {jla['n']:>5} │ {janus_jla['q0']:>9.4f} │ {janus_jla['chi2_red']:>8.4f} │ {janus_jla['offset']:>+8.4f}         ║")
    print(f"    ║  Pantheon+     │ {pantheon['n']:>5} │ {janus_pantheon['q0']:>9.4f} │ {janus_pantheon['chi2_red']:>8.4f} │ {janus_pantheon['offset']:>+8.4f}         ║")
    print(f"    ║  JLA (z<1.3)   │ {len(jla_filt['z']):>5} │ {janus_jla_filt['q0']:>9.4f} │ {janus_jla_filt['chi2_red']:>8.4f} │ {janus_jla_filt['offset']:>+8.4f}         ║")
    print(f"    ║  Panth (z<1.3) │ {len(pan_filt['z']):>5} │ {janus_pan_filt['q0']:>9.4f} │ {janus_pan_filt['chi2_red']:>8.4f} │ {janus_pan_filt['offset']:>+8.4f}         ║")
    print("""    ╚══════════════════════════════════════════════════════════════════╝

    ╔══════════════════════════════════════════════════════════════════╗
    ║                    RESULTATS LAMBDA-CDM                          ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║  Dataset       │ Omega_m │ chi2/dof │ Omega_m (libre)            ║""")
    print(f"    ║  JLA           │ {lcdm_jla['Omega_m']:>7.2f} │ {lcdm_jla['chi2_red']:>8.4f} │ {lcdm_jla_free['Omega_m']:>8.4f}                    ║")
    print(f"    ║  Pantheon+     │ {lcdm_pantheon['Omega_m']:>7.2f} │ {lcdm_pantheon['chi2_red']:>8.4f} │ {lcdm_pantheon_free['Omega_m']:>8.4f}                    ║")
    print("""    ╚══════════════════════════════════════════════════════════════════╝

    ╔══════════════════════════════════════════════════════════════════╗
    ║                    COMPARAISON STATISTIQUE                       ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║  Dataset       │ Delta chi2 │ Delta AIC │ Delta BIC │ Prefere   ║""")

    pref_jla = "JANUS" if stats_jla['delta_AIC'] > 0 else "LCDM"
    pref_pan = "JANUS" if stats_pantheon['delta_AIC'] > 0 else "LCDM"
    print(f"    ║  JLA           │ {stats_jla['delta_chi2']:>+10.2f} │ {stats_jla['delta_AIC']:>+9.2f} │ {stats_jla['delta_BIC']:>+9.2f} │ {pref_jla:<9} ║")
    print(f"    ║  Pantheon+     │ {stats_pantheon['delta_chi2']:>+10.2f} │ {stats_pantheon['delta_AIC']:>+9.2f} │ {stats_pantheon['delta_BIC']:>+9.2f} │ {pref_pan:<9} ║")
    print("""    ╚══════════════════════════════════════════════════════════════════╝

    Interpretation Delta AIC/BIC:
      > +10 : Forte evidence pour JANUS
      +2 a +10 : Evidence moderee pour JANUS
      -2 a +2 : Modeles equivalents
      < -2 : Evidence pour Lambda-CDM

    Figures generees:
      - results/figures/validation_jla_pantheon.png
      - results/figures/janus_vs_lcdm_jla.png
      - results/figures/janus_vs_lcdm_pantheon.png
    """)

    # Diagnostic divergence
    print("=" * 70)
    print("DIAGNOSTIC DIVERGENCE JLA vs PANTHEON+")
    print("=" * 70)

    delta_offset = janus_pantheon['offset'] - janus_jla['offset']
    print(f"""
    Offset JANUS JLA:      {janus_jla['offset']:+.4f}
    Offset JANUS Pantheon: {janus_pantheon['offset']:+.4f}
    Difference offsets:    {delta_offset:+.4f}

    HYPOTHESE: La difference de ~19 mag dans l'offset suggere que:
    - JLA utilise mu = m_B - M_B + corrections (M_B ~ -19.05)
    - Pantheon+ m_b_corr est DEJA le module de distance mu

    -> Les magnitudes ne sont PAS comparables directement!
    -> L'offset absorbe cette difference de calibration.
    -> q0 different car poids relatifs des SNe differents.
    """)

    return {
        'janus_jla': janus_jla,
        'janus_pantheon': janus_pantheon,
        'lcdm_jla': lcdm_jla,
        'lcdm_pantheon': lcdm_pantheon,
        'stats_jla': stats_jla,
        'stats_pantheon': stats_pantheon
    }


if __name__ == "__main__":
    results = main()
