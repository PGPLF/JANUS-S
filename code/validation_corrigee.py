#!/usr/bin/env python3
"""
validation_corrigee.py
Validation CORRIGEE - Utilisation de MU_SH0ES pour Pantheon+

CORRECTION:
- Pantheon+ m_b_corr = mu + M_B (PAS le module de distance)
- Il faut utiliser MU_SH0ES qui EST le module de distance

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

C_LIGHT = 299792.458
H0 = 70.0

# ==============================================================================
# MODELES
# ==============================================================================

class JanusCosmology:
    def __init__(self, H0=70.0):
        self.H0 = H0
        self.c = C_LIGHT

    def distance_modulus(self, z, q0):
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
    def __init__(self, H0=70.0, Omega_m=0.3, Omega_L=0.7):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_L = Omega_L
        self.c = C_LIGHT

    def E(self, z):
        return np.sqrt(self.Omega_m * (1 + z)**3 + self.Omega_L)

    def comoving_distance(self, z):
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
        d_c = self.comoving_distance(z)
        return d_c * (1 + z)

    def distance_modulus(self, z):
        d_L = self.luminosity_distance(z)
        return 5 * np.log10(d_L) + 25


# ==============================================================================
# CHARGEMENT CORRIGE
# ==============================================================================

def load_jla():
    """Charge JLA avec mu calcule"""
    print("\n" + "=" * 60)
    print("CHARGEMENT JLA")
    print("=" * 60)

    filepath = 'data/jla/jla_likelihood_v6/data/jla_lcparams.txt'
    data = pd.read_csv(filepath, sep=r'\s+')

    alpha = 0.141
    beta = 3.101
    M_B = -19.05

    mu = data['mb'] - M_B + alpha * data['x1'] - beta * data['color']
    sigma = np.sqrt(data['dmb']**2 + (alpha * data['dx1'])**2 + (beta * data['dcolor'])**2)

    print(f"  N SNe: {len(data)}")
    print(f"  z: [{data['zcmb'].min():.4f}, {data['zcmb'].max():.4f}]")
    print(f"  mu: [{mu.min():.2f}, {mu.max():.2f}]")
    print(f"  Calcul: mu = mb - M_B + alpha*x1 - beta*color")

    return {
        'z': data['zcmb'].values,
        'mu': mu.values,
        'sigma': sigma.values,
        'n': len(data),
        'name': 'JLA'
    }


def load_pantheon_corrected():
    """
    Charge Pantheon+ avec MU_SH0ES (module de distance CORRECT)
    et MU_SH0ES_ERR_DIAG pour les erreurs
    """
    print("\n" + "=" * 60)
    print("CHARGEMENT PANTHEON+ (CORRIGE)")
    print("=" * 60)

    filepath = 'data/pantheon/Pantheon+SH0ES.dat'
    data = pd.read_csv(filepath, sep=r'\s+')

    print(f"  Entrees brutes: {len(data)}")
    print(f"  SNe uniques: {data['CID'].nunique()}")

    # CORRECTION: Utiliser MU_SH0ES comme module de distance
    print(f"\n  AVANT: m_b_corr range = [{data['m_b_corr'].min():.2f}, {data['m_b_corr'].max():.2f}]")
    print(f"  APRES: MU_SH0ES range = [{data['MU_SH0ES'].min():.2f}, {data['MU_SH0ES'].max():.2f}]")
    print(f"  -> MU_SH0ES est le module de distance correct!")

    # Moyenne ponderee par SNe unique en utilisant MU_SH0ES
    grouped = data.groupby('CID')
    z_list, mu_list, sigma_list = [], [], []

    for cid, group in grouped:
        # Utiliser MU_SH0ES_ERR_DIAG pour les poids
        weights = 1 / group['MU_SH0ES_ERR_DIAG'].values**2
        z_mean = np.average(group['zHD'].values, weights=weights)
        mu_mean = np.average(group['MU_SH0ES'].values, weights=weights)
        sigma_mean = 1 / np.sqrt(np.sum(weights))

        z_list.append(z_mean)
        mu_list.append(mu_mean)
        sigma_list.append(sigma_mean)

    z_arr = np.array(z_list)
    mu_arr = np.array(mu_list)
    sigma_arr = np.array(sigma_list)

    print(f"\n  Apres moyenne ponderee: {len(z_arr)} SNe uniques")
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

    return {'q0': q0, 'offset': offset, 'chi2': chi2_val, 'dof': dof,
            'chi2_red': chi2_val/dof, 'model': model, 'name': name}


def fit_lcdm(z, mu_obs, sigma_mu, name="", Omega_m=0.3):
    model = LambdaCDM(H0=70.0, Omega_m=Omega_m, Omega_L=1-Omega_m)
    mu_th_base = model.distance_modulus(z)

    def chi2(offset):
        mu_th = mu_th_base + offset
        return np.sum(((mu_obs - mu_th) / sigma_mu)**2)

    result = minimize(chi2, [0.0], method='Nelder-Mead')
    offset = result.x[0]
    chi2_val = result.fun
    dof = len(z) - 1

    return {'Omega_m': Omega_m, 'offset': offset, 'chi2': chi2_val, 'dof': dof,
            'chi2_red': chi2_val/dof, 'model': model, 'name': name}


def bootstrap_error(z, mu, sigma, n_boot=100):
    """Erreur bootstrap sur q0"""
    n = len(z)
    q0_samples = []
    for _ in range(n_boot):
        idx = np.random.choice(n, n, replace=True)
        try:
            res = fit_janus(z[idx], mu[idx], sigma[idx])
            if res['chi2_red'] < 10:
                q0_samples.append(res['q0'])
        except:
            pass
    return np.std(q0_samples) if len(q0_samples) > 10 else np.nan


# ==============================================================================
# VISUALISATION
# ==============================================================================

def plot_hubble_comparison(jla, pantheon, janus_jla, janus_pan, output_path):
    """Diagramme de Hubble comparatif"""
    fig, axes = plt.subplots(2, 1, figsize=(14, 12), height_ratios=[2, 1])

    # Hubble diagram
    ax = axes[0]
    ax.errorbar(jla['z'], jla['mu'], yerr=jla['sigma'], fmt='o', ms=3, alpha=0.4,
                color='blue', elinewidth=0.3, label=f"JLA (n={jla['n']})")
    ax.errorbar(pantheon['z'], pantheon['mu'], yerr=pantheon['sigma'], fmt='o', ms=3, alpha=0.4,
                color='red', elinewidth=0.3, label=f"Pantheon+ (n={pantheon['n']})")

    z_model = np.logspace(-2.5, 0.4, 200)

    mu_jla = janus_jla['model'].distance_modulus(z_model, janus_jla['q0']) + janus_jla['offset']
    mu_pan = janus_pan['model'].distance_modulus(z_model, janus_pan['q0']) + janus_pan['offset']

    ax.plot(z_model, mu_jla, 'b-', lw=2.5, label=f"JANUS JLA: q0={janus_jla['q0']:.4f}")
    ax.plot(z_model, mu_pan, 'r--', lw=2.5, label=f"JANUS Panth: q0={janus_pan['q0']:.4f}")

    ax.set_xscale('log')
    ax.set_ylabel('Module de distance mu (mag)', fontsize=12)
    ax.set_title('Diagramme de Hubble - Comparaison JLA et Pantheon+ (CORRIGE)', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Residuals
    ax = axes[1]
    mu_th_jla = janus_jla['model'].distance_modulus(jla['z'], janus_jla['q0']) + janus_jla['offset']
    mu_th_pan = janus_pan['model'].distance_modulus(pantheon['z'], janus_pan['q0']) + janus_pan['offset']

    ax.scatter(jla['z'], jla['mu'] - mu_th_jla, s=5, alpha=0.4, color='blue', label='JLA')
    ax.scatter(pantheon['z'], pantheon['mu'] - mu_th_pan, s=5, alpha=0.4, color='red', label='Pantheon+')
    ax.axhline(0, color='black', lw=1)
    ax.set_xscale('log')
    ax.set_xlabel('Redshift z', fontsize=12)
    ax.set_ylabel('Residus (mag)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-1.5, 1.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure: {output_path}")


def plot_janus_vs_lcdm(data, janus, lcdm, output_path):
    """Comparaison JANUS vs LCDM"""
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), height_ratios=[2, 1])

    ax = axes[0]
    ax.errorbar(data['z'], data['mu'], yerr=data['sigma'], fmt='o', ms=2, alpha=0.3,
                color='gray', elinewidth=0.3, label=f"{data['name']} (n={data['n']})")

    z_model = np.logspace(np.log10(data['z'].min()), np.log10(data['z'].max()), 200)

    mu_janus = janus['model'].distance_modulus(z_model, janus['q0']) + janus['offset']
    mu_lcdm = lcdm['model'].distance_modulus(z_model) + lcdm['offset']

    ax.plot(z_model, mu_janus, 'b-', lw=2.5,
            label=f"JANUS: q0={janus['q0']:.4f}, chi2/dof={janus['chi2_red']:.4f}")
    ax.plot(z_model, mu_lcdm, 'r--', lw=2.5,
            label=f"LCDM: Om={lcdm['Omega_m']:.2f}, chi2/dof={lcdm['chi2_red']:.4f}")

    ax.set_xscale('log')
    ax.set_ylabel('mu (mag)')
    ax.set_title(f'JANUS vs Lambda-CDM - {data["name"]}')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)

    # Residus
    ax = axes[1]
    mu_j = janus['model'].distance_modulus(data['z'], janus['q0']) + janus['offset']
    mu_l = lcdm['model'].distance_modulus(data['z']) + lcdm['offset']

    ax.scatter(data['z'], data['mu'] - mu_j, s=3, alpha=0.4, c='blue', label='JANUS')
    ax.scatter(data['z'], data['mu'] - mu_l, s=3, alpha=0.4, c='red', label='LCDM')
    ax.axhline(0, color='black', lw=1)
    ax.set_xscale('log')
    ax.set_xlabel('z')
    ax.set_ylabel('Residus')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-1, 1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Figure: {output_path}")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("\n" + "#" * 70)
    print("# VALIDATION CORRIGEE - ANALYSE COMPLETE")
    print("# Correction: Utilisation de MU_SH0ES pour Pantheon+")
    print("#" * 70)

    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(project_dir)
    os.makedirs('results/figures', exist_ok=True)

    # Chargement
    jla = load_jla()
    pantheon = load_pantheon_corrected()

    # =========================================================================
    # AJUSTEMENTS JANUS
    # =========================================================================
    print("\n" + "=" * 60)
    print("AJUSTEMENTS JANUS")
    print("=" * 60)

    janus_jla = fit_janus(jla['z'], jla['mu'], jla['sigma'], "JLA")
    print(f"\n  JLA:")
    print(f"    q0 = {janus_jla['q0']:.6f}")
    print(f"    offset = {janus_jla['offset']:.4f}")
    print(f"    chi2/dof = {janus_jla['chi2_red']:.4f}")

    janus_pan = fit_janus(pantheon['z'], pantheon['mu'], pantheon['sigma'], "Pantheon+")
    print(f"\n  Pantheon+:")
    print(f"    q0 = {janus_pan['q0']:.6f}")
    print(f"    offset = {janus_pan['offset']:.4f}")
    print(f"    chi2/dof = {janus_pan['chi2_red']:.4f}")

    # Bootstrap
    print("\n  Bootstrap (100 samples)...")
    err_jla = bootstrap_error(jla['z'], jla['mu'], jla['sigma'], 100)
    err_pan = bootstrap_error(pantheon['z'], pantheon['mu'], pantheon['sigma'], 100)
    print(f"    sigma(q0) JLA: {err_jla:.4f}")
    print(f"    sigma(q0) Pantheon+: {err_pan:.4f}")

    # Comparaison
    delta_q0 = abs(janus_jla['q0'] - janus_pan['q0'])
    print(f"\n  Delta q0 = {delta_q0:.4f}")
    print(f"  Coherence (seuil 0.02): {'OUI' if delta_q0 < 0.02 else 'NON'}")

    # =========================================================================
    # AJUSTEMENTS LAMBDA-CDM
    # =========================================================================
    print("\n" + "=" * 60)
    print("AJUSTEMENTS LAMBDA-CDM")
    print("=" * 60)

    lcdm_jla = fit_lcdm(jla['z'], jla['mu'], jla['sigma'], "JLA", Omega_m=0.3)
    print(f"\n  JLA (Omega_m=0.3 fixe):")
    print(f"    offset = {lcdm_jla['offset']:.4f}")
    print(f"    chi2/dof = {lcdm_jla['chi2_red']:.4f}")

    lcdm_pan = fit_lcdm(pantheon['z'], pantheon['mu'], pantheon['sigma'], "Pantheon+", Omega_m=0.3)
    print(f"\n  Pantheon+ (Omega_m=0.3 fixe):")
    print(f"    offset = {lcdm_pan['offset']:.4f}")
    print(f"    chi2/dof = {lcdm_pan['chi2_red']:.4f}")

    # =========================================================================
    # COMPARAISON STATISTIQUE
    # =========================================================================
    print("\n" + "=" * 60)
    print("COMPARAISON JANUS vs LAMBDA-CDM")
    print("=" * 60)

    for name, janus, lcdm, n in [("JLA", janus_jla, lcdm_jla, jla['n']),
                                   ("Pantheon+", janus_pan, lcdm_pan, pantheon['n'])]:
        k_j, k_l = 2, 1
        aic_j = janus['chi2'] + 2*k_j
        aic_l = lcdm['chi2'] + 2*k_l
        bic_j = janus['chi2'] + k_j*np.log(n)
        bic_l = lcdm['chi2'] + k_l*np.log(n)

        print(f"\n  {name}:")
        print(f"    chi2 JANUS: {janus['chi2']:.2f}, LCDM: {lcdm['chi2']:.2f}")
        print(f"    Delta chi2 = {lcdm['chi2'] - janus['chi2']:+.2f}")
        print(f"    Delta AIC = {aic_l - aic_j:+.2f}")
        print(f"    Delta BIC = {bic_l - bic_j:+.2f}")
        pref = "JANUS" if (aic_l - aic_j) > 0 else "LCDM"
        print(f"    -> Preference: {pref}")

    # =========================================================================
    # SOUS-ECHANTILLONS PANTHEON+
    # =========================================================================
    print("\n" + "=" * 60)
    print("ANALYSE PAR SOUS-ECHANTILLONS (Pantheon+)")
    print("=" * 60)

    for z_max, label in [(0.1, "z < 0.1"), (0.5, "z < 0.5"), (1.0, "z < 1.0"), (1.3, "z < 1.3")]:
        mask = pantheon['z'] < z_max
        if np.sum(mask) > 50:
            res = fit_janus(pantheon['z'][mask], pantheon['mu'][mask], pantheon['sigma'][mask])
            print(f"  {label} (n={np.sum(mask)}): q0 = {res['q0']:.4f}, chi2/dof = {res['chi2_red']:.4f}")

    # =========================================================================
    # FIGURES
    # =========================================================================
    print("\n" + "=" * 60)
    print("GENERATION DES FIGURES")
    print("=" * 60)

    plot_hubble_comparison(jla, pantheon, janus_jla, janus_pan,
                           'results/figures/hubble_jla_pantheon_corrige.png')
    plot_janus_vs_lcdm(jla, janus_jla, lcdm_jla,
                       'results/figures/janus_vs_lcdm_jla_v2.png')
    plot_janus_vs_lcdm(pantheon, janus_pan, lcdm_pan,
                       'results/figures/janus_vs_lcdm_pantheon_v2.png')

    # =========================================================================
    # RESUME
    # =========================================================================
    print("\n" + "=" * 70)
    print("RESUME FINAL (CORRIGE)")
    print("=" * 70)
    print(f"""
    ╔══════════════════════════════════════════════════════════════════╗
    ║                    RESULTATS JANUS (CORRIGE)                     ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║  Dataset       │ N SNe │ q0        │ sigma(q0) │ chi2/dof       ║
    ║  JLA           │ {jla['n']:>5} │ {janus_jla['q0']:>9.4f} │ {err_jla:>9.4f} │ {janus_jla['chi2_red']:>9.4f}      ║
    ║  Pantheon+     │ {pantheon['n']:>5} │ {janus_pan['q0']:>9.4f} │ {err_pan:>9.4f} │ {janus_pan['chi2_red']:>9.4f}      ║
    ║  Reference 2018│   740 │   -0.0870 │    0.0150 │    0.8900      ║
    ╚══════════════════════════════════════════════════════════════════╝

    Delta q0 (JLA - Pantheon+) = {delta_q0:.4f}
    Coherence: {'OUI - Les datasets sont coherents!' if delta_q0 < 0.02 else 'NON - Ecart significatif'}

    ╔══════════════════════════════════════════════════════════════════╗
    ║                    COMPARAISON JANUS vs LCDM                     ║
    ╠══════════════════════════════════════════════════════════════════╣
    ║  Dataset       │ chi2/dof JANUS │ chi2/dof LCDM │ Preference    ║
    ║  JLA           │ {janus_jla['chi2_red']:>14.4f} │ {lcdm_jla['chi2_red']:>13.4f} │ {'JANUS' if janus_jla['chi2_red'] < lcdm_jla['chi2_red'] else 'LCDM':>13} ║
    ║  Pantheon+     │ {janus_pan['chi2_red']:>14.4f} │ {lcdm_pan['chi2_red']:>13.4f} │ {'JANUS' if janus_pan['chi2_red'] < lcdm_pan['chi2_red'] else 'LCDM':>13} ║
    ╚══════════════════════════════════════════════════════════════════╝
    """)

    return {
        'jla': jla, 'pantheon': pantheon,
        'janus_jla': janus_jla, 'janus_pan': janus_pan,
        'lcdm_jla': lcdm_jla, 'lcdm_pan': lcdm_pan,
        'err_jla': err_jla, 'err_pan': err_pan,
        'delta_q0': delta_q0
    }


if __name__ == "__main__":
    results = main()
