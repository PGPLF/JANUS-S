#!/usr/bin/env python3
"""
generate_publication_figures.py
Generation des figures pour la publication V0

Figures:
1. Diagramme de Hubble JLA + JANUS + LCDM
2. Diagramme de Hubble Pantheon+ + JANUS + LCDM
3. Residus comparatifs
4. Comparaison q0 par sous-echantillons
5. Tableau recapitulatif (pour LaTeX)

Auteur: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.optimize import minimize
from scipy.integrate import quad
import os

# Style publication
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

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
    def __init__(self, H0=70.0, Omega_m=0.3):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_L = 1 - Omega_m
        self.c = C_LIGHT

    def E(self, z):
        return np.sqrt(self.Omega_m * (1 + z)**3 + self.Omega_L)

    def distance_modulus(self, z):
        z = np.atleast_1d(z)
        result = np.zeros_like(z, dtype=float)
        for i, zi in enumerate(z):
            integral, _ = quad(lambda zp: 1/self.E(zp), 0, zi)
            d_L = (self.c / self.H0) * integral * (1 + zi)
            result[i] = 5 * np.log10(d_L) + 25
        return result


# ==============================================================================
# CHARGEMENT
# ==============================================================================

def load_data():
    """Charge JLA et Pantheon+"""
    # JLA
    jla_data = pd.read_csv('data/jla/jla_likelihood_v6/data/jla_lcparams.txt', sep=r'\s+')
    alpha, beta, M_B = 0.141, 3.101, -19.05
    jla_mu = jla_data['mb'] - M_B + alpha * jla_data['x1'] - beta * jla_data['color']
    jla_sigma = np.sqrt(jla_data['dmb']**2 + (alpha * jla_data['dx1'])**2 + (beta * jla_data['dcolor'])**2)
    jla = {'z': jla_data['zcmb'].values, 'mu': jla_mu.values, 'sigma': jla_sigma.values, 'n': len(jla_data)}

    # Pantheon+
    pan_data = pd.read_csv('data/pantheon/Pantheon+SH0ES.dat', sep=r'\s+')
    grouped = pan_data.groupby('CID')
    z_list, mu_list, sigma_list = [], [], []
    for cid, group in grouped:
        weights = 1 / group['MU_SH0ES_ERR_DIAG'].values**2
        z_list.append(np.average(group['zHD'].values, weights=weights))
        mu_list.append(np.average(group['MU_SH0ES'].values, weights=weights))
        sigma_list.append(1 / np.sqrt(np.sum(weights)))
    pantheon = {'z': np.array(z_list), 'mu': np.array(mu_list), 'sigma': np.array(sigma_list), 'n': len(z_list)}

    return jla, pantheon


def fit_models(z, mu, sigma):
    """Ajuste JANUS et LCDM"""
    janus = JanusCosmology()
    lcdm = LambdaCDM()

    # JANUS
    def chi2_j(p):
        q0, off = p
        if q0 >= 0 or q0 < -0.5: return 1e10
        mu_th = janus.distance_modulus(z, q0) + off
        if np.any(~np.isfinite(mu_th)): return 1e10
        return np.sum(((mu - mu_th) / sigma)**2)

    res_j = minimize(chi2_j, [-0.1, 0.0], method='Nelder-Mead')
    janus_fit = {'q0': res_j.x[0], 'offset': res_j.x[1], 'chi2': res_j.fun, 'dof': len(z)-2}

    # LCDM
    mu_lcdm_base = lcdm.distance_modulus(z)
    def chi2_l(off):
        return np.sum(((mu - mu_lcdm_base - off) / sigma)**2)
    res_l = minimize(chi2_l, [0.0], method='Nelder-Mead')
    lcdm_fit = {'offset': res_l.x[0], 'chi2': res_l.fun, 'dof': len(z)-1}

    return janus_fit, lcdm_fit, janus, lcdm


# ==============================================================================
# FIGURES
# ==============================================================================

def fig1_hubble_jla(jla, janus_fit, lcdm_fit, janus, lcdm):
    """Figure 1: Diagramme de Hubble JLA"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), height_ratios=[3, 1], sharex=True)
    fig.subplots_adjust(hspace=0.05)

    z, mu, sigma = jla['z'], jla['mu'], jla['sigma']

    # Panel superieur: Hubble diagram
    ax1.errorbar(z, mu, yerr=sigma, fmt='o', ms=2, alpha=0.4, color='gray',
                 elinewidth=0.3, capsize=0, label=f'JLA ({jla["n"]} SNe Ia)')

    z_model = np.logspace(np.log10(z.min()), np.log10(z.max()), 200)
    mu_j = janus.distance_modulus(z_model, janus_fit['q0']) + janus_fit['offset']
    mu_l = lcdm.distance_modulus(z_model) + lcdm_fit['offset']

    ax1.plot(z_model, mu_j, 'b-', lw=2, label=f'JANUS ($q_0$={janus_fit["q0"]:.3f})')
    ax1.plot(z_model, mu_l, 'r--', lw=2, label=r'$\Lambda$CDM ($\Omega_m$=0.3)')

    ax1.set_xscale('log')
    ax1.set_ylabel(r'Distance modulus $\mu$ (mag)')
    ax1.set_title('Hubble Diagram - JLA Dataset')
    ax1.legend(loc='lower right')
    ax1.grid(True, alpha=0.3)

    # Panel inferieur: residus
    mu_j_data = janus.distance_modulus(z, janus_fit['q0']) + janus_fit['offset']
    mu_l_data = lcdm.distance_modulus(z) + lcdm_fit['offset']

    ax2.scatter(z, mu - mu_j_data, s=4, alpha=0.5, c='blue', label='JANUS')
    ax2.scatter(z, mu - mu_l_data, s=4, alpha=0.5, c='red', label=r'$\Lambda$CDM')
    ax2.axhline(0, color='black', lw=0.8)
    ax2.set_xscale('log')
    ax2.set_xlabel('Redshift $z$')
    ax2.set_ylabel('Residuals (mag)')
    ax2.set_ylim(-0.8, 0.8)
    ax2.legend(loc='upper left', ncol=2)
    ax2.grid(True, alpha=0.3)

    plt.savefig('publications/article/figures/fig1_hubble_jla.pdf')
    plt.savefig('publications/article/figures/fig1_hubble_jla.png', dpi=150)
    plt.close()
    print("Figure 1: fig1_hubble_jla.pdf")


def fig2_hubble_pantheon(pantheon, janus_fit, lcdm_fit, janus, lcdm):
    """Figure 2: Diagramme de Hubble Pantheon+"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), height_ratios=[3, 1], sharex=True)
    fig.subplots_adjust(hspace=0.05)

    z, mu, sigma = pantheon['z'], pantheon['mu'], pantheon['sigma']

    ax1.errorbar(z, mu, yerr=sigma, fmt='o', ms=1.5, alpha=0.3, color='gray',
                 elinewidth=0.2, capsize=0, label=f'Pantheon+ ({pantheon["n"]} SNe Ia)')

    z_model = np.logspace(np.log10(z.min()), np.log10(z.max()), 200)
    mu_j = janus.distance_modulus(z_model, janus_fit['q0']) + janus_fit['offset']
    mu_l = lcdm.distance_modulus(z_model) + lcdm_fit['offset']

    ax1.plot(z_model, mu_j, 'b-', lw=2, label=f'JANUS ($q_0$={janus_fit["q0"]:.3f})')
    ax1.plot(z_model, mu_l, 'r--', lw=2, label=r'$\Lambda$CDM ($\Omega_m$=0.3)')

    ax1.set_xscale('log')
    ax1.set_ylabel(r'Distance modulus $\mu$ (mag)')
    ax1.set_title('Hubble Diagram - Pantheon+ Dataset')
    ax1.legend(loc='lower right')
    ax1.grid(True, alpha=0.3)

    mu_j_data = janus.distance_modulus(z, janus_fit['q0']) + janus_fit['offset']
    mu_l_data = lcdm.distance_modulus(z) + lcdm_fit['offset']

    ax2.scatter(z, mu - mu_j_data, s=2, alpha=0.4, c='blue', label='JANUS')
    ax2.scatter(z, mu - mu_l_data, s=2, alpha=0.4, c='red', label=r'$\Lambda$CDM')
    ax2.axhline(0, color='black', lw=0.8)
    ax2.set_xscale('log')
    ax2.set_xlabel('Redshift $z$')
    ax2.set_ylabel('Residuals (mag)')
    ax2.set_ylim(-0.8, 0.8)
    ax2.legend(loc='upper left', ncol=2)
    ax2.grid(True, alpha=0.3)

    plt.savefig('publications/article/figures/fig2_hubble_pantheon.pdf')
    plt.savefig('publications/article/figures/fig2_hubble_pantheon.png', dpi=150)
    plt.close()
    print("Figure 2: fig2_hubble_pantheon.pdf")


def fig3_q0_redshift(pantheon):
    """Figure 3: Evolution de q0 avec le redshift"""
    z, mu, sigma = pantheon['z'], pantheon['mu'], pantheon['sigma']
    janus = JanusCosmology()

    # Sous-echantillons
    z_cuts = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.3]
    results = []

    for i in range(len(z_cuts)-1):
        mask = (z >= z_cuts[i]) & (z < z_cuts[i+1])
        if np.sum(mask) > 30:
            z_sub, mu_sub, sig_sub = z[mask], mu[mask], sigma[mask]

            def chi2(p):
                q0, off = p
                if q0 >= 0 or q0 < -0.5: return 1e10
                mu_th = janus.distance_modulus(z_sub, q0) + off
                if np.any(~np.isfinite(mu_th)): return 1e10
                return np.sum(((mu_sub - mu_th) / sig_sub)**2)

            res = minimize(chi2, [-0.1, 0.0], method='Nelder-Mead')
            z_center = (z_cuts[i] + z_cuts[i+1]) / 2
            results.append({'z': z_center, 'q0': res.x[0], 'n': np.sum(mask)})

    fig, ax = plt.subplots(figsize=(8, 5))

    zs = [r['z'] for r in results]
    q0s = [r['q0'] for r in results]

    ax.plot(zs, q0s, 'bo-', ms=8, lw=2, label='JANUS $q_0(z)$')
    ax.axhline(-0.087, color='green', ls='--', lw=1.5, label='D\'Agostini & Petit 2018')
    ax.axhline(0, color='gray', ls=':', lw=1)

    ax.set_xlabel('Redshift $z$')
    ax.set_ylabel('Deceleration parameter $q_0$')
    ax.set_title('Evolution of $q_0$ with Redshift (Pantheon+)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.35, 0.05)

    plt.savefig('publications/article/figures/fig3_q0_evolution.pdf')
    plt.savefig('publications/article/figures/fig3_q0_evolution.png', dpi=150)
    plt.close()
    print("Figure 3: fig3_q0_evolution.pdf")


def fig4_comparison_bar(jla_j, jla_l, pan_j, pan_l):
    """Figure 4: Comparaison chi2/dof"""
    fig, ax = plt.subplots(figsize=(8, 5))

    x = np.arange(2)
    width = 0.35

    janus_vals = [jla_j['chi2']/jla_j['dof'], pan_j['chi2']/pan_j['dof']]
    lcdm_vals = [jla_l['chi2']/jla_l['dof'], pan_l['chi2']/pan_l['dof']]

    bars1 = ax.bar(x - width/2, janus_vals, width, label='JANUS', color='steelblue')
    bars2 = ax.bar(x + width/2, lcdm_vals, width, label=r'$\Lambda$CDM', color='coral')

    ax.axhline(1.0, color='black', ls='--', lw=1, alpha=0.5)
    ax.set_ylabel(r'$\chi^2$/dof')
    ax.set_title('Model Comparison: JANUS vs $\Lambda$CDM')
    ax.set_xticks(x)
    ax.set_xticklabels(['JLA (740 SNe)', 'Pantheon+ (1543 SNe)'])
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Annotations
    for bar, val in zip(bars1, janus_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{val:.3f}', ha='center', fontsize=9)
    for bar, val in zip(bars2, lcdm_vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{val:.3f}', ha='center', fontsize=9)

    plt.savefig('publications/article/figures/fig4_comparison.pdf')
    plt.savefig('publications/article/figures/fig4_comparison.png', dpi=150)
    plt.close()
    print("Figure 4: fig4_comparison.pdf")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("\n" + "=" * 60)
    print("GENERATION DES FIGURES POUR PUBLICATION")
    print("=" * 60)

    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    os.makedirs('publications/article/figures', exist_ok=True)

    # Charger donnees
    print("\nChargement des donnees...")
    jla, pantheon = load_data()
    print(f"  JLA: {jla['n']} SNe, Pantheon+: {pantheon['n']} SNe")

    # Ajustements
    print("\nAjustements...")
    jla_j, jla_l, janus_jla, lcdm_jla = fit_models(jla['z'], jla['mu'], jla['sigma'])
    pan_j, pan_l, janus_pan, lcdm_pan = fit_models(pantheon['z'], pantheon['mu'], pantheon['sigma'])

    print(f"  JLA JANUS: q0={jla_j['q0']:.4f}, chi2/dof={jla_j['chi2']/jla_j['dof']:.4f}")
    print(f"  JLA LCDM: chi2/dof={jla_l['chi2']/jla_l['dof']:.4f}")
    print(f"  Pantheon+ JANUS: q0={pan_j['q0']:.4f}, chi2/dof={pan_j['chi2']/pan_j['dof']:.4f}")
    print(f"  Pantheon+ LCDM: chi2/dof={pan_l['chi2']/pan_l['dof']:.4f}")

    # Figures
    print("\nGeneration des figures...")
    fig1_hubble_jla(jla, jla_j, jla_l, janus_jla, lcdm_jla)
    fig2_hubble_pantheon(pantheon, pan_j, pan_l, janus_pan, lcdm_pan)
    fig3_q0_redshift(pantheon)
    fig4_comparison_bar(jla_j, jla_l, pan_j, pan_l)

    print("\n" + "=" * 60)
    print("FIGURES GENEREES DANS publications/article/figures/")
    print("=" * 60)


if __name__ == "__main__":
    main()
