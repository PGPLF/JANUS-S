#!/usr/bin/env python3
"""
generate_v2_figures.py
Figures publication pour l'analyse V2 JANUS-S (Pantheon+ + DES-SN5YR)

Figures:
1. Hubble diagram triple panel (Pantheon+, DES, Combined)
2. Parameter comparison (q0, Omega_m, ΔAIC)
3. q0 vs Omega_m synthesis
4. Chi2/dof comparison
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import json
import sys

# Style publication
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.family': 'serif',
})

sys.path.insert(0, str(Path(__file__).parent))
from data_loader import load_pantheon_with_covariance, load_des_sn5yr, combine_datasets
from janus_model import JanusCosmology, LambdaCDM

# Paths
base_dir = Path(__file__).parent.parent
data_dir = base_dir / 'data'
results_dir = base_dir / 'results' / 'v2_analysis'
figures_dir = results_dir / 'figures'
figures_dir.mkdir(parents=True, exist_ok=True)

# Load results
with open(results_dir / 'v2_results.json', 'r') as f:
    results = json.load(f)

print("=" * 70)
print("GÉNÉRATION DES FIGURES PUBLICATION V2")
print("=" * 70)

# Load data
print("\nChargement des données...")
pantheon = load_pantheon_with_covariance(
    str(data_dir / 'pantheon' / 'Pantheon+SH0ES.dat'),
    str(data_dir / 'pantheon' / 'Pantheon+SH0ES_STAT+SYS.cov')
)
des = load_des_sn5yr(
    str(data_dir / 'des_sn5yr' / 'DES-Dovekie_HD.csv'),
    str(data_dir / 'des_sn5yr' / 'STAT+SYS.npz')
)
combined = combine_datasets(pantheon, des, method='priority_pantheon')

# Models
janus = JanusCosmology(H0=70)
lcdm = LambdaCDM(H0=70)

# Colors
color_janus = '#e74c3c'
color_lcdm = '#3498db'
colors_data = ['#1f77b4', '#2ca02c', '#9467bd']

# =============================================================================
# FIGURE 1: HUBBLE DIAGRAM TRIPLE PANEL
# =============================================================================
print("\n[1/4] Hubble Diagram triple panel...")

fig = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 3, height_ratios=[3, 1], hspace=0.05, wspace=0.25)

datasets = [
    ('Pantheon+', pantheon, 'Pantheon+'),
    ('DES-SN5YR', des, 'DES-SN5YR'),
    ('Combined', combined, 'Combined')
]

for idx, (name, data, key) in enumerate(datasets):
    ax_main = fig.add_subplot(gs[0, idx])
    ax_res = fig.add_subplot(gs[1, idx], sharex=ax_main)

    z = data['z']
    mu = data['mu']
    sigma = data['sigma_mu']

    # Get best-fit parameters
    q0 = results[key]['janus']['q0']
    offset_j = results[key]['janus']['offset']
    Om = results[key]['lcdm']['Omega_m']
    offset_l = results[key]['lcdm']['offset']

    # Model curves
    z_model = np.linspace(0.001, z.max() * 1.05, 500)
    mu_janus = janus.distance_modulus(z_model, q0=q0) + offset_j
    mu_lcdm = lcdm.distance_modulus(z_model, Omega_m=Om) + offset_l

    # Subsample for clarity
    np.random.seed(42)
    if len(z) > 500:
        idx_sub = np.random.choice(len(z), 500, replace=False)
        z_plot, mu_plot, sigma_plot = z[idx_sub], mu[idx_sub], sigma[idx_sub]
    else:
        z_plot, mu_plot, sigma_plot = z, mu, sigma

    # Main panel
    ax_main.errorbar(z_plot, mu_plot, yerr=sigma_plot, fmt='o', ms=2, alpha=0.3,
                     color=colors_data[idx], elinewidth=0.5, label='Data')
    ax_main.plot(z_model, mu_janus, '-', color=color_janus, lw=2,
                 label=f'JANUS (q₀={q0:.3f})')
    ax_main.plot(z_model, mu_lcdm, '--', color=color_lcdm, lw=2,
                 label=f'ΛCDM (Ωₘ={Om:.3f})')

    ax_main.set_ylabel('Distance modulus μ' if idx == 0 else '')
    ax_main.set_title(f'{name} (N={data["n_sne"]})')
    ax_main.legend(loc='lower right', fontsize=8)
    ax_main.set_xlim(0, z.max() * 1.05)
    plt.setp(ax_main.get_xticklabels(), visible=False)

    # Residuals panel
    mu_janus_data = janus.distance_modulus(z_plot, q0=q0) + offset_j
    mu_lcdm_data = lcdm.distance_modulus(z_plot, Omega_m=Om) + offset_l

    res_j = mu_plot - mu_janus_data
    res_l = mu_plot - mu_lcdm_data

    ax_res.scatter(z_plot, res_j, s=3, alpha=0.4, c=color_janus, label='JANUS')
    ax_res.scatter(z_plot, res_l, s=3, alpha=0.4, c=color_lcdm, label='ΛCDM')
    ax_res.axhline(0, color='black', ls='-', lw=0.5)
    ax_res.axhline(0.5, color='gray', ls=':', lw=0.5)
    ax_res.axhline(-0.5, color='gray', ls=':', lw=0.5)

    ax_res.set_xlabel('Redshift z')
    ax_res.set_ylabel('Residuals' if idx == 0 else '')
    ax_res.set_ylim(-1.5, 1.5)
    ax_res.set_xlim(0, z.max() * 1.05)
    if idx == 2:
        ax_res.legend(loc='upper right', fontsize=8)

plt.savefig(figures_dir / 'fig1_hubble_v2.png', dpi=300)
plt.savefig(figures_dir / 'fig1_hubble_v2.pdf')
plt.close()
print(f"  ✓ fig1_hubble_v2.png/pdf")

# =============================================================================
# FIGURE 2: PARAMETER COMPARISON
# =============================================================================
print("\n[2/4] Parameter comparison...")

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

datasets_names = ['Pantheon+', 'DES-SN5YR', 'Combined']
x = np.arange(3)

# Panel 1: q0 values
ax = axes[0]
q0_vals = [results[d]['janus']['q0'] for d in datasets_names]
q0_errs = [results[d]['janus']['q0_err'] for d in datasets_names]

ax.errorbar(x, q0_vals, yerr=q0_errs, fmt='o', markersize=12, capsize=6,
            capthick=2, color='black', zorder=3)
for i, (xi, yi, c) in enumerate(zip(x, q0_vals, colors_data)):
    ax.scatter(xi, yi, s=150, c=c, zorder=4, edgecolor='black', linewidth=1.5)

ax.axhline(-0.087, color='red', ls='--', lw=1.5, label='Ref 2018: -0.087')
ax.axhline(0, color='gray', ls=':', lw=1)
ax.set_xticks(x)
ax.set_xticklabels(datasets_names, rotation=15)
ax.set_ylabel('q₀ (JANUS)')
ax.set_title('Deceleration parameter')
ax.legend(loc='lower left')
ax.grid(True, alpha=0.3)

# Panel 2: Omega_m values
ax = axes[1]
Om_vals = [results[d]['lcdm']['Omega_m'] for d in datasets_names]
Om_errs = [results[d]['lcdm']['Omega_m_err'] for d in datasets_names]

ax.errorbar(x, Om_vals, yerr=Om_errs, fmt='o', markersize=12, capsize=6,
            capthick=2, color='black', zorder=3)
for i, (xi, yi, c) in enumerate(zip(x, Om_vals, colors_data)):
    ax.scatter(xi, yi, s=150, c=c, zorder=4, edgecolor='black', linewidth=1.5)

ax.axhline(0.315, color='red', ls='--', lw=1.5, label='Planck 2018: 0.315')
ax.set_xticks(x)
ax.set_xticklabels(datasets_names, rotation=15)
ax.set_ylabel('Ωₘ (ΛCDM)')
ax.set_title('Matter density')
ax.legend(loc='upper right')
ax.grid(True, alpha=0.3)

# Panel 3: Delta AIC
ax = axes[2]
daic_vals = [results[d]['comparison']['delta_aic'] for d in datasets_names]
colors_aic = ['#e74c3c' if d > 0 else '#27ae60' for d in daic_vals]

bars = ax.bar(x, daic_vals, color=colors_aic, alpha=0.7, edgecolor='black', linewidth=1.5)
ax.axhline(0, color='black', ls='-', lw=1)
ax.axhline(10, color='gray', ls='--', lw=1)
ax.axhline(-10, color='gray', ls='--', lw=1)

ax.set_xticks(x)
ax.set_xticklabels(datasets_names, rotation=15)
ax.set_ylabel('ΔAIC (JANUS - ΛCDM)')
ax.set_title('Model comparison')

ax.annotate('ΛCDM preferred ↑', xy=(0.95, 0.95), xycoords='axes fraction',
            ha='right', va='top', fontsize=9, color='#e74c3c')
ax.grid(True, alpha=0.3, axis='y')

# Add values on bars
for bar, val in zip(bars, daic_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
            f'{val:.0f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig(figures_dir / 'fig2_parameters_v2.png', dpi=300)
plt.savefig(figures_dir / 'fig2_parameters_v2.pdf')
plt.close()
print(f"  ✓ fig2_parameters_v2.png/pdf")

# =============================================================================
# FIGURE 3: CHI2/DOF COMPARISON
# =============================================================================
print("\n[3/4] Chi2/dof comparison...")

fig, ax = plt.subplots(figsize=(8, 5))

width = 0.35

chi2_janus = [results[d]['janus']['chi2_dof'] for d in datasets_names]
chi2_lcdm = [results[d]['lcdm']['chi2_dof'] for d in datasets_names]

bars1 = ax.bar(x - width/2, chi2_janus, width, label='JANUS', color=color_janus,
               alpha=0.8, edgecolor='black')
bars2 = ax.bar(x + width/2, chi2_lcdm, width, label='ΛCDM', color=color_lcdm,
               alpha=0.8, edgecolor='black')

ax.axhline(1.0, color='black', ls='--', lw=1.5, label='Perfect fit')
ax.axhline(1.1, color='gray', ls=':', lw=1)
ax.axhline(0.9, color='gray', ls=':', lw=1)

ax.set_xticks(x)
ax.set_xticklabels(datasets_names)
ax.set_ylabel('χ²/dof')
ax.set_title('Goodness of fit: JANUS vs ΛCDM')
ax.legend()
ax.set_ylim(0.8, 1.15)
ax.grid(True, alpha=0.3, axis='y')

for bar, val in zip(bars1, chi2_janus):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{val:.3f}', ha='center', va='bottom', fontsize=9)
for bar, val in zip(bars2, chi2_lcdm):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{val:.3f}', ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(figures_dir / 'fig3_chi2_v2.png', dpi=300)
plt.savefig(figures_dir / 'fig3_chi2_v2.pdf')
plt.close()
print(f"  ✓ fig3_chi2_v2.png/pdf")

# =============================================================================
# FIGURE 4: SYNTHESIS TABLE AS FIGURE
# =============================================================================
print("\n[4/4] Summary table figure...")

fig, ax = plt.subplots(figsize=(10, 4))
ax.axis('off')

# Table data
columns = ['Dataset', 'N SNe', 'JANUS q₀', 'χ²/dof', 'ΛCDM Ωₘ', 'χ²/dof', 'ΔAIC', 'Preferred']
rows = []
for d in datasets_names:
    q0 = results[d]['janus']['q0']
    q0_err = results[d]['janus']['q0_err']
    chi2_j = results[d]['janus']['chi2_dof']
    Om = results[d]['lcdm']['Omega_m']
    Om_err = results[d]['lcdm']['Omega_m_err']
    chi2_l = results[d]['lcdm']['chi2_dof']
    daic = results[d]['comparison']['delta_aic']
    pref = 'ΛCDM' if daic > 0 else 'JANUS'

    n_sne = results[d]['comparison']['n_sne'] if 'n_sne' in results[d]['comparison'] else '?'

    rows.append([
        d,
        str(n_sne),
        f'{q0:.4f} ± {q0_err:.4f}',
        f'{chi2_j:.3f}',
        f'{Om:.4f} ± {Om_err:.4f}',
        f'{chi2_l:.3f}',
        f'{daic:.1f}',
        pref
    ])

table = ax.table(cellText=rows, colLabels=columns, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.8)

# Style header
for i in range(len(columns)):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(color='white', fontweight='bold')

# Style rows
for i in range(1, len(rows) + 1):
    for j in range(len(columns)):
        if i % 2 == 0:
            table[(i, j)].set_facecolor('#D9E2F3')
        # Color the Preferred column
        if j == 7:
            if rows[i-1][7] == 'ΛCDM':
                table[(i, j)].set_text_props(color='#e74c3c', fontweight='bold')
            else:
                table[(i, j)].set_text_props(color='#27ae60', fontweight='bold')

ax.set_title('JANUS-S V2 Analysis Summary\n(Pantheon+ + DES-SN5YR)', fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig(figures_dir / 'fig4_summary_v2.png', dpi=300, bbox_inches='tight')
plt.savefig(figures_dir / 'fig4_summary_v2.pdf', bbox_inches='tight')
plt.close()
print(f"  ✓ fig4_summary_v2.png/pdf")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("FIGURES GÉNÉRÉES")
print("=" * 70)
print(f"""
Dossier: {figures_dir}

1. fig1_hubble_v2.png/pdf
   → Diagramme de Hubble triple panel avec résidus

2. fig2_parameters_v2.png/pdf
   → Comparaison q₀, Ωₘ, ΔAIC entre datasets

3. fig3_chi2_v2.png/pdf
   → Qualité d'ajustement χ²/dof

4. fig4_summary_v2.png/pdf
   → Tableau récapitulatif
""")
print("=" * 70)
