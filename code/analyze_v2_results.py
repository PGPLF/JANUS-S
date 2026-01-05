#!/usr/bin/env python3
"""
analyze_v2_results.py
Analyse approfondie des résultats V2 JANUS vs LCDM

Questions clés:
1. Pourquoi q0 varie entre datasets (-0.013 à -0.096)?
2. D'où vient la tension dans le dataset Combined?
3. Comment se compare q0=-0.096 avec la référence 2018 (-0.087)?
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json
import sys

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
print("ANALYSE DES RÉSULTATS V2")
print("=" * 70)

# =============================================================================
# 1. ANALYSE DE LA VARIATION DE q0
# =============================================================================
print("\n" + "=" * 70)
print("1. VARIATION DE q0 ENTRE DATASETS")
print("=" * 70)

q0_pantheon = results['Pantheon+']['janus']['q0']
q0_des = results['DES-SN5YR']['janus']['q0']
q0_combined = results['Combined']['janus']['q0']

q0_err_pantheon = results['Pantheon+']['janus']['q0_err']
q0_err_des = results['DES-SN5YR']['janus']['q0_err']
q0_err_combined = results['Combined']['janus']['q0_err']

print(f"\nValeurs de q0:")
print(f"  Pantheon+:  {q0_pantheon:.4f} ± {q0_err_pantheon:.4f}")
print(f"  DES-SN5YR:  {q0_des:.4f} ± {q0_err_des:.4f}")
print(f"  Combined:   {q0_combined:.4f} ± {q0_err_combined:.4f}")
print(f"  Ref 2018:   -0.087")

# Tension entre datasets individuels
tension_p_d = abs(q0_pantheon - q0_des) / np.sqrt(q0_err_pantheon**2 + q0_err_des**2)
print(f"\nTension Pantheon+ vs DES: {tension_p_d:.1f}σ")

# Le Combined devrait être entre les deux, mais il est très différent!
q0_weighted_avg = (q0_pantheon/q0_err_pantheon**2 + q0_des/q0_err_des**2) / \
                  (1/q0_err_pantheon**2 + 1/q0_err_des**2)
print(f"\nMoyenne pondérée attendue: {q0_weighted_avg:.4f}")
print(f"Valeur Combined obtenue:   {q0_combined:.4f}")
print(f"Différence: {abs(q0_combined - q0_weighted_avg):.4f}")

# =============================================================================
# 2. ANALYSE DES DISTRIBUTIONS EN REDSHIFT
# =============================================================================
print("\n" + "=" * 70)
print("2. DISTRIBUTIONS EN REDSHIFT")
print("=" * 70)

# Load data
pantheon = load_pantheon_with_covariance(
    str(data_dir / 'pantheon' / 'Pantheon+SH0ES.dat'),
    str(data_dir / 'pantheon' / 'Pantheon+SH0ES_STAT+SYS.cov')
)
des = load_des_sn5yr(
    str(data_dir / 'des_sn5yr' / 'DES-Dovekie_HD.csv'),
    str(data_dir / 'des_sn5yr' / 'STAT+SYS.npz')
)
combined = combine_datasets(pantheon, des, method='priority_pantheon')

print(f"\nPantheon+:")
print(f"  N = {pantheon['n_sne']}")
print(f"  z: [{pantheon['z'].min():.4f}, {pantheon['z'].max():.4f}]")
print(f"  z median: {np.median(pantheon['z']):.4f}")
print(f"  z mean: {np.mean(pantheon['z']):.4f}")

print(f"\nDES-SN5YR:")
print(f"  N = {des['n_sne']}")
print(f"  z: [{des['z'].min():.4f}, {des['z'].max():.4f}]")
print(f"  z median: {np.median(des['z']):.4f}")
print(f"  z mean: {np.mean(des['z']):.4f}")

print(f"\nCombined:")
print(f"  N = {combined['n_sne']}")
print(f"  z: [{combined['z'].min():.4f}, {combined['z'].max():.4f}]")
print(f"  z median: {np.median(combined['z']):.4f}")
print(f"  z mean: {np.mean(combined['z']):.4f}")

# =============================================================================
# 3. FIGURE: DISTRIBUTIONS EN REDSHIFT
# =============================================================================
fig, axes = plt.subplots(1, 3, figsize=(14, 4))

axes[0].hist(pantheon['z'], bins=50, alpha=0.7, color='blue', edgecolor='black')
axes[0].axvline(np.median(pantheon['z']), color='red', ls='--', label=f'median={np.median(pantheon["z"]):.3f}')
axes[0].set_xlabel('Redshift z')
axes[0].set_ylabel('Count')
axes[0].set_title(f'Pantheon+ (N={pantheon["n_sne"]})')
axes[0].legend()

axes[1].hist(des['z'], bins=50, alpha=0.7, color='green', edgecolor='black')
axes[1].axvline(np.median(des['z']), color='red', ls='--', label=f'median={np.median(des["z"]):.3f}')
axes[1].set_xlabel('Redshift z')
axes[1].set_title(f'DES-SN5YR (N={des["n_sne"]})')
axes[1].legend()

axes[2].hist(combined['z'], bins=50, alpha=0.7, color='purple', edgecolor='black')
axes[2].axvline(np.median(combined['z']), color='red', ls='--', label=f'median={np.median(combined["z"]):.3f}')
axes[2].set_xlabel('Redshift z')
axes[2].set_title(f'Combined (N={combined["n_sne"]})')
axes[2].legend()

plt.tight_layout()
plt.savefig(figures_dir / 'redshift_distributions.png', dpi=150)
plt.close()
print(f"\nSaved: {figures_dir / 'redshift_distributions.png'}")

# =============================================================================
# 4. ANALYSE: POURQUOI LE COMBINED EST SI DIFFÉRENT?
# =============================================================================
print("\n" + "=" * 70)
print("3. ANALYSE DE LA TENSION DANS COMBINED")
print("=" * 70)

# Hypothèse: Les deux datasets ont des calibrations différentes
# Le "offset" absorbe cette différence dans les fits individuels
# Mais quand on combine, la matrice de covariance ne capture pas cette tension

offset_p_janus = results['Pantheon+']['janus']['offset']
offset_d_janus = results['DES-SN5YR']['janus']['offset']
offset_c_janus = results['Combined']['janus']['offset']

offset_p_lcdm = results['Pantheon+']['lcdm']['offset']
offset_d_lcdm = results['DES-SN5YR']['lcdm']['offset']
offset_c_lcdm = results['Combined']['lcdm']['offset']

print(f"\nOffsets JANUS:")
print(f"  Pantheon+: {offset_p_janus:.4f}")
print(f"  DES-SN5YR: {offset_d_janus:.4f}")
print(f"  Combined:  {offset_c_janus:.4f}")
print(f"  Différence P+ vs DES: {offset_p_janus - offset_d_janus:.4f}")

print(f"\nOffsets LCDM:")
print(f"  Pantheon+: {offset_p_lcdm:.4f}")
print(f"  DES-SN5YR: {offset_d_lcdm:.4f}")
print(f"  Combined:  {offset_c_lcdm:.4f}")
print(f"  Différence P+ vs DES: {offset_p_lcdm - offset_d_lcdm:.4f}")

# La grande différence d'offset (~0.09-0.12 mag) entre les datasets
# indique une tension de calibration!

print("\n" + "-" * 50)
print("DIAGNOSTIC: Différence d'offset ~0.09-0.12 mag")
print("Ceci indique une TENSION DE CALIBRATION entre datasets")
print("-" * 50)

# =============================================================================
# 5. CHI2 PAR DATASET DANS LE FIT COMBINED
# =============================================================================
print("\n" + "=" * 70)
print("4. CHI² COMPARAISON")
print("=" * 70)

chi2_p_janus = results['Pantheon+']['janus']['chi2']
chi2_d_janus = results['DES-SN5YR']['janus']['chi2']
chi2_c_janus = results['Combined']['janus']['chi2']

chi2_p_lcdm = results['Pantheon+']['lcdm']['chi2']
chi2_d_lcdm = results['DES-SN5YR']['lcdm']['chi2']
chi2_c_lcdm = results['Combined']['lcdm']['chi2']

print(f"\nχ² JANUS:")
print(f"  Pantheon+: {chi2_p_janus:.1f} (dof={1701-2}={1699})")
print(f"  DES-SN5YR: {chi2_d_janus:.1f} (dof={1820-2}={1818})")
print(f"  Combined:  {chi2_c_janus:.1f} (dof={3182-2}={3180})")
print(f"  Sum individuel: {chi2_p_janus + chi2_d_janus:.1f}")
print(f"  Excès Combined: {chi2_c_janus - (chi2_p_janus + chi2_d_janus):.1f}")

print(f"\nχ² LCDM:")
print(f"  Pantheon+: {chi2_p_lcdm:.1f}")
print(f"  DES-SN5YR: {chi2_d_lcdm:.1f}")
print(f"  Combined:  {chi2_c_lcdm:.1f}")
print(f"  Sum individuel: {chi2_p_lcdm + chi2_d_lcdm:.1f}")
print(f"  Excès Combined: {chi2_c_lcdm - (chi2_p_lcdm + chi2_d_lcdm):.1f}")

# =============================================================================
# 6. FIGURE: COMPARAISON q0 AVEC RÉFÉRENCE
# =============================================================================
fig, ax = plt.subplots(figsize=(8, 6))

datasets = ['Pantheon+', 'DES-SN5YR', 'Combined', 'Ref 2018']
q0_values = [q0_pantheon, q0_des, q0_combined, -0.087]
q0_errors = [q0_err_pantheon, q0_err_des, q0_err_combined, 0.01]  # assuming 0.01 for ref
colors = ['blue', 'green', 'purple', 'red']

x = np.arange(len(datasets))
ax.errorbar(x, q0_values, yerr=q0_errors, fmt='o', markersize=10,
            capsize=5, capthick=2, color='black')

for i, (xi, yi, c) in enumerate(zip(x, q0_values, colors)):
    ax.scatter(xi, yi, s=150, c=c, zorder=5, edgecolor='black', linewidth=2)

ax.axhline(0, color='gray', ls=':', alpha=0.5)
ax.axhline(-0.087, color='red', ls='--', alpha=0.5, label='Ref 2018: q₀=-0.087')

ax.set_xticks(x)
ax.set_xticklabels(datasets)
ax.set_ylabel('q₀ (deceleration parameter)', fontsize=12)
ax.set_title('JANUS q₀ across datasets', fontsize=14)
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(figures_dir / 'q0_comparison.png', dpi=150)
plt.close()
print(f"\nSaved: {figures_dir / 'q0_comparison.png'}")

# =============================================================================
# 7. CONCLUSIONS
# =============================================================================
print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)

print("""
1. VARIATION DE q₀:
   - Pantheon+ et DES-SN5YR donnent des q₀ proches de 0 (~-0.01 à -0.03)
   - Le Combined donne q₀ ≈ -0.096, proche de la référence 2018 (-0.087)
   - Cette variation suggère une sensibilité au coverage en redshift

2. TENSION DE CALIBRATION:
   - Grande différence d'offset entre Pantheon+ et DES (~0.09-0.12 mag)
   - Ceci indique une tension systématique entre les deux surveys
   - Le fit Combined doit "moyenner" ces différences

3. INTERPRÉTATION:
   - q₀ proche de 0 (Pantheon+, DES seuls) → expansion quasi-linéaire
   - q₀ ≈ -0.1 (Combined, Ref 2018) → décélération modérée
   - La valeur dépend fortement de la combinaison des données

4. RECOMMANDATIONS:
   - Investiguer la source de la tension de calibration
   - Ajouter un paramètre d'offset relatif entre surveys
   - Considérer une analyse avec matrices de covariance croisées
""")

print("=" * 70)
print("Analyse terminée")
print("=" * 70)
