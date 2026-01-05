#!/usr/bin/env python3
"""
mcmc_diagnostics.py
Diagnostics MCMC pour l'analyse V2 JANUS-S

Vérifie:
1. Autocorrélation et temps de corrélation
2. Effective Sample Size (ESS)
3. Gelman-Rubin R-hat (convergence inter-chaînes)
4. Trace plots visuels
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pickle
import warnings

# Style
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.family': 'serif',
})

base_dir = Path(__file__).parent.parent
checkpoints_dir = base_dir / 'results' / 'v2_analysis' / 'checkpoints'
figures_dir = base_dir / 'results' / 'v2_analysis' / 'figures'
figures_dir.mkdir(parents=True, exist_ok=True)

# =============================================================================
# FONCTIONS DIAGNOSTICS
# =============================================================================

def autocorr_function(x, max_lag=None):
    """Calcule la fonction d'autocorrélation."""
    n = len(x)
    if max_lag is None:
        max_lag = n // 2
    x = x - np.mean(x)
    result = np.correlate(x, x, mode='full')
    result = result[n-1:n+max_lag]
    return result / result[0]

def integrated_autocorr_time(x, c=5.0):
    """
    Estime le temps d'autocorrélation intégré.
    Méthode de Goodman & Weare (2010).
    """
    n = len(x)
    max_lag = n // 2

    acf = autocorr_function(x, max_lag)

    # Trouver où l'autocorrélation devient négligeable
    taus = 2.0 * np.cumsum(acf) - 1.0

    # Critère d'arrêt: M < c * tau
    for M in range(1, max_lag):
        if M >= c * taus[M]:
            return taus[M]

    # Si pas de convergence, retourner estimation prudente
    return taus[max_lag-1]

def effective_sample_size(chain, tau=None):
    """Calcule le nombre effectif d'échantillons."""
    n_steps, n_walkers = chain.shape[:2]
    n_total = n_steps * n_walkers

    if tau is None:
        # Estimer tau sur la moyenne des walkers
        mean_chain = np.mean(chain, axis=1)
        tau = integrated_autocorr_time(mean_chain)

    return n_total / tau

def gelman_rubin(chains):
    """
    Calcule le diagnostic R-hat de Gelman-Rubin.

    chains: array de shape (n_steps, n_chains) ou (n_steps, n_chains, n_params)
    Retourne R-hat pour chaque paramètre.
    """
    if chains.ndim == 2:
        chains = chains[:, :, np.newaxis]

    n_steps, n_chains, n_params = chains.shape

    # Utiliser la seconde moitié des chaînes
    chains = chains[n_steps//2:, :, :]
    n = chains.shape[0]

    rhats = []
    for p in range(n_params):
        # Moyenne de chaque chaîne
        chain_means = np.mean(chains[:, :, p], axis=0)
        # Moyenne globale
        global_mean = np.mean(chain_means)

        # Variance inter-chaînes (B)
        B = n * np.var(chain_means, ddof=1)

        # Variance intra-chaîne (W)
        chain_vars = np.var(chains[:, :, p], axis=0, ddof=1)
        W = np.mean(chain_vars)

        # Estimateur de variance
        var_est = (1 - 1/n) * W + (1/n) * B

        # R-hat
        rhat = np.sqrt(var_est / W) if W > 0 else np.nan
        rhats.append(rhat)

    return np.array(rhats)

# =============================================================================
# CHARGEMENT DES CHECKPOINTS
# =============================================================================

print("=" * 70)
print("DIAGNOSTICS MCMC - ANALYSE V2")
print("=" * 70)

checkpoint_files = sorted(checkpoints_dir.glob('checkpoint_*.pkl'))
print(f"\nCheckpoints trouvés: {len(checkpoint_files)}")

all_diagnostics = {}

for ckpt_file in checkpoint_files:
    with open(ckpt_file, 'rb') as f:
        ckpt = pickle.load(f)

    name = f"{ckpt['model']}_{ckpt['dataset']}"
    print(f"\n{'='*50}")
    print(f"Analyse: {name}")
    print(f"{'='*50}")

    chains = ckpt['chains']  # (n_steps, n_walkers, n_params)
    log_prob = ckpt['log_prob']  # (n_steps, n_walkers)
    param_names = ckpt['param_names']
    n_steps, n_walkers, n_params = chains.shape

    print(f"  Steps: {n_steps}, Walkers: {n_walkers}, Params: {n_params}")
    print(f"  Paramètres: {param_names}")

    diagnostics = {'name': name, 'model': ckpt['model'], 'dataset': ckpt['dataset']}

    # 1. Autocorrélation et ESS par paramètre
    print(f"\n  Autocorrélation et ESS:")
    taus = []
    ess_values = []

    for p in range(n_params):
        # Moyenne sur les walkers pour estimer tau
        mean_chain = np.mean(chains[:, :, p], axis=1)
        try:
            tau = integrated_autocorr_time(mean_chain)
        except:
            tau = n_steps  # Valeur prudente si échec

        ess = (n_steps * n_walkers) / tau
        taus.append(tau)
        ess_values.append(ess)

        print(f"    {param_names[p]}: tau={tau:.1f}, ESS={ess:.0f}")

    diagnostics['tau'] = taus
    diagnostics['ess'] = ess_values

    # 2. Gelman-Rubin R-hat (walkers comme chaînes séparées)
    print(f"\n  Gelman-Rubin R-hat:")
    rhats = gelman_rubin(chains)
    for p, rhat in enumerate(rhats):
        status = "OK" if rhat < 1.1 else "ATTENTION" if rhat < 1.2 else "PROBLEME"
        print(f"    {param_names[p]}: R-hat={rhat:.4f} [{status}]")

    diagnostics['rhat'] = rhats.tolist()

    # 3. Acceptance rate approximatif (changements dans log_prob)
    n_changes = np.sum(np.diff(log_prob, axis=0) != 0)
    acceptance = n_changes / ((n_steps - 1) * n_walkers)
    print(f"\n  Taux d'acceptation approx: {acceptance:.1%}")
    diagnostics['acceptance'] = acceptance

    # 4. Log-prob final
    final_log_prob = np.mean(log_prob[-100:, :])
    print(f"  Log-prob moyen (100 derniers): {final_log_prob:.2f}")
    diagnostics['final_log_prob'] = final_log_prob

    all_diagnostics[name] = diagnostics

# =============================================================================
# FIGURE 1: TRACE PLOTS
# =============================================================================
print("\n" + "=" * 70)
print("Génération des figures...")
print("=" * 70)

fig, axes = plt.subplots(4, 3, figsize=(15, 12))
fig.suptitle('MCMC Trace Plots - V2 Analysis', fontsize=14, fontweight='bold')

datasets = ['Pantheon+', 'DES-SN5YR', 'Combined']
models = ['janus', 'lcdm']
colors = {'janus': '#e74c3c', 'lcdm': '#3498db'}
param_labels = {'janus': ['q₀', 'offset'], 'lcdm': ['Ωₘ', 'offset']}

for col, dataset in enumerate(datasets):
    for row, model in enumerate(models):
        name = f"{model}_{dataset}"
        ckpt_file = checkpoints_dir / f"checkpoint_{model}_{dataset.lower().replace('-', '_').replace('+', '_')}.pkl"

        if not ckpt_file.exists():
            continue

        with open(ckpt_file, 'rb') as f:
            ckpt = pickle.load(f)

        chains = ckpt['chains']

        # Plot premier paramètre (q0 ou Omega_m)
        ax = axes[row*2, col]
        for w in range(min(8, chains.shape[1])):  # 8 walkers pour lisibilité
            ax.plot(chains[:, w, 0], alpha=0.5, lw=0.5, color=colors[model])
        ax.set_ylabel(param_labels[model][0])
        if row == 0:
            ax.set_title(dataset)
        ax.axhline(np.mean(chains[-500:, :, 0]), color='black', ls='--', lw=1)

        # Plot offset
        ax = axes[row*2 + 1, col]
        for w in range(min(8, chains.shape[1])):
            ax.plot(chains[:, w, 1], alpha=0.5, lw=0.5, color=colors[model])
        ax.set_ylabel('offset')
        if row == 1:
            ax.set_xlabel('Step')
        ax.axhline(np.mean(chains[-500:, :, 1]), color='black', ls='--', lw=1)

# Add model labels
axes[0, 0].annotate('JANUS', xy=(-0.3, 0.5), xycoords='axes fraction',
                    fontsize=12, fontweight='bold', rotation=90, va='center')
axes[2, 0].annotate('ΛCDM', xy=(-0.3, 0.5), xycoords='axes fraction',
                    fontsize=12, fontweight='bold', rotation=90, va='center')

plt.tight_layout()
plt.savefig(figures_dir / 'mcmc_trace_plots.png', dpi=300)
plt.savefig(figures_dir / 'mcmc_trace_plots.pdf')
plt.close()
print("  ✓ mcmc_trace_plots.png/pdf")

# =============================================================================
# FIGURE 2: AUTOCORRÉLATION
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(12, 6))
fig.suptitle('Autocorrelation Functions', fontsize=14, fontweight='bold')

for col, dataset in enumerate(datasets):
    for row, model in enumerate(models):
        ckpt_file = checkpoints_dir / f"checkpoint_{model}_{dataset.lower().replace('-', '_').replace('+', '_')}.pkl"

        if not ckpt_file.exists():
            continue

        with open(ckpt_file, 'rb') as f:
            ckpt = pickle.load(f)

        chains = ckpt['chains']
        mean_chain = np.mean(chains[:, :, 0], axis=1)  # Paramètre principal

        ax = axes[row, col]
        max_lag = min(200, len(mean_chain) // 4)
        acf = autocorr_function(mean_chain, max_lag)

        ax.plot(acf, color=colors[model], lw=1.5)
        ax.axhline(0, color='gray', ls='--', lw=0.5)
        ax.axhline(0.1, color='gray', ls=':', lw=0.5)
        ax.axhline(-0.1, color='gray', ls=':', lw=0.5)

        ax.set_ylabel('ACF' if col == 0 else '')
        ax.set_xlabel('Lag' if row == 1 else '')
        if row == 0:
            ax.set_title(dataset)
        ax.set_ylim(-0.3, 1.1)
        ax.text(0.95, 0.95, f'{model.upper()}', transform=ax.transAxes,
                ha='right', va='top', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig(figures_dir / 'mcmc_autocorr.png', dpi=300)
plt.savefig(figures_dir / 'mcmc_autocorr.pdf')
plt.close()
print("  ✓ mcmc_autocorr.png/pdf")

# =============================================================================
# FIGURE 3: R-HAT SUMMARY
# =============================================================================
fig, ax = plt.subplots(figsize=(10, 5))

x_labels = []
rhat_values = []
colors_bar = []

for dataset in datasets:
    for model in models:
        name = f"{model}_{dataset}"
        if name in all_diagnostics:
            for p, rhat in enumerate(all_diagnostics[name]['rhat']):
                param = param_labels[model][p]
                x_labels.append(f"{dataset}\n{model}\n{param}")
                rhat_values.append(rhat)
                colors_bar.append(colors[model])

x = np.arange(len(x_labels))
bars = ax.bar(x, rhat_values, color=colors_bar, alpha=0.7, edgecolor='black')

ax.axhline(1.0, color='green', ls='-', lw=2, label='Perfect (R-hat=1)')
ax.axhline(1.1, color='orange', ls='--', lw=1.5, label='Acceptable (R-hat<1.1)')
ax.axhline(1.2, color='red', ls='--', lw=1.5, label='Warning (R-hat>1.2)')

ax.set_xticks(x)
ax.set_xticklabels(x_labels, fontsize=7)
ax.set_ylabel('Gelman-Rubin R-hat')
ax.set_title('Convergence Diagnostic: R-hat by Parameter')
ax.legend(loc='upper right')
ax.set_ylim(0.95, 1.25)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(figures_dir / 'mcmc_rhat.png', dpi=300)
plt.savefig(figures_dir / 'mcmc_rhat.pdf')
plt.close()
print("  ✓ mcmc_rhat.png/pdf")

# =============================================================================
# RAPPORT FINAL
# =============================================================================
print("\n" + "=" * 70)
print("RAPPORT DIAGNOSTICS MCMC")
print("=" * 70)

print("\n┌" + "─" * 68 + "┐")
print("│{:^68}│".format("RÉSUMÉ DES DIAGNOSTICS"))
print("├" + "─" * 20 + "┬" + "─" * 15 + "┬" + "─" * 15 + "┬" + "─" * 15 + "┤")
print("│{:^20}│{:^15}│{:^15}│{:^15}│".format("Modèle/Dataset", "R-hat max", "ESS min", "Status"))
print("├" + "─" * 20 + "┼" + "─" * 15 + "┼" + "─" * 15 + "┼" + "─" * 15 + "┤")

all_ok = True
for name, diag in all_diagnostics.items():
    rhat_max = max(diag['rhat'])
    ess_min = min(diag['ess'])

    if rhat_max < 1.1 and ess_min > 100:
        status = "✓ OK"
    elif rhat_max < 1.2 and ess_min > 50:
        status = "~ Acceptable"
    else:
        status = "✗ Check"
        all_ok = False

    short_name = f"{diag['model']}/{diag['dataset'][:8]}"
    print("│{:^20}│{:^15.4f}│{:^15.0f}│{:^15}│".format(short_name, rhat_max, ess_min, status))

print("└" + "─" * 20 + "┴" + "─" * 15 + "┴" + "─" * 15 + "┴" + "─" * 15 + "┘")

print("\n" + "=" * 70)
if all_ok:
    print("CONCLUSION: Toutes les chaînes ont convergé correctement ✓")
else:
    print("CONCLUSION: Certaines chaînes nécessitent attention")
print("=" * 70)

print(f"""
Figures générées:
  - mcmc_trace_plots.png/pdf : Évolution des paramètres
  - mcmc_autocorr.png/pdf    : Fonctions d'autocorrélation
  - mcmc_rhat.png/pdf        : Diagnostic R-hat

Critères utilisés:
  - R-hat < 1.1 : Convergence acceptable
  - ESS > 100   : Échantillonnage suffisant
""")
