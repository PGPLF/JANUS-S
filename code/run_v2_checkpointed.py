#!/usr/bin/env python3
"""
run_v2_checkpointed.py
Analysis V2 with robust checkpointing system

Features:
- Batch-by-batch MCMC execution with checkpoints
- Automatic resume from last checkpoint
- Progress saved every N steps
- Can handle crashes and resume

Author: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json
import pickle
import sys
from datetime import datetime
from typing import Dict, Optional, Tuple
import emcee
from tqdm import tqdm

# Local imports
sys.path.insert(0, str(Path(__file__).parent))
from data_loader import (load_pantheon_with_covariance, load_des_sn5yr,
                         combine_datasets)
from janus_model import JanusCosmology, LambdaCDM


class CheckpointedMCMC:
    """
    MCMC runner with checkpoint support for crash recovery
    """

    def __init__(self, z: np.ndarray, mu_obs: np.ndarray,
                 covariance: np.ndarray, H0: float = 70.0,
                 checkpoint_dir: str = 'checkpoints',
                 precision: np.ndarray = None):
        """
        Initialize with data and checkpoint directory

        Parameters
        ----------
        z : array
            Redshifts
        mu_obs : array
            Observed distance moduli
        covariance : array
            Covariance matrix
        H0 : float
            Hubble constant
        checkpoint_dir : str
            Directory for checkpoints
        precision : array, optional
            Precision matrix (inverse covariance). If provided, used directly
            for chi2 calculation (more efficient)
        """
        self.z = z
        self.mu_obs = mu_obs
        self.cov = covariance
        self.n_data = len(z)
        self.H0 = H0

        # Use precision matrix directly if available (more efficient)
        if precision is not None:
            self.precision = precision
            self.use_precision = True
            print(f"  Using precision matrix directly ({precision.shape})")
        else:
            self.precision = None
            self.use_precision = False
            # Pre-compute Cholesky of covariance
            try:
                self.L_cov = np.linalg.cholesky(covariance)
            except np.linalg.LinAlgError:
                print("WARNING: Covariance regularization needed")
                reg = 1e-8 * np.eye(len(covariance))
                self.L_cov = np.linalg.cholesky(covariance + reg)

        # Models
        self.janus = JanusCosmology(H0=H0)
        self.lcdm = LambdaCDM(H0=H0)

        # Checkpoint directory
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def _compute_chi2(self, residuals: np.ndarray) -> float:
        """Compute chi2 using precision or covariance"""
        if self.use_precision:
            # chi2 = r^T @ P @ r (direct with precision matrix)
            return residuals @ self.precision @ residuals
        else:
            # chi2 = r^T @ C^{-1} @ r (via Cholesky)
            y = np.linalg.solve(self.L_cov, residuals)
            return np.dot(y, y)

    def log_likelihood_janus(self, theta: np.ndarray) -> float:
        """Log-likelihood for JANUS"""
        q0, offset = theta
        q0_min = -1.0 / (2.0 * self.z.max() + 0.1)

        if not (q0_min < q0 < 0.0) or not (-10 < offset < 10):
            return -np.inf

        mu_theory = self.janus.distance_modulus(self.z, q0=q0)
        if not np.all(np.isfinite(mu_theory)):
            return -np.inf

        residuals = self.mu_obs - mu_theory - offset
        chi2 = self._compute_chi2(residuals)
        return -0.5 * chi2

    def log_likelihood_lcdm(self, theta: np.ndarray) -> float:
        """Log-likelihood for Lambda-CDM"""
        Omega_m, offset = theta

        if not (0.01 < Omega_m < 0.99) or not (-10 < offset < 10):
            return -np.inf

        # Create LCDM instance with this Omega_m
        lcdm = LambdaCDM(H0=self.H0, Omega_m=Omega_m, Omega_L=1.0-Omega_m)
        mu_theory = lcdm.distance_modulus(self.z)
        residuals = self.mu_obs - mu_theory - offset
        chi2 = self._compute_chi2(residuals)
        return -0.5 * chi2

    def get_checkpoint_path(self, model: str, dataset: str) -> Path:
        """Get checkpoint file path"""
        safe_name = dataset.lower().replace('+', '_').replace(' ', '_').replace('-', '_')
        return self.checkpoint_dir / f"checkpoint_{model}_{safe_name}.pkl"

    def save_checkpoint(self, checkpoint_path: Path, state: Dict):
        """Save checkpoint to disk"""
        # Save to temp file first, then rename (atomic)
        temp_path = checkpoint_path.with_suffix('.tmp')
        with open(temp_path, 'wb') as f:
            pickle.dump(state, f)
        temp_path.rename(checkpoint_path)
        print(f"  Checkpoint saved: {checkpoint_path.name}")

    def load_checkpoint(self, checkpoint_path: Path) -> Optional[Dict]:
        """Load checkpoint if exists"""
        if checkpoint_path.exists():
            with open(checkpoint_path, 'rb') as f:
                state = pickle.load(f)
            print(f"  Checkpoint loaded: {state['completed_steps']} steps completed")
            return state
        return None

    def run_checkpointed(self, model: str, dataset_name: str,
                         n_walkers: int = 32,
                         n_steps_total: int = 2000,
                         batch_size: int = 100,
                         n_burn: int = 500) -> Dict:
        """
        Run MCMC with checkpointing

        Parameters
        ----------
        model : str
            'janus' or 'lcdm'
        dataset_name : str
            Name for checkpoint file
        n_walkers : int
            Number of walkers
        n_steps_total : int
            Total steps to run
        batch_size : int
            Save checkpoint every batch_size steps
        n_burn : int
            Burn-in steps (for final analysis)

        Returns
        -------
        results : dict
            MCMC results
        """
        print(f"\n{'='*60}")
        print(f"MCMC: {model.upper()} on {dataset_name}")
        print(f"{'='*60}")
        print(f"Config: {n_walkers} walkers, {n_steps_total} steps, batch={batch_size}")

        # Setup model
        if model == 'janus':
            log_prob = self.log_likelihood_janus
            param_names = ['q0', 'offset']
            p0_center = np.array([-0.05, 0.0])
            bounds = [(-0.5, 0.0), (-5, 5)]
        elif model == 'lcdm':
            log_prob = self.log_likelihood_lcdm
            param_names = ['Omega_m', 'offset']
            p0_center = np.array([0.3, 0.0])
            bounds = [(0.01, 0.99), (-5, 5)]
        else:
            raise ValueError(f"Unknown model: {model}")

        n_dim = len(param_names)
        checkpoint_path = self.get_checkpoint_path(model, dataset_name)

        # Try to load existing checkpoint
        checkpoint = self.load_checkpoint(checkpoint_path)

        if checkpoint is not None:
            # Resume from checkpoint
            completed_steps = checkpoint['completed_steps']
            all_chains = checkpoint['chains']  # Shape: (completed_steps, n_walkers, n_dim)
            all_log_prob = checkpoint['log_prob']  # Shape: (completed_steps, n_walkers)
            p0 = checkpoint['last_position']  # Last walker positions

            remaining_steps = n_steps_total - completed_steps
            if remaining_steps <= 0:
                print(f"  Already completed {completed_steps} steps, skipping MCMC")
            else:
                print(f"  Resuming from step {completed_steps}, {remaining_steps} remaining")
        else:
            # Fresh start
            completed_steps = 0
            all_chains = []
            all_log_prob = []
            remaining_steps = n_steps_total

            # Initialize walkers
            p0 = p0_center + 0.01 * np.random.randn(n_walkers, n_dim)
            for i, (low, high) in enumerate(bounds):
                p0[:, i] = np.clip(p0[:, i], low + 0.001, high - 0.001)

            print(f"  Starting fresh MCMC run")

        # Run remaining batches
        if remaining_steps > 0:
            sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_prob)

            n_batches = (remaining_steps + batch_size - 1) // batch_size

            for batch_idx in range(n_batches):
                steps_this_batch = min(batch_size, remaining_steps - batch_idx * batch_size)

                print(f"\n  Batch {batch_idx + 1}/{n_batches}: {steps_this_batch} steps")

                # Run this batch
                sampler.reset()
                state = sampler.run_mcmc(p0, steps_this_batch, progress=True)

                # Get results from this batch
                batch_chain = sampler.get_chain()  # (steps, n_walkers, n_dim)
                batch_log_prob = sampler.get_log_prob()  # (steps, n_walkers)

                # Append to accumulated results
                if isinstance(all_chains, list):
                    all_chains = batch_chain
                    all_log_prob = batch_log_prob
                else:
                    all_chains = np.concatenate([all_chains, batch_chain], axis=0)
                    all_log_prob = np.concatenate([all_log_prob, batch_log_prob], axis=0)

                # Update position for next batch
                p0 = state.coords
                completed_steps += steps_this_batch

                # Save checkpoint
                checkpoint_state = {
                    'model': model,
                    'dataset': dataset_name,
                    'n_walkers': n_walkers,
                    'n_dim': n_dim,
                    'param_names': param_names,
                    'completed_steps': completed_steps,
                    'chains': all_chains,
                    'log_prob': all_log_prob,
                    'last_position': p0,
                    'timestamp': datetime.now().isoformat()
                }
                self.save_checkpoint(checkpoint_path, checkpoint_state)

                # Print progress
                acc = np.mean(sampler.acceptance_fraction)
                print(f"    Acceptance: {acc:.2%}")

        # Analyze final results
        print(f"\n  Analyzing {completed_steps} total steps...")

        # Flatten chains after burn-in
        if isinstance(all_chains, list):
            all_chains = np.array(all_chains)
            all_log_prob = np.array(all_log_prob)

        actual_burn = min(n_burn, len(all_chains) - 100)
        samples = all_chains[actual_burn:].reshape(-1, n_dim)
        log_probs_flat = all_log_prob[actual_burn:].flatten()

        # Best-fit
        best_idx = np.argmax(log_probs_flat)
        best_fit = samples[best_idx]

        # Percentiles
        percentiles = np.percentile(samples, [16, 50, 84], axis=0)
        median = percentiles[1]
        sigma_low = median - percentiles[0]
        sigma_high = percentiles[2] - median

        # Chi2
        chi2 = -2 * log_prob(best_fit)
        chi2_dof = chi2 / (self.n_data - n_dim)

        results = {
            'model': model,
            'dataset': dataset_name,
            'param_names': param_names,
            'samples': samples,
            'best_fit': best_fit,
            'median': median,
            'sigma_low': sigma_low,
            'sigma_high': sigma_high,
            'chi2': chi2,
            'chi2_dof': chi2_dof,
            'n_steps': completed_steps,
            'n_burn': actual_burn,
            'n_walkers': n_walkers,
            'n_effective': len(samples)
        }

        # Print results
        print(f"\n  Results:")
        for i, name in enumerate(param_names):
            print(f"    {name} = {median[i]:.4f} +{sigma_high[i]:.4f} -{sigma_low[i]:.4f}")
        print(f"    chi2/dof = {chi2_dof:.4f}")

        return results


def test_data_loading():
    """Test that all datasets load correctly"""
    print("="*60)
    print("TEST: Data Loading")
    print("="*60)

    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / 'data'

    # Test Pantheon+
    print("\n1. Loading Pantheon+...")
    try:
        pantheon = load_pantheon_with_covariance(
            data_path=str(data_dir / 'pantheon' / 'Pantheon+SH0ES.dat'),
            cov_path=str(data_dir / 'pantheon' / 'Pantheon+SH0ES_STAT+SYS.cov')
        )
        print(f"   OK: {pantheon['n_sne']} SNe, z=[{pantheon['z'].min():.3f}, {pantheon['z'].max():.3f}]")
    except Exception as e:
        print(f"   FAILED: {e}")
        return None, None

    # Test DES-SN5YR
    print("\n2. Loading DES-SN5YR...")
    try:
        des = load_des_sn5yr(
            data_path=str(data_dir / 'des_sn5yr' / 'DES-Dovekie_HD.csv'),
            cov_path=str(data_dir / 'des_sn5yr' / 'STAT+SYS.npz')
        )
        print(f"   OK: {des['n_sne']} SNe, z=[{des['z'].min():.3f}, {des['z'].max():.3f}]")
    except Exception as e:
        print(f"   FAILED: {e}")
        return pantheon, None

    # Test combination
    print("\n3. Combining datasets...")
    try:
        combined = combine_datasets(pantheon, des, method='priority_pantheon')
        print(f"   OK: {combined['n_sne']} SNe combined")
    except Exception as e:
        print(f"   FAILED: {e}")

    return pantheon, des


def run_quick_test(pantheon_data: Dict, des_data: Dict,
                   n_steps: int = 200, batch_size: int = 50):
    """Run quick MCMC test on DES data"""
    print("\n" + "="*60)
    print("TEST: Quick MCMC on DES-SN5YR")
    print("="*60)

    base_dir = Path(__file__).parent.parent
    checkpoint_dir = base_dir / 'results' / 'v2_test' / 'checkpoints'

    # Use DES data for test - use precision matrix if available
    precision = des_data.get('precision', None)
    runner = CheckpointedMCMC(
        z=des_data['z'],
        mu_obs=des_data['mu'],
        covariance=des_data['covariance'],
        precision=precision,
        checkpoint_dir=str(checkpoint_dir)
    )

    # Quick JANUS test
    res_janus = runner.run_checkpointed(
        model='janus',
        dataset_name='DES_test',
        n_walkers=16,
        n_steps_total=n_steps,
        batch_size=batch_size,
        n_burn=50
    )

    # Quick LCDM test
    res_lcdm = runner.run_checkpointed(
        model='lcdm',
        dataset_name='DES_test',
        n_walkers=16,
        n_steps_total=n_steps,
        batch_size=batch_size,
        n_burn=50
    )

    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    print(f"JANUS: q0 = {res_janus['median'][0]:.4f} +/- {(res_janus['sigma_low'][0]+res_janus['sigma_high'][0])/2:.4f}")
    print(f"       chi2/dof = {res_janus['chi2_dof']:.3f}")
    print(f"LCDM:  Om = {res_lcdm['median'][0]:.4f} +/- {(res_lcdm['sigma_low'][0]+res_lcdm['sigma_high'][0])/2:.4f}")
    print(f"       chi2/dof = {res_lcdm['chi2_dof']:.3f}")

    delta_chi2 = res_janus['chi2'] - res_lcdm['chi2']
    print(f"\nDelta chi2 (JANUS - LCDM) = {delta_chi2:.1f}")

    return res_janus, res_lcdm


def run_full_v2_analysis(n_steps: int = 2000, batch_size: int = 200):
    """Run complete V2 analysis with checkpoints"""
    print("="*70)
    print("JANUS V2 FULL ANALYSIS")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*70)

    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / 'data'
    results_dir = base_dir / 'results' / 'v2_analysis'
    checkpoint_dir = results_dir / 'checkpoints'
    figures_dir = results_dir / 'figures'

    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Load all data
    print("\nLoading datasets...")
    pantheon = load_pantheon_with_covariance(
        data_path=str(data_dir / 'pantheon' / 'Pantheon+SH0ES.dat'),
        cov_path=str(data_dir / 'pantheon' / 'Pantheon+SH0ES_STAT+SYS.cov')
    )

    des = load_des_sn5yr(
        data_path=str(data_dir / 'des_sn5yr' / 'DES-Dovekie_HD.csv'),
        cov_path=str(data_dir / 'des_sn5yr' / 'STAT+SYS.npz')
    )

    combined = combine_datasets(pantheon, des, method='priority_pantheon')

    all_results = {}

    # Analysis on each dataset
    datasets = [
        ('Pantheon+', pantheon),
        ('DES-SN5YR', des),
        ('Combined', combined)
    ]

    for name, data in datasets:
        runner = CheckpointedMCMC(
            z=data['z'],
            mu_obs=data['mu'],
            covariance=data['covariance'],
            checkpoint_dir=str(checkpoint_dir)
        )

        res_janus = runner.run_checkpointed(
            model='janus',
            dataset_name=name,
            n_steps_total=n_steps,
            batch_size=batch_size,
            n_burn=500
        )

        res_lcdm = runner.run_checkpointed(
            model='lcdm',
            dataset_name=name,
            n_steps_total=n_steps,
            batch_size=batch_size,
            n_burn=500
        )

        # Compute comparison metrics
        delta_chi2 = res_janus['chi2'] - res_lcdm['chi2']
        aic_janus = res_janus['chi2'] + 4  # 2 params
        aic_lcdm = res_lcdm['chi2'] + 4
        delta_aic = aic_janus - aic_lcdm

        bic_janus = res_janus['chi2'] + 2 * np.log(data['n_sne'])
        bic_lcdm = res_lcdm['chi2'] + 2 * np.log(data['n_sne'])
        delta_bic = bic_janus - bic_lcdm

        all_results[name] = {
            'janus': {
                'q0': float(res_janus['median'][0]),
                'q0_err': float((res_janus['sigma_low'][0] + res_janus['sigma_high'][0]) / 2),
                'offset': float(res_janus['median'][1]),
                'chi2': float(res_janus['chi2']),
                'chi2_dof': float(res_janus['chi2_dof'])
            },
            'lcdm': {
                'Omega_m': float(res_lcdm['median'][0]),
                'Omega_m_err': float((res_lcdm['sigma_low'][0] + res_lcdm['sigma_high'][0]) / 2),
                'offset': float(res_lcdm['median'][1]),
                'chi2': float(res_lcdm['chi2']),
                'chi2_dof': float(res_lcdm['chi2_dof'])
            },
            'comparison': {
                'delta_chi2': float(delta_chi2),
                'delta_aic': float(delta_aic),
                'delta_bic': float(delta_bic),
                'n_sne': data['n_sne']
            }
        }

    # Save results
    results_file = results_dir / 'v2_results.json'
    with open(results_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to: {results_file}")

    # Print summary
    print("\n" + "="*70)
    print("V2 ANALYSIS SUMMARY")
    print("="*70)
    print(f"\n{'Dataset':<15} {'N':<6} {'q0 (JANUS)':<18} {'Om (LCDM)':<18} {'dAIC':<8}")
    print("-"*70)

    for name, res in all_results.items():
        n = res['comparison']['n_sne']
        q0 = f"{res['janus']['q0']:.4f}+/-{res['janus']['q0_err']:.4f}"
        Om = f"{res['lcdm']['Omega_m']:.4f}+/-{res['lcdm']['Omega_m_err']:.4f}"
        daic = f"{res['comparison']['delta_aic']:.1f}"
        print(f"{name:<15} {n:<6} {q0:<18} {Om:<18} {daic:<8}")

    return all_results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='JANUS V2 Analysis with Checkpoints')
    parser.add_argument('--test', action='store_true', help='Run quick test only')
    parser.add_argument('--full', action='store_true', help='Run full analysis')
    parser.add_argument('--steps', type=int, default=200, help='Total MCMC steps (test mode)')
    parser.add_argument('--batch', type=int, default=50, help='Checkpoint batch size')
    args = parser.parse_args()

    # Always test data loading first
    pantheon, des = test_data_loading()

    if pantheon is None:
        print("\nERROR: Could not load Pantheon+ data")
        sys.exit(1)

    if des is None:
        print("\nERROR: Could not load DES-SN5YR data")
        sys.exit(1)

    if args.full:
        # Full analysis
        run_full_v2_analysis(n_steps=2000, batch_size=200)
    else:
        # Quick test (default)
        run_quick_test(pantheon, des, n_steps=args.steps, batch_size=args.batch)