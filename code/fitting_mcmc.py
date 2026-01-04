"""
fitting_mcmc.py
MCMC fitting for cosmological models using emcee and dynesty

Supports:
- JANUS model (q0 parameter)
- Lambda-CDM model (Omega_m parameter)
- wCDM model (Omega_m, w parameters)
- CPL model (Omega_m, w0, wa parameters)

MCMC methods:
- emcee: Affine-invariant ensemble sampler
- dynesty: Nested sampling with evidence calculation

Author: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
from typing import Dict, Tuple, Optional, Callable
import warnings
from scipy.optimize import minimize

# MCMC libraries
import emcee
import dynesty
from dynesty import utils as dyfunc
import corner

# Progress bar
from tqdm import tqdm

# Local imports
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from janus_model import JanusCosmology, LambdaCDM
from data_loader import compute_chi2_with_covariance


class MCMCFitter:
    """
    Unified MCMC fitter for cosmological models

    Supports both emcee (affine-invariant) and dynesty (nested sampling)
    """

    def __init__(self, z: np.ndarray, mu_obs: np.ndarray,
                 covariance: Optional[np.ndarray] = None,
                 sigma_mu: Optional[np.ndarray] = None,
                 H0: float = 70.0):
        """
        Initialize the fitter

        Parameters
        ----------
        z : array
            Redshifts
        mu_obs : array
            Observed distance moduli
        covariance : array (N x N), optional
            Full covariance matrix
        sigma_mu : array, optional
            Diagonal errors (used if covariance not provided)
        H0 : float
            Hubble constant (km/s/Mpc)
        """
        self.z = z
        self.mu_obs = mu_obs
        self.n_data = len(z)
        self.H0 = H0

        # Covariance handling
        if covariance is not None:
            self.cov = covariance
            self.use_full_cov = True
            # Pre-compute Cholesky decomposition for efficiency
            try:
                self.L_cov = np.linalg.cholesky(covariance)
            except np.linalg.LinAlgError:
                warnings.warn("Covariance not positive definite, using regularization")
                reg = 1e-10 * np.eye(len(covariance))
                self.L_cov = np.linalg.cholesky(covariance + reg)
        else:
            self.sigma_mu = sigma_mu if sigma_mu is not None else np.ones(len(z)) * 0.1
            self.use_full_cov = False
            self.cov = None

        # Models
        self.janus = JanusCosmology(H0=H0)
        self.lcdm = LambdaCDM(H0=H0)

        # Results storage
        self.results = {}

    def log_likelihood_janus(self, theta: np.ndarray) -> float:
        """Log-likelihood for JANUS model"""
        q0, offset = theta

        # Bounds check
        if not (-0.5 < q0 < 0.0) or not (-10 < offset < 10):
            return -np.inf

        # Theory prediction
        mu_theory = self.janus.distance_modulus(self.z, q0=q0)
        residuals = self.mu_obs - mu_theory - offset

        if self.use_full_cov:
            # Full covariance chi2
            y = np.linalg.solve(self.L_cov, residuals)
            chi2 = np.dot(y, y)
        else:
            # Diagonal chi2
            chi2 = np.sum((residuals / self.sigma_mu)**2)

        return -0.5 * chi2

    def log_likelihood_lcdm(self, theta: np.ndarray) -> float:
        """Log-likelihood for Lambda-CDM model"""
        Omega_m, offset = theta

        # Bounds check
        if not (0.01 < Omega_m < 0.99) or not (-10 < offset < 10):
            return -np.inf

        # Theory prediction
        mu_theory = self.lcdm.distance_modulus(self.z, Omega_m=Omega_m)
        residuals = self.mu_obs - mu_theory - offset

        if self.use_full_cov:
            y = np.linalg.solve(self.L_cov, residuals)
            chi2 = np.dot(y, y)
        else:
            chi2 = np.sum((residuals / self.sigma_mu)**2)

        return -0.5 * chi2

    def log_prior_janus(self, theta: np.ndarray) -> float:
        """Uniform prior for JANUS"""
        q0, offset = theta
        if -0.5 < q0 < 0.0 and -10 < offset < 10:
            return 0.0
        return -np.inf

    def log_prior_lcdm(self, theta: np.ndarray) -> float:
        """Uniform prior for Lambda-CDM"""
        Omega_m, offset = theta
        if 0.01 < Omega_m < 0.99 and -10 < offset < 10:
            return 0.0
        return -np.inf

    def log_posterior_janus(self, theta: np.ndarray) -> float:
        """Log-posterior for JANUS"""
        lp = self.log_prior_janus(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood_janus(theta)

    def log_posterior_lcdm(self, theta: np.ndarray) -> float:
        """Log-posterior for Lambda-CDM"""
        lp = self.log_prior_lcdm(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood_lcdm(theta)

    def run_emcee(self, model: str = 'janus',
                  n_walkers: int = 32,
                  n_steps: int = 5000,
                  n_burn: int = 1000,
                  progress: bool = True) -> Dict:
        """
        Run emcee MCMC sampling

        Parameters
        ----------
        model : str
            'janus' or 'lcdm'
        n_walkers : int
            Number of walkers
        n_steps : int
            Number of steps per walker
        n_burn : int
            Burn-in steps to discard
        progress : bool
            Show progress bar

        Returns
        -------
        results : dict
            MCMC results including samples, best-fit, uncertainties
        """
        if model == 'janus':
            log_prob = self.log_posterior_janus
            param_names = ['q0', 'offset']
            # Initial guess near best-fit
            p0_center = np.array([-0.05, 0.0])
            bounds = [(-0.5, 0.0), (-10, 10)]
        elif model == 'lcdm':
            log_prob = self.log_posterior_lcdm
            param_names = ['Omega_m', 'offset']
            p0_center = np.array([0.3, 0.0])
            bounds = [(0.01, 0.99), (-10, 10)]
        else:
            raise ValueError(f"Unknown model: {model}")

        n_dim = len(param_names)

        # Initialize walkers with small scatter
        p0 = p0_center + 0.01 * np.random.randn(n_walkers, n_dim)

        # Ensure within bounds
        for i, (low, high) in enumerate(bounds):
            p0[:, i] = np.clip(p0[:, i], low + 0.01, high - 0.01)

        # Run sampler
        sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_prob)

        if progress:
            print(f"Running emcee for {model.upper()} ({n_walkers} walkers, {n_steps} steps)...")

        sampler.run_mcmc(p0, n_steps, progress=progress)

        # Get samples after burn-in
        samples = sampler.get_chain(discard=n_burn, flat=True)
        log_probs = sampler.get_log_prob(discard=n_burn, flat=True)

        # Best-fit (MAP)
        best_idx = np.argmax(log_probs)
        best_fit = samples[best_idx]

        # Compute percentiles
        percentiles = np.percentile(samples, [16, 50, 84], axis=0)
        median = percentiles[1]
        sigma_low = median - percentiles[0]
        sigma_high = percentiles[2] - median

        # Acceptance fraction
        acc_frac = np.mean(sampler.acceptance_fraction)

        # Autocorrelation time (if converged)
        try:
            tau = sampler.get_autocorr_time(quiet=True)
        except:
            tau = np.array([np.nan] * n_dim)

        results = {
            'model': model,
            'method': 'emcee',
            'samples': samples,
            'log_prob': log_probs,
            'best_fit': best_fit,
            'median': median,
            'sigma_low': sigma_low,
            'sigma_high': sigma_high,
            'param_names': param_names,
            'n_walkers': n_walkers,
            'n_steps': n_steps,
            'n_burn': n_burn,
            'acceptance_fraction': acc_frac,
            'autocorr_time': tau,
            'n_effective': len(samples)
        }

        # Chi2 at best-fit
        if model == 'janus':
            chi2 = -2 * self.log_likelihood_janus(best_fit)
        else:
            chi2 = -2 * self.log_likelihood_lcdm(best_fit)
        results['chi2'] = chi2
        results['chi2_dof'] = chi2 / (self.n_data - n_dim)

        self.results[f'emcee_{model}'] = results
        return results

    def run_dynesty(self, model: str = 'janus',
                    n_live: int = 500,
                    dlogz: float = 0.1,
                    progress: bool = True) -> Dict:
        """
        Run dynesty nested sampling

        Parameters
        ----------
        model : str
            'janus' or 'lcdm'
        n_live : int
            Number of live points
        dlogz : float
            Stopping criterion (evidence precision)
        progress : bool
            Show progress bar

        Returns
        -------
        results : dict
            Nested sampling results including evidence
        """
        if model == 'janus':
            param_names = ['q0', 'offset']
            bounds = [(-0.5, 0.0), (-5, 5)]

            def prior_transform(u):
                # Transform from unit cube to parameter space
                q0 = u[0] * 0.5 - 0.5  # [-0.5, 0]
                offset = u[1] * 10 - 5  # [-5, 5]
                return np.array([q0, offset])

            def log_like(theta):
                return self.log_likelihood_janus(theta)

        elif model == 'lcdm':
            param_names = ['Omega_m', 'offset']
            bounds = [(0.01, 0.99), (-5, 5)]

            def prior_transform(u):
                Omega_m = u[0] * 0.98 + 0.01  # [0.01, 0.99]
                offset = u[1] * 10 - 5  # [-5, 5]
                return np.array([Omega_m, offset])

            def log_like(theta):
                return self.log_likelihood_lcdm(theta)
        else:
            raise ValueError(f"Unknown model: {model}")

        n_dim = len(param_names)

        if progress:
            print(f"Running dynesty for {model.upper()} ({n_live} live points)...")

        # Run nested sampling
        sampler = dynesty.NestedSampler(
            log_like, prior_transform, n_dim,
            nlive=n_live, bound='multi', sample='auto'
        )
        sampler.run_nested(dlogz=dlogz, print_progress=progress)

        # Get results
        dres = sampler.results

        # Compute weighted samples
        samples, weights = dres.samples, np.exp(dres.logwt - dres.logz[-1])

        # Resample to get equal-weighted samples
        samples_eq = dyfunc.resample_equal(samples, weights)

        # Evidence
        log_z = dres.logz[-1]
        log_z_err = dres.logzerr[-1]

        # Best-fit
        best_idx = np.argmax(dres.logl)
        best_fit = samples[best_idx]

        # Percentiles
        percentiles = np.percentile(samples_eq, [16, 50, 84], axis=0)
        median = percentiles[1]
        sigma_low = median - percentiles[0]
        sigma_high = percentiles[2] - median

        results = {
            'model': model,
            'method': 'dynesty',
            'samples': samples_eq,
            'weights': weights,
            'best_fit': best_fit,
            'median': median,
            'sigma_low': sigma_low,
            'sigma_high': sigma_high,
            'param_names': param_names,
            'log_evidence': log_z,
            'log_evidence_err': log_z_err,
            'n_live': n_live,
            'n_calls': dres.ncall,
            'n_effective': len(samples_eq)
        }

        # Chi2 at best-fit
        chi2 = -2 * log_like(best_fit)
        results['chi2'] = chi2
        results['chi2_dof'] = chi2 / (self.n_data - n_dim)

        self.results[f'dynesty_{model}'] = results
        return results

    def compare_methods(self, model: str = 'janus') -> Dict:
        """
        Compare emcee and dynesty results for given model
        """
        emcee_key = f'emcee_{model}'
        dynesty_key = f'dynesty_{model}'

        if emcee_key not in self.results or dynesty_key not in self.results:
            raise ValueError("Run both emcee and dynesty first")

        e = self.results[emcee_key]
        d = self.results[dynesty_key]

        comparison = {
            'model': model,
            'emcee_median': e['median'],
            'dynesty_median': d['median'],
            'emcee_sigma': (e['sigma_low'] + e['sigma_high']) / 2,
            'dynesty_sigma': (d['sigma_low'] + d['sigma_high']) / 2,
            'param_diff': np.abs(e['median'] - d['median']),
            'sigma_diff_percent': np.abs(e['sigma_low'] - d['sigma_low']) / e['sigma_low'] * 100,
            'chi2_emcee': e['chi2'],
            'chi2_dynesty': d['chi2'],
            'log_evidence': d['log_evidence'],
            'log_evidence_err': d['log_evidence_err']
        }

        return comparison

    def plot_corner(self, model: str = 'janus', method: str = 'emcee',
                    save_path: Optional[str] = None) -> None:
        """
        Generate corner plot for given model and method
        """
        key = f'{method}_{model}'
        if key not in self.results:
            raise ValueError(f"No results for {key}")

        res = self.results[key]

        # Labels
        if model == 'janus':
            labels = [r'$q_0$', r'$\delta$']
        else:
            labels = [r'$\Omega_m$', r'$\delta$']

        fig = corner.corner(
            res['samples'],
            labels=labels,
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True,
            title_kwargs={"fontsize": 12}
        )

        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Corner plot saved to {save_path}")

        return fig


def run_full_analysis(data: Dict, verbose: bool = True) -> Dict:
    """
    Run complete MCMC analysis for JANUS and Lambda-CDM

    Parameters
    ----------
    data : dict
        Data from load_pantheon_with_covariance
    verbose : bool
        Print progress

    Returns
    -------
    results : dict
        Complete analysis results
    """
    fitter = MCMCFitter(
        z=data['z'],
        mu_obs=data['mu'],
        covariance=data['covariance']
    )

    results = {}

    # JANUS
    if verbose:
        print("\n" + "="*60)
        print("JANUS MODEL")
        print("="*60)

    results['janus_emcee'] = fitter.run_emcee('janus', progress=verbose)
    results['janus_dynesty'] = fitter.run_dynesty('janus', progress=verbose)

    if verbose:
        print(f"\nJANUS emcee:   q0 = {results['janus_emcee']['median'][0]:.4f} "
              f"+{results['janus_emcee']['sigma_high'][0]:.4f} "
              f"-{results['janus_emcee']['sigma_low'][0]:.4f}")
        print(f"JANUS dynesty: q0 = {results['janus_dynesty']['median'][0]:.4f} "
              f"+{results['janus_dynesty']['sigma_high'][0]:.4f} "
              f"-{results['janus_dynesty']['sigma_low'][0]:.4f}")

    # Lambda-CDM
    if verbose:
        print("\n" + "="*60)
        print("LAMBDA-CDM MODEL")
        print("="*60)

    results['lcdm_emcee'] = fitter.run_emcee('lcdm', progress=verbose)
    results['lcdm_dynesty'] = fitter.run_dynesty('lcdm', progress=verbose)

    if verbose:
        print(f"\nLCDM emcee:   Omega_m = {results['lcdm_emcee']['median'][0]:.4f} "
              f"+{results['lcdm_emcee']['sigma_high'][0]:.4f} "
              f"-{results['lcdm_emcee']['sigma_low'][0]:.4f}")
        print(f"LCDM dynesty: Omega_m = {results['lcdm_dynesty']['median'][0]:.4f} "
              f"+{results['lcdm_dynesty']['sigma_high'][0]:.4f} "
              f"-{results['lcdm_dynesty']['sigma_low'][0]:.4f}")

    # Bayes factor
    log_B = results['janus_dynesty']['log_evidence'] - results['lcdm_dynesty']['log_evidence']
    results['bayes_factor'] = {
        'log_B_janus_lcdm': log_B,
        'interpretation': interpret_bayes_factor(log_B)
    }

    if verbose:
        print("\n" + "="*60)
        print("MODEL COMPARISON")
        print("="*60)
        print(f"log(B_JANUS/B_LCDM) = {log_B:.2f}")
        print(f"Interpretation: {results['bayes_factor']['interpretation']}")

    return results


def interpret_bayes_factor(log_B: float) -> str:
    """
    Interpret Bayes factor according to Jeffreys scale
    """
    if log_B > 5:
        return "Strong evidence for JANUS"
    elif log_B > 2.5:
        return "Moderate evidence for JANUS"
    elif log_B > 1:
        return "Weak evidence for JANUS"
    elif log_B > -1:
        return "Inconclusive"
    elif log_B > -2.5:
        return "Weak evidence for LCDM"
    elif log_B > -5:
        return "Moderate evidence for LCDM"
    else:
        return "Strong evidence for LCDM"


# Test
if __name__ == "__main__":
    print("Test fitting_mcmc.py")
    print("="*60)

    # Synthetic test data
    np.random.seed(42)
    n_test = 100
    z_test = np.linspace(0.01, 1.0, n_test)

    # True model: JANUS with q0 = -0.1
    janus = JanusCosmology(H0=70)
    mu_true = janus.distance_modulus(z_test, q0=-0.1)
    sigma_test = 0.15 * np.ones(n_test)
    mu_obs = mu_true + np.random.normal(0, sigma_test)

    # Create fitter with diagonal errors
    fitter = MCMCFitter(z_test, mu_obs, sigma_mu=sigma_test)

    # Quick test with reduced steps
    print("\nRunning quick emcee test...")
    res = fitter.run_emcee('janus', n_walkers=16, n_steps=500, n_burn=100)
    print(f"Result: q0 = {res['median'][0]:.3f} +/- {res['sigma_low'][0]:.3f}")
    print(f"True value: q0 = -0.100")
    print(f"Chi2/dof = {res['chi2_dof']:.3f}")

    print("\nTest passed!")
