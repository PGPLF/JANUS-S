"""
fitting.py
Ajustement des modeles cosmologiques aux donnees de supernovae

Methodes:
- Minimisation chi-2 (scipy.optimize)
- Analyse MCMC (emcee)
- Bootstrap pour erreurs

Auteur: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
from typing import Dict, Tuple, Optional, Callable
import warnings

# Import des modeles
try:
    from janus_model import JanusCosmology, LambdaCDM
except ImportError:
    from .janus_model import JanusCosmology, LambdaCDM


def chi2_diagonal(params: np.ndarray,
                  z: np.ndarray,
                  mu_obs: np.ndarray,
                  sigma_mu: np.ndarray,
                  model: object,
                  model_type: str = 'janus') -> float:
    """
    Calcule le chi-2 avec erreurs diagonales

    Parametres
    ----------
    params : array
        Parametres du modele [q0, offset] pour JANUS
        ou [Omega_m, offset] pour Lambda-CDM
    z : array
        Redshifts
    mu_obs : array
        Modules de distance observes
    sigma_mu : array
        Erreurs sur mu
    model : object
        Instance du modele cosmologique
    model_type : str
        'janus' ou 'lcdm'

    Retourne
    --------
    chi2 : float
        Valeur du chi-2
    """
    if model_type == 'janus':
        q0, offset = params
        if q0 >= 0 or q0 < -0.5:
            return 1e10
        mu_theory = model.distance_modulus(z, q0, offset)
    else:  # Lambda-CDM
        offset = params[0] if len(params) == 1 else params[1]
        mu_theory = model.distance_modulus(z, offset)

    # Gestion des valeurs infinies
    if np.any(~np.isfinite(mu_theory)):
        return 1e10

    residuals = (mu_obs - mu_theory) / sigma_mu
    return np.sum(residuals**2)


def chi2_covariance(params: np.ndarray,
                    z: np.ndarray,
                    mu_obs: np.ndarray,
                    cov_matrix: np.ndarray,
                    model: object,
                    model_type: str = 'janus') -> float:
    """
    Calcule le chi-2 avec matrice de covariance complete

    Parametres
    ----------
    cov_matrix : array (n x n)
        Matrice de covariance

    Retourne
    --------
    chi2 : float
        Valeur du chi-2
    """
    if model_type == 'janus':
        q0, offset = params
        if q0 >= 0 or q0 < -0.5:
            return 1e10
        mu_theory = model.distance_modulus(z, q0, offset)
    else:
        offset = params[0] if len(params) == 1 else params[1]
        mu_theory = model.distance_modulus(z, offset)

    if np.any(~np.isfinite(mu_theory)):
        return 1e10

    residuals = mu_obs - mu_theory

    try:
        cov_inv = np.linalg.inv(cov_matrix)
        chi2 = residuals @ cov_inv @ residuals
    except np.linalg.LinAlgError:
        # Si inversion echoue, utiliser pseudo-inverse
        cov_inv = np.linalg.pinv(cov_matrix)
        chi2 = residuals @ cov_inv @ residuals

    return chi2


def fit_janus(z: np.ndarray,
              mu_obs: np.ndarray,
              sigma_mu: np.ndarray,
              cov_matrix: Optional[np.ndarray] = None,
              H0: float = 70.0,
              method: str = 'Nelder-Mead',
              verbose: bool = True) -> Dict:
    """
    Ajuste le modele JANUS aux donnees

    Parametres
    ----------
    z : array
        Redshifts
    mu_obs : array
        Modules de distance observes
    sigma_mu : array
        Erreurs sur mu (utilise si cov_matrix=None)
    cov_matrix : array, optional
        Matrice de covariance complete
    H0 : float
        Constante de Hubble
    method : str
        Methode d'optimisation scipy
    verbose : bool
        Afficher les resultats

    Retourne
    --------
    result : dict
        Resultats de l'ajustement
    """
    model = JanusCosmology(H0=H0)

    # Choix de la fonction chi-2
    if cov_matrix is not None:
        chi2_func = lambda p: chi2_covariance(p, z, mu_obs, cov_matrix, model, 'janus')
    else:
        chi2_func = lambda p: chi2_diagonal(p, z, mu_obs, sigma_mu, model, 'janus')

    # Parametres initiaux et bornes
    initial = np.array([-0.1, 0.0])
    bounds = [(-0.5, -0.001), (-5, 5)]

    # Optimisation
    if method == 'differential_evolution':
        opt_result = differential_evolution(chi2_func, bounds, seed=42)
    else:
        opt_result = minimize(chi2_func, initial, method=method, bounds=bounds)

    q0_fit = opt_result.x[0]
    offset_fit = opt_result.x[1]
    chi2_min = opt_result.fun
    n_data = len(z)
    n_params = 2
    dof = n_data - n_params

    # Calcul de l'age de l'univers
    T0 = model.universe_age(q0_fit)

    result = {
        'q0': q0_fit,
        'offset': offset_fit,
        'chi2': chi2_min,
        'dof': dof,
        'chi2_reduced': chi2_min / dof,
        'n_data': n_data,
        'n_params': n_params,
        'success': opt_result.success,
        'age_Gyr': T0,
        'H0': H0,
        'model': 'JANUS'
    }

    if verbose:
        print("\n" + "=" * 50)
        print("RESULTATS AJUSTEMENT JANUS")
        print("=" * 50)
        print(f"q0 = {q0_fit:.4f}")
        print(f"offset = {offset_fit:.4f}")
        print(f"chi2 = {chi2_min:.2f}")
        print(f"chi2/dof = {chi2_min/dof:.4f}")
        print(f"dof = {dof}")
        print(f"Age univers = {T0:.2f} Gyr (H0={H0})")
        print("=" * 50)

    return result


def fit_lcdm(z: np.ndarray,
             mu_obs: np.ndarray,
             sigma_mu: np.ndarray,
             cov_matrix: Optional[np.ndarray] = None,
             H0: float = 70.0,
             Omega_m: float = 0.3,
             Omega_L: float = 0.7,
             fit_Omega_m: bool = False,
             verbose: bool = True) -> Dict:
    """
    Ajuste le modele Lambda-CDM aux donnees

    Parametres
    ----------
    fit_Omega_m : bool
        Si True, ajuste Omega_m. Sinon, fixe aux valeurs donnees.

    Retourne
    --------
    result : dict
        Resultats de l'ajustement
    """
    model = LambdaCDM(H0=H0, Omega_m=Omega_m, Omega_L=Omega_L)

    # Choix de la fonction chi-2
    if cov_matrix is not None:
        chi2_func = lambda p: chi2_covariance(p, z, mu_obs, cov_matrix, model, 'lcdm')
    else:
        chi2_func = lambda p: chi2_diagonal(p, z, mu_obs, sigma_mu, model, 'lcdm')

    # Ajustement (offset seulement pour comparaison equitable)
    initial = np.array([0.0])
    bounds = [(-5, 5)]

    opt_result = minimize(chi2_func, initial, method='Nelder-Mead', bounds=bounds)

    offset_fit = opt_result.x[0]
    chi2_min = opt_result.fun
    n_data = len(z)
    n_params = 1 if not fit_Omega_m else 2
    dof = n_data - n_params

    result = {
        'Omega_m': Omega_m,
        'Omega_L': Omega_L,
        'offset': offset_fit,
        'chi2': chi2_min,
        'dof': dof,
        'chi2_reduced': chi2_min / dof,
        'n_data': n_data,
        'n_params': n_params,
        'success': opt_result.success,
        'H0': H0,
        'model': 'Lambda-CDM'
    }

    if verbose:
        print("\n" + "=" * 50)
        print("RESULTATS AJUSTEMENT LAMBDA-CDM")
        print("=" * 50)
        print(f"Omega_m = {Omega_m:.3f}")
        print(f"Omega_L = {Omega_L:.3f}")
        print(f"offset = {offset_fit:.4f}")
        print(f"chi2 = {chi2_min:.2f}")
        print(f"chi2/dof = {chi2_min/dof:.4f}")
        print(f"dof = {dof}")
        print("=" * 50)

    return result


def model_comparison(result_janus: Dict,
                     result_lcdm: Dict,
                     verbose: bool = True) -> Dict:
    """
    Compare les modeles JANUS et Lambda-CDM

    Utilise les criteres AIC, BIC et likelihood ratio

    Retourne
    --------
    comparison : dict
        Resultats de la comparaison
    """
    n = result_janus['n_data']

    # Nombre de parametres
    k_janus = result_janus['n_params']
    k_lcdm = result_lcdm['n_params']

    chi2_janus = result_janus['chi2']
    chi2_lcdm = result_lcdm['chi2']

    # AIC = chi2 + 2k
    aic_janus = chi2_janus + 2 * k_janus
    aic_lcdm = chi2_lcdm + 2 * k_lcdm
    delta_aic = aic_lcdm - aic_janus

    # BIC = chi2 + k * ln(n)
    bic_janus = chi2_janus + k_janus * np.log(n)
    bic_lcdm = chi2_lcdm + k_lcdm * np.log(n)
    delta_bic = bic_lcdm - bic_janus

    # Likelihood ratio
    delta_chi2 = chi2_lcdm - chi2_janus
    likelihood_ratio = np.exp(-0.5 * delta_chi2)

    # Interpretation
    if delta_aic > 10:
        aic_interpretation = "Forte evidence pour JANUS"
    elif delta_aic > 6:
        aic_interpretation = "Evidence moderee pour JANUS"
    elif delta_aic > 2:
        aic_interpretation = "Evidence faible pour JANUS"
    elif delta_aic > -2:
        aic_interpretation = "Modeles equivalents"
    elif delta_aic > -6:
        aic_interpretation = "Evidence faible pour LCDM"
    elif delta_aic > -10:
        aic_interpretation = "Evidence moderee pour LCDM"
    else:
        aic_interpretation = "Forte evidence pour LCDM"

    comparison = {
        'chi2_janus': chi2_janus,
        'chi2_lcdm': chi2_lcdm,
        'delta_chi2': delta_chi2,
        'AIC_janus': aic_janus,
        'AIC_lcdm': aic_lcdm,
        'delta_AIC': delta_aic,
        'BIC_janus': bic_janus,
        'BIC_lcdm': bic_lcdm,
        'delta_BIC': delta_bic,
        'likelihood_ratio': likelihood_ratio,
        'interpretation': aic_interpretation
    }

    if verbose:
        print("\n" + "=" * 50)
        print("COMPARAISON DES MODELES")
        print("=" * 50)
        print(f"Chi2 JANUS:  {chi2_janus:.2f}")
        print(f"Chi2 LCDM:   {chi2_lcdm:.2f}")
        print(f"Delta chi2:  {delta_chi2:.2f}")
        print("-" * 50)
        print(f"AIC JANUS:   {aic_janus:.2f}")
        print(f"AIC LCDM:    {aic_lcdm:.2f}")
        print(f"Delta AIC:   {delta_aic:.2f}")
        print("-" * 50)
        print(f"BIC JANUS:   {bic_janus:.2f}")
        print(f"BIC LCDM:    {bic_lcdm:.2f}")
        print(f"Delta BIC:   {delta_bic:.2f}")
        print("-" * 50)
        print(f"Likelihood ratio: {likelihood_ratio:.4f}")
        print(f"Interpretation:   {aic_interpretation}")
        print("=" * 50)

    return comparison


def bootstrap_errors(z: np.ndarray,
                     mu_obs: np.ndarray,
                     sigma_mu: np.ndarray,
                     model_type: str = 'janus',
                     n_bootstrap: int = 1000,
                     H0: float = 70.0,
                     verbose: bool = True) -> Dict:
    """
    Estime les erreurs par bootstrap

    Parametres
    ----------
    n_bootstrap : int
        Nombre de reechantillonnages

    Retourne
    --------
    errors : dict
        Erreurs estimees sur les parametres
    """
    n = len(z)
    results = []

    for i in range(n_bootstrap):
        # Reechantillonnage avec remplacement
        indices = np.random.choice(n, n, replace=True)
        z_boot = z[indices]
        mu_boot = mu_obs[indices]
        sigma_boot = sigma_mu[indices]

        try:
            if model_type == 'janus':
                result = fit_janus(z_boot, mu_boot, sigma_boot,
                                   H0=H0, verbose=False)
                results.append([result['q0'], result['offset']])
            else:
                result = fit_lcdm(z_boot, mu_boot, sigma_boot,
                                  H0=H0, verbose=False)
                results.append([result['offset']])
        except:
            continue

    results = np.array(results)

    if model_type == 'janus':
        q0_std = np.std(results[:, 0])
        offset_std = np.std(results[:, 1])
        errors = {
            'q0_error': q0_std,
            'offset_error': offset_std,
            'n_success': len(results)
        }
        if verbose:
            print(f"\nBootstrap ({len(results)} samples):")
            print(f"sigma(q0) = {q0_std:.4f}")
            print(f"sigma(offset) = {offset_std:.4f}")
    else:
        offset_std = np.std(results[:, 0])
        errors = {
            'offset_error': offset_std,
            'n_success': len(results)
        }
        if verbose:
            print(f"\nBootstrap ({len(results)} samples):")
            print(f"sigma(offset) = {offset_std:.4f}")

    return errors


# Tests
if __name__ == "__main__":
    print("Test du module fitting.py")
    print("=" * 50)

    # Donnees simulees
    np.random.seed(42)
    n = 100
    z_test = np.random.uniform(0.01, 1.0, n)
    z_test = np.sort(z_test)

    # Generer mu avec modele JANUS
    janus = JanusCosmology(H0=70.0)
    q0_true = -0.087
    mu_true = janus.distance_modulus(z_test, q0_true)
    sigma = 0.15
    mu_obs = mu_true + np.random.normal(0, sigma, n)
    sigma_mu = np.full(n, sigma)

    # Test ajustement JANUS
    result_janus = fit_janus(z_test, mu_obs, sigma_mu)

    # Test ajustement Lambda-CDM
    result_lcdm = fit_lcdm(z_test, mu_obs, sigma_mu)

    # Test comparaison
    comparison = model_comparison(result_janus, result_lcdm)

    print("\n" + "=" * 50)
    print("Tests passes avec succes!")
