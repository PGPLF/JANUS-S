"""
data_loader.py
Chargement et preparation des donnees de supernovae

Datasets supportes:
- JLA (Joint Light-curve Analysis) - Betoule et al. (2014)
- Pantheon+ - Brout et al. (2022)

Auteur: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple, Optional
import warnings


def load_jla_data(filepath: str = 'data/jla/jla_lcparams.txt',
                  verbose: bool = True) -> Dict:
    """
    Charge le dataset JLA (740 supernovae de type Ia)

    Parametres
    ----------
    filepath : str
        Chemin vers le fichier jla_lcparams.txt
    verbose : bool
        Afficher les informations de chargement

    Retourne
    --------
    data : dict
        Dictionnaire contenant:
        - z : redshifts (z_cmb)
        - mb : magnitudes apparentes (m_B^*)
        - x1 : parametres d'etirement (X_1)
        - color : couleurs (C)
        - dmb, dx1, dcolor : erreurs associees
        - n_sne : nombre de supernovae
        - names : noms des supernovae
    """
    path = Path(filepath)

    if not path.exists():
        raise FileNotFoundError(
            f"Fichier non trouve: {filepath}\n"
            "Telechargez les donnees depuis: "
            "https://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html"
        )

    # Lecture du fichier
    try:
        data = pd.read_csv(filepath, delim_whitespace=True)
    except Exception as e:
        raise IOError(f"Erreur de lecture du fichier: {e}")

    # Verification des colonnes attendues
    required_cols = ['zcmb', 'mb', 'x1', 'color', 'dmb']
    for col in required_cols:
        if col not in data.columns:
            # Essayer des noms alternatifs
            alt_names = {
                'zcmb': ['z_cmb', 'z', 'redshift'],
                'mb': ['m_b', 'mB', 'mag'],
                'x1': ['X1', 'stretch'],
                'color': ['c', 'C'],
                'dmb': ['e_mb', 'err_mb', 'sigma_mb']
            }
            found = False
            for alt in alt_names.get(col, []):
                if alt in data.columns:
                    data = data.rename(columns={alt: col})
                    found = True
                    break
            if not found:
                warnings.warn(f"Colonne '{col}' non trouvee dans le fichier")

    result = {
        'z': data['zcmb'].values if 'zcmb' in data.columns else None,
        'mb': data['mb'].values if 'mb' in data.columns else None,
        'x1': data['x1'].values if 'x1' in data.columns else None,
        'color': data['color'].values if 'color' in data.columns else None,
        'dmb': data['dmb'].values if 'dmb' in data.columns else None,
        'dx1': data.get('dx1', pd.Series(np.zeros(len(data)))).values,
        'dcolor': data.get('dcolor', pd.Series(np.zeros(len(data)))).values,
        'n_sne': len(data),
        'names': data.get('name', pd.Series(range(len(data)))).values,
        'dataset': 'JLA',
        'reference': 'Betoule et al. (2014)'
    }

    if verbose:
        print(f"Dataset JLA charge: {result['n_sne']} supernovae")
        print(f"Plage de redshift: z = [{result['z'].min():.3f}, {result['z'].max():.3f}]")

    return result


def load_pantheon_data(filepath: str = 'data/pantheon/Pantheon+SH0ES.dat',
                       verbose: bool = True) -> Dict:
    """
    Charge le dataset Pantheon+ (1550 supernovae de type Ia)

    Parametres
    ----------
    filepath : str
        Chemin vers le fichier Pantheon+SH0ES.dat
    verbose : bool
        Afficher les informations de chargement

    Retourne
    --------
    data : dict
        Dictionnaire avec les memes cles que load_jla_data
    """
    path = Path(filepath)

    if not path.exists():
        raise FileNotFoundError(
            f"Fichier non trouve: {filepath}\n"
            "Telechargez les donnees depuis: "
            "https://github.com/PantheonPlusSH0ES"
        )

    # Lecture du fichier Pantheon+
    try:
        data = pd.read_csv(filepath, delim_whitespace=True, comment='#')
    except Exception as e:
        raise IOError(f"Erreur de lecture du fichier: {e}")

    # Mapping des colonnes Pantheon+ vers format standard
    col_mapping = {
        'zHD': 'zcmb',  # Redshift Hubble diagram
        'mB': 'mb',
        'mBERR': 'dmb',
        'x1': 'x1',
        'x1ERR': 'dx1',
        'c': 'color',
        'cERR': 'dcolor'
    }

    for old, new in col_mapping.items():
        if old in data.columns:
            data = data.rename(columns={old: new})

    result = {
        'z': data['zcmb'].values if 'zcmb' in data.columns else data['zHD'].values,
        'mb': data['mb'].values if 'mb' in data.columns else data['mB'].values,
        'x1': data['x1'].values if 'x1' in data.columns else np.zeros(len(data)),
        'color': data['color'].values if 'color' in data.columns else data['c'].values,
        'dmb': data['dmb'].values if 'dmb' in data.columns else data['mBERR'].values,
        'dx1': data.get('dx1', pd.Series(np.zeros(len(data)))).values,
        'dcolor': data.get('dcolor', pd.Series(np.zeros(len(data)))).values,
        'n_sne': len(data),
        'names': data.get('CID', data.get('name', pd.Series(range(len(data))))).values,
        'dataset': 'Pantheon+',
        'reference': 'Brout et al. (2022)'
    }

    if verbose:
        print(f"Dataset Pantheon+ charge: {result['n_sne']} supernovae")
        print(f"Plage de redshift: z = [{result['z'].min():.3f}, {result['z'].max():.3f}]")

    return result


def compute_distance_modulus(data: Dict,
                             alpha: float = 0.141,
                             beta: float = 3.101,
                             M_B: float = -19.05) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calcule le module de distance standardise

    mu = m_B - M_B + alpha * x1 - beta * C

    Parametres
    ----------
    data : dict
        Donnees chargees par load_jla_data ou load_pantheon_data
    alpha : float
        Coefficient d'etirement (defaut: 0.141, Betoule 2014)
    beta : float
        Coefficient de couleur (defaut: 3.101, Betoule 2014)
    M_B : float
        Magnitude absolue de reference (defaut: -19.05)

    Retourne
    --------
    mu : array
        Modules de distance standardises
    sigma_mu : array
        Erreurs sur les modules de distance

    Notes
    -----
    Les parametres alpha, beta, M_B sont fixes aux valeurs best-fit
    de Betoule et al. (2014) pour Lambda-CDM, permettant une
    comparaison directe avec le modele JANUS.
    """
    # Module de distance
    mu = data['mb'] - M_B + alpha * data['x1'] - beta * data['color']

    # Propagation des erreurs (approximation diagonale)
    sigma_mu = np.sqrt(
        data['dmb']**2 +
        (alpha * data['dx1'])**2 +
        (beta * data['dcolor'])**2
    )

    return mu, sigma_mu


def load_covariance_matrix(cov_dir: str = 'data/jla/covmat',
                           n_sne: int = 740) -> np.ndarray:
    """
    Charge la matrice de covariance complete pour JLA

    Parametres
    ----------
    cov_dir : str
        Repertoire contenant les fichiers de covariance
    n_sne : int
        Nombre de supernovae

    Retourne
    --------
    cov : array (n_sne x n_sne)
        Matrice de covariance complete
    """
    cov_dir = Path(cov_dir)

    # Fichiers de covariance a sommer
    cov_files = [
        'C_stat.dat',
        'C_cal.dat',
        'C_model.dat',
        'C_bias.dat',
        'C_dust.dat',
        'C_host.dat',
        'C_nonia.dat',
        'C_pecvel.dat'
    ]

    cov = np.zeros((n_sne, n_sne))

    for fname in cov_files:
        fpath = cov_dir / fname
        if fpath.exists():
            try:
                cov_part = np.loadtxt(fpath)
                if cov_part.shape == (n_sne, n_sne):
                    cov += cov_part
                else:
                    # Format flatten
                    cov += cov_part.reshape(n_sne, n_sne)
            except Exception as e:
                warnings.warn(f"Erreur lecture {fname}: {e}")
        else:
            warnings.warn(f"Fichier non trouve: {fpath}")

    return cov


def load_pantheon_covariance(filepath: str = 'data/pantheon/Pantheon+SH0ES_STAT+SYS.cov',
                              verbose: bool = True) -> np.ndarray:
    """
    Charge la matrice de covariance Pantheon+ (stat + sys)

    Format du fichier:
    - Ligne 1 : dimension N (1701)
    - Lignes suivantes : N*N valeurs de la matrice aplatie

    Parametres
    ----------
    filepath : str
        Chemin vers le fichier .cov
    verbose : bool
        Afficher les informations de chargement

    Retourne
    --------
    cov : array (N x N)
        Matrice de covariance complete
    """
    path = Path(filepath)

    if not path.exists():
        raise FileNotFoundError(
            f"Fichier covariance non trouve: {filepath}\n"
            "Telechargez depuis: https://github.com/PantheonPlusSH0ES/DataRelease"
        )

    # Lecture du fichier
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Premiere ligne = dimension
    n_dim = int(lines[0].strip())

    # Lire toutes les valeurs suivantes
    values = []
    for line in lines[1:]:
        values.extend([float(x) for x in line.strip().split()])

    # Verifier le nombre de valeurs
    expected = n_dim * n_dim
    if len(values) != expected:
        raise ValueError(
            f"Nombre de valeurs incorrect: {len(values)} vs {expected} attendues"
        )

    # Reshape en matrice carree
    cov = np.array(values).reshape(n_dim, n_dim)

    if verbose:
        print(f"Matrice covariance Pantheon+ chargee: {n_dim}x{n_dim}")
        print(f"Diagonale: [{cov[0,0]:.4f}, ..., {cov[-1,-1]:.4f}]")
        print(f"Trace: {np.trace(cov):.2f}")

    return cov


def load_pantheon_with_covariance(data_path: str = 'data/pantheon/Pantheon+SH0ES.dat',
                                   cov_path: str = 'data/pantheon/Pantheon+SH0ES_STAT+SYS.cov',
                                   use_shoes: bool = True,
                                   verbose: bool = True) -> Dict:
    """
    Charge Pantheon+ avec matrice de covariance complete

    Parametres
    ----------
    data_path : str
        Chemin vers le fichier de donnees
    cov_path : str
        Chemin vers le fichier de covariance
    use_shoes : bool
        Utiliser MU_SH0ES (modules de distance calibres) au lieu de mB
    verbose : bool
        Afficher les informations

    Retourne
    --------
    data : dict
        Dictionnaire avec donnees + covariance
    """
    path = Path(data_path)

    if not path.exists():
        raise FileNotFoundError(f"Fichier non trouve: {data_path}")

    # Charger les donnees brutes
    df = pd.read_csv(data_path, delim_whitespace=True, comment='#')

    # Verifier que MU_SH0ES est present si demande
    if use_shoes and 'MU_SH0ES' not in df.columns:
        raise ValueError("Colonne MU_SH0ES non trouvee dans le fichier")

    # Charger la covariance
    cov = load_pantheon_covariance(cov_path, verbose=verbose)

    # Verifier coherence dimensions
    if len(df) != cov.shape[0]:
        raise ValueError(
            f"Incoherence: {len(df)} SNe vs covariance {cov.shape[0]}x{cov.shape[0]}"
        )

    # Extraire les erreurs diagonales de la covariance
    sigma_mu = np.sqrt(np.diag(cov))

    result = {
        'z': df['zHD'].values,
        'mu': df['MU_SH0ES'].values if use_shoes else None,
        'sigma_mu': sigma_mu,
        'mb': df['mB'].values if 'mB' in df.columns else None,
        'x1': df['x1'].values if 'x1' in df.columns else None,
        'color': df['c'].values if 'c' in df.columns else None,
        'host_mass': df['HOST_LOGMASS'].values if 'HOST_LOGMASS' in df.columns else None,
        'covariance': cov,
        'n_sne': len(df),
        'names': df['CID'].values if 'CID' in df.columns else np.arange(len(df)),
        'dataset': 'Pantheon+',
        'reference': 'Brout et al. (2022), Scolnic et al. (2022)'
    }

    if verbose:
        print(f"\nDataset Pantheon+ charge: {result['n_sne']} supernovae")
        print(f"Plage redshift: z = [{result['z'].min():.4f}, {result['z'].max():.4f}]")
        if use_shoes:
            print(f"Plage mu (MU_SH0ES): [{result['mu'].min():.2f}, {result['mu'].max():.2f}]")
        print(f"Covariance: {cov.shape[0]}x{cov.shape[1]}, rang = {np.linalg.matrix_rank(cov)}")

    return result


def compute_chi2_with_covariance(mu_obs: np.ndarray,
                                  mu_theory: np.ndarray,
                                  cov: np.ndarray,
                                  offset: float = 0.0) -> float:
    """
    Calcule chi2 avec matrice de covariance complete

    chi2 = (mu_obs - mu_th - offset)^T @ Cov^{-1} @ (mu_obs - mu_th - offset)

    Parametres
    ----------
    mu_obs : array
        Modules de distance observes
    mu_theory : array
        Modules de distance theoriques
    cov : array (N x N)
        Matrice de covariance
    offset : float
        Offset a soustraire (nuisance parameter)

    Retourne
    --------
    chi2 : float
        Valeur du chi2
    """
    residuals = mu_obs - mu_theory - offset

    # Inversion via Cholesky pour stabilite numerique
    try:
        L = np.linalg.cholesky(cov)
        y = np.linalg.solve(L, residuals)
        chi2 = np.dot(y, y)
    except np.linalg.LinAlgError:
        # Fallback: pseudo-inverse si matrice singuliere
        warnings.warn("Matrice singuliere, utilisation pseudo-inverse")
        cov_inv = np.linalg.pinv(cov)
        chi2 = residuals @ cov_inv @ residuals

    return chi2


def filter_data(data: Dict,
                z_min: float = 0.0,
                z_max: float = np.inf,
                quality_cut: bool = True) -> Dict:
    """
    Filtre les donnees selon les criteres specifies

    Parametres
    ----------
    data : dict
        Donnees a filtrer
    z_min, z_max : float
        Limites en redshift
    quality_cut : bool
        Appliquer les coupes de qualite standard

    Retourne
    --------
    filtered : dict
        Donnees filtrees
    """
    mask = (data['z'] >= z_min) & (data['z'] <= z_max)

    if quality_cut:
        # Coupes standard (a ajuster selon dataset)
        mask &= np.isfinite(data['mb'])
        mask &= np.isfinite(data['x1'])
        mask &= np.isfinite(data['color'])
        mask &= data['dmb'] > 0
        mask &= np.abs(data['x1']) < 3  # Etirement raisonnable
        mask &= np.abs(data['color']) < 0.3  # Couleur raisonnable

    filtered = {
        'z': data['z'][mask],
        'mb': data['mb'][mask],
        'x1': data['x1'][mask],
        'color': data['color'][mask],
        'dmb': data['dmb'][mask],
        'dx1': data['dx1'][mask],
        'dcolor': data['dcolor'][mask],
        'n_sne': mask.sum(),
        'names': data['names'][mask] if isinstance(data['names'], np.ndarray) else None,
        'dataset': data.get('dataset', 'Unknown'),
        'reference': data.get('reference', '')
    }

    return filtered


# Tests
if __name__ == "__main__":
    print("Test du module data_loader.py")
    print("=" * 50)

    # Test avec donnees simulees
    n_test = 100
    test_data = {
        'z': np.random.uniform(0.01, 1.0, n_test),
        'mb': np.random.normal(24, 0.5, n_test),
        'x1': np.random.normal(0, 1, n_test),
        'color': np.random.normal(0, 0.1, n_test),
        'dmb': np.abs(np.random.normal(0.1, 0.02, n_test)),
        'dx1': np.abs(np.random.normal(0.1, 0.02, n_test)),
        'dcolor': np.abs(np.random.normal(0.02, 0.005, n_test)),
        'n_sne': n_test,
        'names': np.arange(n_test),
        'dataset': 'Test'
    }

    mu, sigma_mu = compute_distance_modulus(test_data)
    print(f"\nDonnees test: {n_test} SNe")
    print(f"mu moyen: {mu.mean():.2f} +/- {sigma_mu.mean():.3f}")

    # Test filtrage
    filtered = filter_data(test_data, z_min=0.1, z_max=0.5)
    print(f"Apres filtrage: {filtered['n_sne']} SNe")

    print("\n" + "=" * 50)
    print("Tests passes avec succes!")
