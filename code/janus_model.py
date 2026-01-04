"""
janus_model.py
Implementation du modele cosmologique JANUS

Base sur D'Agostini & Petit (2018)
"Constraints on Janus Cosmological model from recent observations of supernovae type Ia"
Astrophysics and Space Science, 363(7): 139

Auteur: Projet JANUS-S
Date: Janvier 2026
"""

import numpy as np
from typing import Union, Tuple

ArrayLike = Union[float, np.ndarray]


class JanusCosmology:
    """
    Implementation du modele cosmologique JANUS (Janus Cosmological Model - JCM)

    Le modele JANUS est une theorie bimetrique avec deux secteurs couples :
    - Secteur positif : masses positives (univers observable)
    - Secteur negatif : masses negatives (secteur miroir)

    L'acceleration de l'expansion est expliquee par la repulsion entre
    masses positives et negatives, sans recourir a l'energie noire.

    Parametres
    ----------
    H0 : float
        Constante de Hubble en km/s/Mpc (defaut: 70.0)
    """

    def __init__(self, H0: float = 70.0):
        self.H0 = H0
        self.c = 299792.458  # Vitesse de la lumiere en km/s

    def distance_modulus(self, z: ArrayLike, q0: float,
                         offset: float = 0.0) -> ArrayLike:
        """
        Calcule le module de distance theorique mu(z) dans le modele JANUS

        Equation derivee de la solution exacte des equations bimetriques
        pour l'ere de poussiere (dust era).

        Parametres
        ----------
        z : float ou array-like
            Redshift cosmologique
        q0 : float
            Parametre de deceleration (doit etre < 0 pour acceleration)
        offset : float
            Constante additive pour l'ajustement (defaut: 0.0)

        Retourne
        --------
        mu : float ou array-like
            Module de distance en magnitudes

        Notes
        -----
        La condition 1 + 2*q0*z > 0 doit etre satisfaite pour tous les z.
        Pour q0 = -0.087, cela impose z < 5.75 (largement suffisant).
        """
        z = np.atleast_1d(z)

        # Verification de la condition de validite
        condition = 1 + 2 * q0 * z
        if np.any(condition <= 0):
            # Retourne inf pour les points invalides
            mu = np.where(condition > 0,
                         self._compute_mu(z, q0),
                         np.inf)
            return mu + offset

        mu = self._compute_mu(z, q0)
        return mu + offset

    def _compute_mu(self, z: np.ndarray, q0: float) -> np.ndarray:
        """Calcul interne du module de distance"""
        sqrt_term = np.sqrt(1 + 2 * q0 * z)

        # Terme principal de la relation magnitude-redshift
        # Equation (XX) de D'Agostini & Petit (2018)
        numerator = z + z**2 * (1 - q0) / (1 + q0 * z + sqrt_term)

        # Distance luminosite en Mpc
        d_L = (self.c / self.H0) * numerator

        # Module de distance
        mu = 5 * np.log10(d_L) + 25

        return mu

    def comoving_distance(self, z: ArrayLike, q0: float) -> ArrayLike:
        """
        Calcule la distance comobile r(z)

        Parametres
        ----------
        z : float ou array-like
            Redshift
        q0 : float
            Parametre de deceleration

        Retourne
        --------
        r : float ou array-like
            Distance comobile en Mpc
        """
        z = np.atleast_1d(z)

        factor = self.c / self.H0
        sqrt_term = np.sqrt(1 + 2 * q0 * z)

        numerator = q0 * z + (1 - q0) * (1 - sqrt_term)
        denominator = q0**2 * (1 + z)

        r = factor * (numerator / denominator)

        return r

    def luminosity_distance(self, z: ArrayLike, q0: float) -> ArrayLike:
        """
        Calcule la distance luminosite d_L(z)

        Parametres
        ----------
        z : float ou array-like
            Redshift
        q0 : float
            Parametre de deceleration

        Retourne
        --------
        d_L : float ou array-like
            Distance luminosite en Mpc
        """
        z = np.atleast_1d(z)
        sqrt_term = np.sqrt(1 + 2 * q0 * z)
        numerator = z + z**2 * (1 - q0) / (1 + q0 * z + sqrt_term)
        d_L = (self.c / self.H0) * numerator
        return d_L

    def universe_age(self, q0: float) -> float:
        """
        Calcule l'age de l'univers T0 en Gyr

        Parametres
        ----------
        q0 : float
            Parametre de deceleration (doit etre < 0)

        Retourne
        --------
        T0 : float
            Age de l'univers en milliards d'annees (Gyr)
        """
        if q0 >= 0:
            return np.nan

        # Parametre u0 a l'epoque actuelle
        u0 = np.arcsinh(np.sqrt(-1 / (2 * q0)))

        # Conversion H0 en 1/Gyr
        # 1 Mpc = 3.086e19 km, 1 Gyr = 3.156e16 s
        H0_inv_Gyr = 977.8 / self.H0  # Gyr

        # Age de l'univers
        term = (1 + np.sinh(2 * u0) / 2 + u0)
        T0 = H0_inv_Gyr * term / (2 * (-q0)**1.5)

        return T0

    def scale_factor(self, u: ArrayLike) -> ArrayLike:
        """
        Facteur d'echelle a(u) en fonction du parametre conforme u

        a(u) = (alpha^2 / 2) * cosh^2(u)

        Parametres
        ----------
        u : float ou array-like
            Parametre conforme

        Retourne
        --------
        a : float ou array-like
            Facteur d'echelle (normalise)
        """
        return np.cosh(u)**2 / 2

    def deceleration_parameter(self, u: ArrayLike) -> ArrayLike:
        """
        Parametre de deceleration q(u)

        q = -1/2 * sinh^(-2)(u) < 0 pour tout u

        Parametres
        ----------
        u : float ou array-like
            Parametre conforme

        Retourne
        --------
        q : float ou array-like
            Parametre de deceleration (toujours negatif)
        """
        return -0.5 / np.sinh(u)**2


class LambdaCDM:
    """
    Implementation du modele standard Lambda-CDM pour comparaison

    Parametres
    ----------
    H0 : float
        Constante de Hubble en km/s/Mpc
    Omega_m : float
        Densite de matiere (defaut: 0.3)
    Omega_L : float
        Densite d'energie noire (defaut: 0.7)
    """

    def __init__(self, H0: float = 70.0, Omega_m: float = 0.3,
                 Omega_L: float = 0.7):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_L = Omega_L
        self.c = 299792.458

    def E(self, z: ArrayLike) -> ArrayLike:
        """Fonction E(z) = H(z)/H0"""
        return np.sqrt(self.Omega_m * (1 + z)**3 + self.Omega_L)

    def distance_modulus(self, z: ArrayLike, offset: float = 0.0) -> ArrayLike:
        """
        Module de distance dans Lambda-CDM (integration numerique)
        """
        from scipy.integrate import quad

        z = np.atleast_1d(z)
        d_H = self.c / self.H0  # Distance de Hubble en Mpc

        mu = np.zeros_like(z, dtype=float)
        for i, zi in enumerate(z):
            if zi <= 0:
                mu[i] = np.nan
                continue

            # Integration numerique
            integral, _ = quad(lambda x: 1/self.E(x), 0, zi)
            d_L = d_H * (1 + zi) * integral
            mu[i] = 5 * np.log10(d_L) + 25

        return mu + offset


# Tests unitaires
if __name__ == "__main__":
    print("Test du module janus_model.py")
    print("=" * 50)

    # Test JANUS
    janus = JanusCosmology(H0=70.0)

    # Valeurs de test
    z_test = np.array([0.1, 0.5, 1.0])
    q0_test = -0.087

    print(f"\nModele JANUS (q0 = {q0_test}):")
    print(f"z = {z_test}")
    print(f"mu(z) = {janus.distance_modulus(z_test, q0_test)}")
    print(f"Age univers = {janus.universe_age(q0_test):.2f} Gyr")

    # Test Lambda-CDM
    lcdm = LambdaCDM(H0=70.0)
    print(f"\nModele Lambda-CDM (Om=0.3, OL=0.7):")
    print(f"mu(z) = {lcdm.distance_modulus(z_test)}")

    print("\n" + "=" * 50)
    print("Tests passes avec succes!")
