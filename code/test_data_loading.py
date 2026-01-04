#!/usr/bin/env python3
"""
test_data_loading.py
Test de validation du chargement des datasets JLA et Pantheon+

Execute: python code/test_data_loading.py
"""

import sys
import os

# Ajouter le repertoire code au path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd


def test_jla():
    """Test chargement JLA"""
    print("\n" + "=" * 60)
    print("TEST CHARGEMENT JLA")
    print("=" * 60)

    filepath = 'data/jla/jla_likelihood_v6/data/jla_lcparams.txt'

    if not os.path.exists(filepath):
        print(f"ERREUR: Fichier non trouve: {filepath}")
        return False

    try:
        data = pd.read_csv(filepath, delim_whitespace=True)

        print(f"Fichier: {filepath}")
        print(f"Nombre de supernovae: {len(data)}")
        print(f"Colonnes: {list(data.columns[:10])}...")

        # Verifications
        assert len(data) == 740, f"Attendu 740 SNe, trouve {len(data)}"
        assert 'zcmb' in data.columns, "Colonne 'zcmb' manquante"
        assert 'mb' in data.columns, "Colonne 'mb' manquante"
        assert 'x1' in data.columns, "Colonne 'x1' manquante"
        assert 'color' in data.columns, "Colonne 'color' manquante"

        # Statistiques
        print(f"\nStatistiques redshift (z):")
        print(f"  Min: {data['zcmb'].min():.4f}")
        print(f"  Max: {data['zcmb'].max():.4f}")
        print(f"  Moyenne: {data['zcmb'].mean():.4f}")

        print(f"\nStatistiques magnitude (mb):")
        print(f"  Min: {data['mb'].min():.2f}")
        print(f"  Max: {data['mb'].max():.2f}")
        print(f"  Moyenne: {data['mb'].mean():.2f}")

        print("\n✅ TEST JLA: SUCCES")
        return True

    except Exception as e:
        print(f"ERREUR: {e}")
        return False


def test_pantheon():
    """Test chargement Pantheon+"""
    print("\n" + "=" * 60)
    print("TEST CHARGEMENT PANTHEON+")
    print("=" * 60)

    filepath = 'data/pantheon/Pantheon+SH0ES.dat'

    if not os.path.exists(filepath):
        print(f"ERREUR: Fichier non trouve: {filepath}")
        return False

    try:
        data = pd.read_csv(filepath, delim_whitespace=True)

        print(f"Fichier: {filepath}")
        print(f"Nombre de supernovae: {len(data)}")
        print(f"Colonnes: {list(data.columns[:10])}...")

        # Verifications
        assert len(data) >= 1500, f"Attendu ~1700 SNe, trouve {len(data)}"
        assert 'zHD' in data.columns, "Colonne 'zHD' manquante"
        assert 'm_b_corr' in data.columns, "Colonne 'm_b_corr' manquante"

        # Statistiques
        print(f"\nStatistiques redshift (zHD):")
        print(f"  Min: {data['zHD'].min():.4f}")
        print(f"  Max: {data['zHD'].max():.4f}")
        print(f"  Moyenne: {data['zHD'].mean():.4f}")

        print(f"\nStatistiques magnitude corrigee (m_b_corr):")
        print(f"  Min: {data['m_b_corr'].min():.2f}")
        print(f"  Max: {data['m_b_corr'].max():.2f}")
        print(f"  Moyenne: {data['m_b_corr'].mean():.2f}")

        # Compter les SNe uniques (certaines ont plusieurs observations)
        unique_sne = data['CID'].nunique()
        print(f"\nSupernovae uniques: {unique_sne}")

        print("\n✅ TEST PANTHEON+: SUCCES")
        return True

    except Exception as e:
        print(f"ERREUR: {e}")
        return False


def main():
    """Execution des tests"""
    print("\n" + "#" * 60)
    print("# VALIDATION DU CHARGEMENT DES DONNEES")
    print("# Projet JANUS-S")
    print("#" * 60)

    # Changer vers le repertoire du projet
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.chdir(project_dir)
    print(f"\nRepertoire: {os.getcwd()}")

    # Tests
    jla_ok = test_jla()
    pantheon_ok = test_pantheon()

    # Resume
    print("\n" + "=" * 60)
    print("RESUME")
    print("=" * 60)
    print(f"JLA:      {'✅ OK' if jla_ok else '❌ ECHEC'}")
    print(f"Pantheon+: {'✅ OK' if pantheon_ok else '❌ ECHEC'}")

    if jla_ok and pantheon_ok:
        print("\n🎉 TOUS LES TESTS PASSES - Phase 1 validee!")
        return 0
    else:
        print("\n⚠️  Certains tests ont echoue")
        return 1


if __name__ == "__main__":
    sys.exit(main())
