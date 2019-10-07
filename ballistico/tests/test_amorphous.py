"""
Unit and regression test for the ballistico package.
"""

# Import package, test suite, and other packages as needed
from finitedifference.finitedifference import FiniteDifference
import numpy as np
from ballistico.phonons import Phonons
import ballistico.conductivity as bac
import ase.units as units
import shutil


def create_phonons(tmpdir):
    tmp_folder = str(tmpdir.dirpath())
    print(tmp_folder)
    shutil.rmtree(tmp_folder, ignore_errors=True)

    # Create a finite difference object
    finite_difference = FiniteDifference.import_from_dlpoly_folder(folder='si-amorphous')

    # # Create a phonon object
    phonons = Phonons(finite_difference=finite_difference,
                      is_classic=True,
                      temperature=300,
                      folder=tmp_folder,
                      sigma_in= 0.05 / 4.135,
                      broadening_shape='triangle')
    return phonons


def test_first_gamma(tmpdir):
    phonons = create_phonons(tmpdir)
    THZTOMEV = units.J * units._hbar * 2 * np.pi * 1e15
    np.testing.assert_approx_equal(phonons.gamma[3] * THZTOMEV / (2 * np.pi), 22.451, significant=3)


def test_second_gamma(tmpdir):
    phonons = create_phonons(tmpdir)
    THZTOMEV = units.J * units._hbar * 2 * np.pi * 1e15
    np.testing.assert_approx_equal(phonons.gamma[4] * THZTOMEV / (2 * np.pi), 23.980, significant=3)


def test_qhgk_conductivity(tmpdir):
    phonons = create_phonons(tmpdir)
    cond = bac.conductivity(phonons, method='qhgk').sum(axis=0)
    cond = np.abs(np.mean(cond.diagonal()))
    np.testing.assert_approx_equal(cond, 0.99, significant=2)

