import numpy as np

## globals ##
c = 2.998e8  # [m/s], speed of light
h = 6.626e-34  # [J*s], planck's constant
m_to_nm = 1e9  # multiply to convert m to nm
nm_to_m = 1 / m_to_nm  # multiply to convert nm to m
nm_to_A = 10 # multiply to convert nm to Angstrom
A_to_nm = 1 / nm_to_A # multiply to convert Angstrom to nm
J_to_eV = 6.241509e18 # multiply to convert J to eV
eV_to_J = 1 / J_to_eV # multiply to convert eV to J


def photonE_to_wavelen_A(energy_eV):
    # lambda = hc/E [returns value in Angstrom]
    energy_J = energy_eV * eV_to_J
    lamb = (h * c / energy_J) * m_to_nm * nm_to_A
    return lamb

def photonE_to_wavelen_nm(energy_eV):
    energy_J = energy_eV * eV_to_J
    lamb = (h * c / energy_J) * m_to_nm
    return lamb

def wavelen_A_to_photonE(lamb_A):
    # E = hc/lambda [returns value in eV]
    E_J = h * c / (lamb_A * A_to_nm * nm_to_m)
    E_eV = E_J * J_to_eV
    return E_eV

def wavelen_nm_to_photonE(lamb_nm):
    # returns value in eV
    E_J = h * c / (lamb_nm * nm_to_m)
    E_eV = E_J * J_to_eV
    return E_eV