from KappaDiffractometer import KappaDiffractometer
from tools.angle_calc_utils import *
from plotting.plot_lattice import *

def plot_with_kapparot():
    phi, kappa, omega = np.deg2rad(20.34), np.deg2rad(60), np.deg2rad(20.34)
    bravais = np.array([[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]) # si

    plot_lattice_rotation(bravais, phi, kappa, omega)