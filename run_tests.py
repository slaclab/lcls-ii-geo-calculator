import numpy as np
from KappaDiffractometer import KappaDiffractometer
from DiffracSample import DiffracSample
from tools.angle_calc_utils import *
from tests.angle_calc_unittests import *
from tests.test_unitcell_plotter import *


def main():
    #kappa_rot_test()
    #kappa_vs_eulerian_si()
    plot_with_kapparot()


if __name__ == '__main__':
    main()