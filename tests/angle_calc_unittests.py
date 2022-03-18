# unit tests for various angle calculations
from KappaDiffractometer import KappaDiffractometer
from tools.angle_calc_utils import *



def kappa_rot_test():
    # tests the kappa rotation matrix from angle calc utils
    #

    kdiff = KappaDiffractometer()

    # test case 1
    # test conversion done with https://7id.xray.aps.anl.gov/calculators/kappa.html
    #phi, kappa, omega = -4.21, 13.066, -4.21
    #eu_theta, eu_chi, eu_phi = 0., 10., 0.

    phi, kappa, omega = np.deg2rad(25.721), np.deg2rad(32.824),  np.deg2rad(77.721)
    eu_theta, eu_chi, eu_phi = np.deg2rad(15), np.deg2rad(25), np.deg2rad(67)

    kappa_rotmat = kdiff.kappa_rotation_matrix(phi, kappa, omega)
    ref_rotmat = euler_rotationmat(eu_phi, eu_theta, eu_chi)

    err = np.linalg.norm(kappa_rotmat - ref_rotmat)
    print('first test case for kappa rotation matrix - calculated R_kappa : ')
    print(kappa_rotmat)
    print('calculated R_euler :')
    print(ref_rotmat)
    print('L2 norm error is ', err)

    return err

def kappa_vs_eulerian_si():
    kdiff = KappaDiffractometer()

    kphi = np.deg2rad(-65.163)
    kkap = np.deg2rad(32.824)
    komg = np.deg2rad(-100.721)

    R_kappa = kdiff.kappa_rotation_matrix(kphi, kkap, komg)

    ephi = np.deg2rad(90.000)
    etheta = np.deg2rad(58.480)
    echi = np.deg2rad(-94.763)

    R_euler = euler_rotationmat(ephi, etheta, echi)

    print ("\n comparing Kappa diffractometer + si result with Karthik's Euler diff + si")
    print('found kappa angles phi, kappa, omega: ', [kphi, kkap, komg])
    print('rotation matrix:')
    print(R_kappa)
    print("from Karthik, found euler angles phi, theta, chi: ", [ephi, etheta, echi])
    print("rotation matrix:")
    print(R_euler)

    err = np.linalg.norm(R_kappa - R_euler)
    print("error = ", err)
