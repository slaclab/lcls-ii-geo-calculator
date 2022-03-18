import matplotlib.pyplot as plt
import numpy as np
from KappaDiffractometer import KappaDiffractometer
from DiffracSample import DiffracSample
from plotting.plot_solved_system import *
import scipy
from scipy import optimize

## ALL CALCULATIONS IN RADIANS,         ##
## DEGREES ONLY FOR DISPLAY/UI PURPOSES ##

def main():
    diff = KappaDiffractometer()
    sample = DiffracSample()

def testrun_incidence_exit_angles():
    diff = KappaDiffractometer()
    sample = DiffracSample()

    # spoof user inputs
    cif_path = "CIFfiles/AMS_DATA.cif" # reading silicon file info
    sample.readCIF(cif_path)

    # sample normal/inplane vectors
    sample.set_sample_normal(np.array([1, 1, 1]))
    print("Sample normal: ", sample.normal_vect, "\nSample in-plane: ", sample.inplane_vect)


    print('sample at ', sample.sample_rotation_angle, ' about z axis in goniometer fixture')
    alpha = np.radians(5)

    # input bragg plane hkl array
    print('starting bragg HKL indices as [4, 2, 0]')
    sample.bragg_hkl = np.array([4., 2., 0.])

    # calculate the reference rotation matrix for a sample with incidence angle alpha
    ref_rotation = sample.find_rotation_grazing_incidence(alpha, -10)

    spin_deltas(sample, ref_rotation)

    # translate euler rotation to kappa angles
    diff.set_kappa_from_eulerian(ref_rotation)

    # calculate kp, detector angle (?)
    kp = np.dot(ref_rotation, sample.Ghkl) + sample.k0 # check if dot product correct
    nu = np.arctan2(kp[1], kp[0]) * 180/np.pi # azimuthal

    delt = np.arcsin(kp[2] * sample.lambda0)*180/np.pi # Polar detector angle (small)

    # bragg angle
    thetaB = np.arcsin(sample.lambda0 * np.linalg.norm(sample.bg*sample.bragg_hkl/2)) * 180/np.pi
    tth = 2 * thetaB

    print("Calculated kappa diffractometer stage angles using incidence angle alpha=12:")
    print("phi, kappa, omega = ", np.array([diff.phi, diff.kappa, diff.omega]))
    print("HKL ", sample.bragg_hkl)
    print("detector angles nu ", diff.nu, ", delta ", delt)
    plot_diffrac_system(sample, diff, alpha, nu, delt)

    # do the exit angle thing
    beta = np.radians(10)
    exitref_rotation = sample.find_rotation_grazing_exit(beta, 0)
    spin_deltas(sample, exitref_rotation)
    diff.set_kappa_from_eulerian(exitref_rotation)

    kp = np.dot(exitref_rotation, sample.Ghkl) + sample.k0
    kp_norm = kp / np.linalg.norm(kp)
    nu = np.arctan2(kp[1], kp[0]) * 180/np.pi
    #delt = 90 - np.arccos(kp_norm[2] * sample.lambda0) * 180/np.pi
    delt = np.arcsin(kp_norm[2])*180/np.pi # Polar detector angle (small)

    # bragg angle
    thetaB = np.arcsin(sample.lambda0 * np.linalg.norm(sample.bg*sample.bragg_hkl/2)) * 180/np.pi
    tth = 2 * thetaB

    print("Calculated kappa diffractometer stage angles using exit angle beta=20:")
    print("phi, kappa, omega = ", np.array([diff.phi, diff.kappa, diff.omega]))
    print("HKL ", sample.bragg_hkl)
    print("detector angles nu ", diff.nu, ", delta ", delt)
    plot_diffrac_system(sample, diff, alpha, nu, delt)

def spin_deltas(sample, res_rotation):
    # rotates the sample about the beam axis to check for effect on detector angle
    extra_rotations = np.linspace(0, 2*np.pi, 180)
    deltas = np.zeros_like(extra_rotations)
    
    for i in range(len(extra_rotations)):
        r = extra_rotations[i]
        rmat = rotationmat3D(r, sample.beam_axis) @ res_rotation
        kp = np.dot(rmat, sample.Ghkl) + sample.k0  # check if dot product correct
        kp_norm = kp / np.linalg.norm(kp)
        deltas[i] = np.arcsin(kp_norm[2]) * 180/np.pi
        #deltas[i] = 90.0 - np.arccos(kp_norm[2] * sample.lambda0) * 180 / np.pi  # Polar detector angle (small)
    plt.figure()
    plt.plot(np.rad2deg(extra_rotations), deltas)
    plt.title("Detector angle vs additional beam ax rotation")
    plt.show()
if __name__ == "__main__":
    testrun_incidence_exit_angles()