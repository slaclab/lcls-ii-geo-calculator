##  Class definition for LCLS-ii test sample ##

import scipy
from scipy import optimize
from KappaDiffractometer import KappaDiffractometer
from tools.angle_calc_utils import *
from tools.conversions import *
from tools.cif_file_tools import *


class DiffracSample:
    # ALL CALCULATIONS SHALL USE RADIANS
    # degrees only allowed at points in user interface

    def __init__(self):
        self.beam_axis = np.array([1, 0, 0])  # beam at x-axis unless overridden
        self.normal_vect = None
        self.inplane_vect = None
        self.lambda0 = 15.4  # xray wavelength (set using beam energy) (todo check units)
        self.k0 = self.get_k0()  # Incident xray wave-vector
        self.at = None  # np.array([[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]]) # realspace lattice vectors
        self.bg = None  # np.linalg.inv(self.at) # reciprocal space lattice vectors
        self.bragg_hkl = None
        self.lattice_filepath = None
        self.sg = None
        self.sample_rotation_angle = 0  # sample rotation in z when held in goniometer, all instrument angles at zero
        self.sample_rmat = np.identity(3)  # sample rotation matrix
        self.alpha0 = None  # incident angle between wavevector and sample normal (NOT BRAGG PLANE)
        self.Ghkl = None  # reciprocal/rotated hkl
        # unit cell lengths [a,b,c] in Angstrom
        self.ucell_a = None
        self.ucell_b = None
        self.ucell_c = None
        # unit cell angles [alpha, beta, gamma] in radians
        self.ucell_alpha = None
        self.ucell_beta = None
        self.ucell_gamma = None

    def readCIF(self, path):

        self.lattice_filepath = path
        try:
            cif_dict = man_read_cif_to_dict(path)
            # save realspace lattice vectors
            self.at, unit_cell_lengths, unit_cell_angles = cif_symm_to_bravais_lattice(cif_dict)
            [self.ucell_a, self.ucell_b, self.ucell_c] = unit_cell_lengths
            [self.ucell_alpha, self.ucell_beta, self.ucell_gamma] = unit_cell_angles
            self.bg = np.linalg.inv(self.at)  # calc reciprocal space lattice vectors

        except Exception:
            print('error - could not read cif file')

    def readVASP(self, path):
        self.lattice_filepath = path
        print("VASP reader not implemented!")

    def try_update_lattice(self):
        # this function is called when a single lattice param is updated
        # (e.g. via gui); should update internal variables only after
        # the full set of numbers are specified (i.e. do not update if only
        # some params have been input at a given time; wait for full set)

        necessary_variables = [self.ucell_a, self.ucell_b, self.ucell_c,
                               self.ucell_alpha, self.ucell_beta, self.ucell_gamma]

        if None in necessary_variables or 0 in necessary_variables:
            return
        else:
            self.at = bravais_matrix(self.ucell_a, self.ucell_b, self.ucell_c,
                                     self.ucell_alpha, self.ucell_beta, self.ucell_gamma)
            self.bg = np.linalg.inv(self.at)

    def find_rotation_grazing_incidence(self, alpha, x0):
        # find_rotation_grazing_incidence finds sample and outgoing wave-vector angles.
        # Finds the rotation matrix to satisfy the Bragg condition at an angle of
        # incidence (alpha) with respect to the sample surface (self.normal_vect)

        lambda0 = self.lambda0
        at = self.at
        bg = self.bg
        sam_rot = self.sample_rmat

        # vectors in reciprocal space
        Ghkl = sam_rot @ bg @ self.bragg_hkl
        self.Ghkl = Ghkl
        k0 = self.k0

        # vectors in real space
        sample_normal = sam_rot @ self.normal_vect

        # angle between incident direction and sample normal
        alph1 = np.arccos(
            np.dot(sample_normal, k0) / (np.linalg.norm(k0) * np.linalg.norm(sample_normal)))

        # rotate around this axis to satisfy the incidence condition
        incidence_rot_axis = np.cross(k0, sample_normal)
        M_incidence = rotationmat3D(-alph1 + np.pi / 2 - alpha, incidence_rot_axis)
        # TODO check sign on alpha

        # and the rotated vectors in the grazing configuration are:
        rotGhkl = M_incidence @ Ghkl
        rotNsample = M_incidence @ sample_normal

        # Find rotation for Bragg condition
        # zeros of this function are the desired condition
        def obj(ph):
            MBragg = rotationmat3D(ph, rotNsample)
            return (np.linalg.norm(k0 + MBragg @ rotGhkl) - np.linalg.norm(k0))

        res = scipy.optimize.fsolve(obj, x0, xtol=1e-20)
        print(obj(res))
        print(res)

        # 3d rotation matrix before adjustment for delta
        semifinal_rotation = np.dot(rotationmat3D(res, rotNsample), M_incidence)

        # adjust rotation for minimum delta angle
        def detector_angle(addtl_beam_rot_angle):
            # rotation about the axis of the incoming beam maintains same experimental result
            # but modifies delta
            adj_rotation_mat = rotationmat3D(addtl_beam_rot_angle, self.beam_axis) @ semifinal_rotation
            kp = np.dot(adj_rotation_mat, Ghkl) + k0
            delt = np.arcsin(kp[2] * self.lambda0)

            return np.abs(delt)

        # TODO add formal bounds based on max delta angle
        # and throw error if not met?

        angle_adjustment = scipy.optimize.minimize(detector_angle, np.pi / 10)
        print("delta angle ", angle_adjustment.fun * 180 / np.pi)

        # get the final rotation
        return rotationmat3D(angle_adjustment.x, self.beam_axis) @ semifinal_rotation

    def find_rotation_grazing_exit(self, exit_angle, x0):
        #   For a given set of miller indices 'hkl', and an exit angle,
        #   find the rotation in 3d space s.t.
        #   the Bragg reflection is satisfied

        # HKL, k0 vectors in reciprocal space
        Ghkl = self.sample_rmat @ self.bg @ self.bragg_hkl
        k0 = 1 / self.lambda0 * self.beam_axis

        # bragg angle
        thetaB = np.arcsin(self.lambda0 * np.linalg.norm(Ghkl) / 2)

        # scattering angle
        twothetaB = 2 * thetaB

        # angle between Ghkl, incident direction
        k0Ghkl = np.dot(k0, Ghkl)
        gamma1 = np.arccos(k0Ghkl / np.linalg.norm(k0Ghkl))

        ## Bragg condition for tentative initial config:
        # rotate around axis normal to Ghkl and k0
        newaxis = np.cross(k0, Ghkl)

        # bragg condition rotation matrix
        MBragg = rotationmat3D(-gamma1 + np.pi / 2 - thetaB, newaxis)

        # rotate intermediate vectors
        rotGhkl = np.dot(MBragg, Ghkl)
        rot_sampleN = np.dot(MBragg, np.dot(self.sample_rmat, self.normal_vect))
        rot_sampleN_norm = rot_sampleN / np.linalg.norm(rot_sampleN)

        # scattered vector
        kp = rotGhkl + k0

        ## rotate to satisfy exit angle requirement
        # zero of this function is the desired angle
        def obj(ph):
            MM = rotationmat3D(ph, rotGhkl)
            rot_n = np.dot(MM, rot_sampleN_norm)
            return np.sin(exit_angle) - np.dot(kp, rot_n) / np.linalg.norm(kp)

        res = scipy.optimize.fsolve(obj, x0, xtol=1e-20)
        print(obj(res))
        print(res)

        # TODO add error checking (if ph=NaN in result, throw e)

        # 3d rotation matrix before adjustment for delta
        semifinal_rotation = np.dot(rotationmat3D(res, rot_sampleN_norm), MBragg)

        # adjust rotation for minimum delta angle
        def detector_angle(addtl_beam_rot_angle):
            # rotation about the axis of the incoming beam maintains same experimental result
            # but modifies delta

            adj_rotation_mat = np.dot(rotationmat3D(addtl_beam_rot_angle, self.beam_axis), semifinal_rotation)
            kp = np.dot(adj_rotation_mat, Ghkl) + k0
            delt = np.arcsin(kp[2] * self.lambda0)
            return np.abs(delt)

        # TODO add formal bounds based on max delta angle
        # and throw error if not met?

        angle_adjustment = scipy.optimize.minimize(detector_angle, 0.001, bounds=[(-np.pi, np.pi)])
        print("delta angle ", angle_adjustment.fun)

        # get the final rotation
        return rotationmat3D(angle_adjustment.x, self.beam_axis) @ semifinal_rotation

    def set_sample_normal(self, normal):
        # allows to set the sample normal vector and auto-calculates a perpendicular in-plane vector

        # normalize the given vector
        self.normal_vect = normal / np.linalg.norm(normal)

        if normal[0] == 0:
            self.inplane_vect = np.array([1, 0, 0])
        elif normal[1] == 0:
            self.inplane_vect = np.array([0, 1, 0])
        elif normal[2] == 0:
            self.inplane_vect = np.array([0, 0, 1])
        else:
            self.inplane_vect = np.array([1, 1, -1 * (normal[0] + normal[1]) / normal[2]])

    def set_beam_energy(self, E):
        # convert beam energy [eV] to wavelength [A]
        self.lambda0 = photonE_to_wavelen_A(E)

    def set_sample_rotation(self, angle=0, axis=np.array([0, 0, 1])):
        self.sample_rotation_angle = angle
        self.sample_rmat = rotationmat3D(angle, axis)

    def get_k0(self):
        return 1 / self.lambda0 * self.beam_axis