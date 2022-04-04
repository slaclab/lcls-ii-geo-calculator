## class definition file for 'diffractometer' obj
# hold/set instrument angles etc
# convert between eulerian sample position and kappa diffractometer angles
# values in degrees for readability
import scipy
from scipy import optimize
from tools.angle_calc_utils import *

class KappaDiffractometer:

    def __init__(self):
        # initialize all angles to zero. ANGLES IN RADIANS
        self.phi = 0  # base stage
        self.kappa = 0  # middle stage
        self.omega = 0 # upper stage
        self.kappa_arm_angle = np.deg2rad(50) # angle between phi and kappa axes

        self.nu = 0 # detector polar angle (sweep)
        self.delta = 0 # detector azimuth angle (tilt) (DiffracSample class handles constraints on this during optimization)

        self.max_delta = np.deg2rad(2) # magnitude, degrees; max achievable delta angle
        self.max_kappa = np.deg2rad(100) # magnitude, degrees; max achievable kappa angle

        self.ref_rotation = None # store full 3dspace rotation matrix for current 'state'

    def set_angles(self, angle_array, radians=True):
        # set kappa diffractometer using a list of angle values in radians [phi, kappa, omega, theta]
        # if radians False, assume angle values in degrees
        if radians:
            self.phi, self.kappa, self.omega = angle_array
        else:
            self.phi, self.kappa, self.omega = np.deg2rad(angle_array)

    def set_kappa_arm_angle(self, angle, radians=True):
        # set the angle of the kappa arm
        # this is fixed in equipment
        if radians:
            self.kappa_arm_angle = angle
        else:
            self.kappa_arm_angle = np.deg2rad(angle)



    def set_kappa_from_eulerian(self, reference_rotation):
        # given a goal endpoint in the form of a rotation matrix,
        # calculate the stage angles of the kappa diffractometer
        # obey instrument constraints of max available kappa

        def kappa_compare(params):
            # objective function compares the 3d rotation angle achieved via instrument angles
            # with the desired 3d rotation matrix
            phi, kappa, omega = params
            # objective value is norm^2 error
            return np.linalg.norm(self.kappa_rotation_matrix(phi, kappa, omega) - reference_rotation) ** 2

        # constrain on kappa:
        def kappa_constraint_fn(params):
            kappa = params[1]
            return self.max_kappa - abs(kappa)
        # currently not using constraint on kappa
        cons = {'type': 'ineq', 'fun': kappa_constraint_fn}

        kappa_solution = scipy.optimize.minimize(kappa_compare, np.array([1,1,1]), tol=1e-20).x
        self.phi = kappa_solution[0]
        self.kappa = kappa_solution[1]
        self.omega = kappa_solution[2]

        print("Kappa solution: ", kappa_solution)


    def kappa_rotation_matrix(self, phi, kappa, omega):
        # input candidate kappa stage angles phi (base), kappa (mid), omega (top)
        # return standard rotation matrix to apply to sample from diffractometer position
        # need to check consistency with coordinate system

        # from http://www.nonius.nl/cad4/manuals/user/chapter02.html, rotation matrix is equivalent ot
        # Z(Omega) * Y(-alpha) * Z(Kappa) * Y(+alpha) * Z(Phi)

        yaxis = np.transpose([0., 1., 0.])
        zaxis = np.transpose([0., 0., 1.])
        xaxis = np.transpose([1., 0., 0.])

        kappa_axis = np.dot(rotationmat3D(self.kappa_arm_angle, -xaxis), yaxis)
        R_kappa = rotationmat3D(kappa, kappa_axis)

        #arm_R1 = rotationmat3D(-self.kappa_arm_angle, yaxis)
        #arm_R2 = rotationmat3D(self.kappa_arm_angle, yaxis)
        #R_kappa = np.dot(arm_R1, np.dot(rotationmat3D(kappa, zaxis), arm_R2))

        R_omega = rotationmat3D(omega, -yaxis) # note changed sign for positive rotation about y
        R_phi = rotationmat3D(phi, -yaxis)

        # phi -> kappa -> omega
        rotmat = np.dot(R_omega, np.dot(R_kappa, R_phi))
        # todo test phi axis rotation
        # test - when theta zero, phi and omega axis the same

        # no longer including detector rotation - make sure to include elsewhere:
        ##rotmat = np.dot(rotationmat3D(self.theta, zaxis), rotmat)
        return rotmat

    def calc_detector_angles(self, sample_kp, sample_lambda0):

        self.nu = np.arctan2(sample_kp[1], sample_kp[0]) # radians
        self.delta = np.arcsin(sample_kp[2] * sample_lambda0)
        #return self.nu, self.delta