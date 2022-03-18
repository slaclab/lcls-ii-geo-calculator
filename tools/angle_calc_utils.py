## utility functions for performing angle calculations
import numpy as np

def bravais_matrix(a, b, c, alpha, beta, gamma):
    # calculate realspace characteristic matrix for a lattice determined by
    # unit cell lengths and angles
    bravais = np.array([[a, 0, 0],
                        [b * np.cos(gamma), b * np.sin(gamma), 0],
                        [c * np.cos(beta), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
                         c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) -
                                     np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2) / np.sin(gamma)]])
    return bravais

def rotationmat3D(r, Axis):
    # creates a rotation matrix such that R * x
    # operates on x by rotating x around the origin r RADIANS around line
    # connecting the origin to the point "Axis"

    #  normalize axis of rotation
    L = np.linalg.norm(Axis)
    if L < np.finfo(float).eps:
        raise ValueError('axis direction must be non-zero vector')
    Axis = Axis / L

    # intermediates
    u = Axis[0]
    v = Axis[1]
    w = Axis[2]
    u2 = u ** 2
    v2 = v ** 2
    w2 = w ** 2
    c = np.cos(r)
    s = np.sin(r)

    R = np.eye(3)
    # assemble rotation matrix
    R[0, 0] = u2 + (v2 + w2) * c
    R[0, 1] = u * v * (1 - c) - w * s
    R[0, 2] = u * w * (1 - c) + v * s
    R[1, 0] = u * v * (1 - c) + w * s
    R[1, 1] = v2 + (u2 + w2) * c
    R[1, 2] = v * w * (1 - c) - u * s
    R[2, 0] = u * w * (1 - c) - v * s
    R[2, 1] = v * w * (1 - c) + u * s
    R[2, 2] = w2 + (u2 + v2) * c
    return R

def euler_rotationmat(phi, theta, chi):
    # returns the 3d rotiation matrix representing rotation of euler angles by x-convention
    # phi about -z -> chi about x -> theta about -z
    # currently borrowing axis conventions from matlab huber_4c

    yaxis = np.transpose([0., 1., 0.])
    zaxis = np.transpose([0., 0., 1.])
    xaxis = np.transpose([1., 0., 0.])

    Mphi = rotationmat3D(phi, -zaxis)
    Mchi = rotationmat3D(chi, xaxis)
    Mtheta = rotationmat3D(theta, -zaxis)

    euler_rot = np.dot(np.dot(Mtheta, Mchi), Mphi)
    # euler_rot = Mtheta * Mchi * Mphi
    return euler_rot



def huber_sixcircle_matrix(phi, theta, chi, mu):
    # todo function block comment

    yaxis = np.transpose([0., 1., 0.])
    zaxis = np.transpose([0., 0., 1.])
    xaxis = np.transpose([1., 0., 0.])

    # order of rotations: phi(y), chi(x), eta(y)

    Mtheta = rotationmat3D(theta, -yaxis)
    Mchi = rotationmat3D(chi, xaxis)
    Mphi = rotationmat3D(phi, -yaxis)
    Mmu = rotationmat3D(mu, zaxis)
    rot = Mmu @ Mtheta @ Mchi @ Mphi
    return rot


def huber_compare(angles, r_final):
    # maybe need to move this out b/c it's for a call to scipy optimize
    phi, theta, chi = angles
    return (np.linalg.norm(huber_sixcircle_matrix(phi, theta, chi, 0) - r_final))



