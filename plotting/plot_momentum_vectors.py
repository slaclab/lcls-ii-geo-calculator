import numpy as np
from tools.angle_calc_utils import rotationmat3D
MAGNITUDE_SCALE = 1.5

def plot_momentum_transfer(sample, unitcell_center, ax):
    plot_beamline_in(sample, unitcell_center, ax)
    plot_scattered_radiation(sample, unitcell_center, ax)
    return ax

def plot_beamline_in(sample, unitcell_center, ax):
    # plot xray radiation vector down beam into sample
    max_cell_length = max(sample.ucell_a, sample.ucell_b, sample.ucell_c)
    vmag = MAGNITUDE_SCALE * max_cell_length
    v = sample.beam_axis * vmag
    v1 = unitcell_center - v
    vx = [unitcell_center[0], v1[0]]
    vy = [unitcell_center[1], v1[1]]
    vz = [unitcell_center[2], v1[2]]
    ax.plot(vx, vy, vz)

    return ax

def plot_scattered_radiation(sample, unitcell_center, ax):
    # plot line of scattered xray radiation out
    max_cell_length = max(sample.ucell_a, sample.ucell_b, sample.ucell_c)
    vmag = MAGNITUDE_SCALE * max_cell_length
    v = vmag * (np.linalg.norm(sample.kp))
    v2 = unitcell_center + v / np.linalg.norm(v)
    vx = [unitcell_center[0], v2[0]]
    vy = [unitcell_center[1], v2[1]]
    vz = [unitcell_center[2], v2[2]]
    ax.plot(vx, vy, vz)
    return ax

def plot_detector_target(sample, unitcell_center, nu, delta, ax):
    # plot representation of detector position/orientation for visual inspection
    max_cell_length = max(sample.ucell_a, sample.ucell_b, sample.ucell_c)
    offset = [2 * max_cell_length, 2.5 * max_cell_length]
    detector_inline = np.transpose(np.array([unitcell_center - offset[0]*sample.beam_axis,
                                unitcell_center - offset[1]*sample.beam_axis]))
    # transformation axes assume beam is along lab x axis
    delta_rotation = rotationmat3D(delta, np.array([0,1,0]))
    nu_rotation = rotationmat3D(nu, np.array([0,0,1]))
    detector_inline = np.dot(nu_rotation, np.dot(delta_rotation, detector_inline))
    # detector inline is a sequence of coordinate triplets, not a numpy vector

    detector_zaxis = np.dot(nu_rotation, np.dot(delta_rotation, np.array([0,0,1])))
    detector_centeraxis = np.array(detector_inline[:,1] - detector_inline[:,0])
    detector_yaxis = np.dot(rotationmat3D(np.pi/2, detector_centeraxis), detector_zaxis)
    detector_zaxis = detector_zaxis / np.linalg.norm(detector_zaxis) * 0.5
    detector_yaxis = detector_yaxis / np.linalg.norm(detector_yaxis) * 0.5
    cross_z = [[detector_inline[0,0], detector_inline[0,0]+detector_zaxis[0]],
               [detector_inline[1,0], detector_inline[1,0]+detector_zaxis[1]],
               [detector_inline[2,0], detector_inline[2,0]+detector_zaxis[2]]]
    cross_y = [[detector_inline[0, 0], detector_inline[0, 0] + detector_yaxis[0]],
               [detector_inline[1, 0], detector_inline[1, 0] + detector_yaxis[1]],
               [detector_inline[2, 0], detector_inline[2, 0] + detector_yaxis[2]]]

    ax.plot(detector_inline[0,:], detector_inline[1,:], detector_inline[2,:], color='k')
    ax.plot(cross_y[0], cross_y[1], cross_y[2], color='k')
    ax.plot(cross_z[0], cross_z[1], cross_z[2], color='k')

    return ax

def plot_milleridx_plane(sample, unitcell_center, ax):
    h = sample.bragg_hkl[0]; k = sample.bragg_hkl[1]; l = sample.bragg_hkl[2]
    plane_spacing = [0,0,0]
    miller_indices = [h, k, l]
    cell_lengths = [sample.ucell_a, sample.ucell_b, sample.ucell_c]
    for i in range(3):
        if miller_indices[i] != 0:
            plane_spacing[i] = cell_lengths[i] / miller_indices[i]
    [a_spacing, b_spacing, c_spacing] = plane_spacing

    o_corner = unitcell_center - 0.5 * np.array(cell_lengths)




