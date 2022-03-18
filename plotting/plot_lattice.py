import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.spatial import ConvexHull
from KappaDiffractometer import KappaDiffractometer
from tools.angle_calc_utils import *
from DiffracSample import DiffracSample
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def plot_lattice_rotation(at, phi, kappa, omega, sample_rot=None,sample_rot_axis=np.array([0,0,1]), kappa_rot_mat=None, show=True, hide_orig=False):

    # sample rotation
    if sample_rot is not None:
        rmat = rotationmat3D(sample_rot, sample_rot_axis)
    else:
        rmat = eye(3)

    # instrument rotation
    kdiff = KappaDiffractometer()
    kmat = kdiff.kappa_rotation_matrix(phi, kappa, omega)
    if kappa_rot_mat is not None:
        # have option to override and use a rotation matrix instead of angles
        kmat = kappa_rot_mat

    # final (total) sample rotation matrix
    rmat = np.dot(kmat, rmat)




    # take linear combinations of bravais vectors
    b1, b2, b3 = at[0, :], at[1, :], at[2, :] # using rows todo check
    o = np.zeros(3)
    # get vertices by combining vectors
    # doing two stacked unit cells
    v = [o, b1, b2, b1+b2, b3, b1+b3, b2+b3, b1+b2+b3, b3*2, b1+b3*2, b2+b3*2, b1+b2+b3*2]
    vrot = [vi.dot(rmat) for vi in v]
    coords = np.array(v)
    rot_coords = np.array(vrot)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # plot surfaces parallel to b1, b2 vectors and perpendicular to b3
    for j in [0, 4, 8]:
        x, y, z = coords[j:j+4, 0], coords[j:j+4, 1], coords[j:j+4, 2]
        xr, yr, zr = rot_coords[j:j+4, 0], rot_coords[j:j+4, 1], rot_coords[j:j+4, 2]
        if not hide_orig:
            surf = ax.plot_trisurf(x, y, z, color='r', alpha=.4, linewidth=0)
        surf_r = ax.plot_trisurf(xr, yr, zr, color='b', alpha=.4, linewidth=0)

    # plot lines along b3 vector direction passing through
    # b1 and b2 indices
    for k in range(0, 4):
        x = np.array((coords[k, 0], coords[k+8, 0]))
        xr = np.array((rot_coords[k, 0], rot_coords[k+8, 0]))
        y = np.array((coords[k, 1], coords[k + 8, 1]))
        yr = np.array((rot_coords[k, 1], rot_coords[k + 8, 1]))
        z = np.array((coords[k, 2], coords[k + 8, 2]))
        zr = np.array((rot_coords[k, 2], rot_coords[k + 8, 2]))
        if k==0:
            if not hide_orig:
                ax.plot(x, y, z, color='r', label='0-rotation lattice')
            ax.plot(xr, yr, zr, color='b', label='Rotated lattice')
        else:
            if not hide_orig:
                ax.plot(x, y, z, color='r')
            ax.plot(xr, yr, zr, color='b')

    # plot lattice index points
    i = 0
    for vi in v:
        randcolor = '#%06X' % randint(0, 0xFFFFFF)
        if not hide_orig:
            ax.scatter(vi[0], vi[1], vi[2], marker = 'o', color=randcolor)
        vr = vrot[i]
        ax.scatter(vr[0], vr[1], vr[2], marker='v', color=randcolor)
        i += 1

    # plot beamline (now at origin instead of center of unrotated lattice)
    beam_x = np.array([0,b1[0]*2])
    beam_y = np.array([0,0])
    beam_z = np.array([0, 0])
    ax.plot(beam_x, beam_y, beam_z, color='g', linestyle='dashed', label='Beam line')

    # formatting
    set_axes_equal(ax)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    if show:
        plt.show()
    else:
        return fig, ax

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])