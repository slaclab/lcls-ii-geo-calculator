import numpy as np
import matplotlib.pyplot as plt
from random import randint
from tools.angle_calc_utils import bravais_matrix
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

# general construction notes:
# gamma is angle between a, b
# alpha is angle between b, c
# beta is angle between a, c

UNIT_CUBE_VERTICES = np.array([[0, 0, 0],
                               [1, 0, 0],
                               [0, 1, 0],
                               [0, 0, 1],
                               [1, 1, 0],
                               [0, 1, 1],
                               [1, 0, 1],
                               [1, 1, 1]])
UNIT_CUBE_FACES = [[0, 1, 4, 2, 0],
                 [0, 1, 6, 3, 0],
                 [1, 4, 7, 6, 1],
                 [0, 2, 5, 3, 0],
                 [2, 4, 7, 5, 2],
                 [3, 5, 7, 6, 3]]

def plot_lattice_rotation(at, rotationmat3d, prev_rotationmat3d = None, show=True, hide_orig=False, ax=None):

    # for stacking or comparing initial/previous and final positions
    if prev_rotationmat3d is not None:
        rmat = np.dot(rotationmat3d, prev_rotationmat3d)
    else:
        rmat = rotationmat3d


    # take linear combinations of bravais vectors
    b1, b2, b3 = at[0, :], at[1, :], at[2, :] # using rows todo check
    o = np.zeros(3)
    # get vertices by combining vectors
    # figure consists of two stacked unit cells
    v = [o, b1, b2, b1+b2, b3, b1+b3, b2+b3, b1+b2+b3, b3*2, b1+b3*2, b2+b3*2, b1+b2+b3*2]
    vrot = [vi.dot(rmat) for vi in v]
    if prev_rotationmat3d is not None:
        v = [vi.dot(prev_rotationmat3d) for vi in v]
    coords = np.array(v)
    rot_coords = np.array(vrot)

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
    else:
        fig = None

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

def get_unitcell_hull(cell_lengths, cell_angles):
    # return as numpy array shape (8, 3) the vertices of the convex hull of
    # the lattice structure's unit cell (ignores internal points)
    [a, b, c] = cell_lengths
    [alpha, beta, gamma] = cell_angles

    # create numpy array of lattice 'point cloud' for a single cell
    bravais = bravais_matrix(a, b, c, alpha, beta, gamma)

    # this order of the dot product keeps the xy plane flat in 3d space
    unit_cell_vectors = np.dot(UNIT_CUBE_VERTICES, bravais)
    return unit_cell_vectors

def get_unitcell_hull_centered(cell_lengths, cell_angles):
    # return unit cell hull with center mass on origin
    unit_cell = get_unitcell_hull(cell_lengths, cell_angles)
    center = get_unitcell_center(unit_cell)

    return unit_cell - center


def plot_unitcell(unit_cell_vertices, rotation_mat = np.eye(3), ax=None, show=False):
    # plot face planes and outlines of unit cell external hull
    # optional rotation inputs to show in rotated position
    if ax is None:
        fig = plt.figure()
        ax = Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(ax)
    else:
        fig = None

    # vertices is shape (8,3)
    unit_cell_vertices = np.transpose(np.dot(rotation_mat, np.transpose(unit_cell_vertices)))
    x = unit_cell_vertices[:,0]
    y = unit_cell_vertices[:,1]
    z = unit_cell_vertices[:,2]

    unit_cell_vertices = unit_cell_vertices.tolist()
    p3d = [[unit_cell_vertices[UNIT_CUBE_FACES[ix][iy]]
            for iy in range(len(UNIT_CUBE_FACES[0]))]
            for ix in range(len(UNIT_CUBE_FACES))]
    ax.scatter(x, y, z)
    face_collection = Poly3DCollection(p3d, linewidths=1, alpha=0.2)
    line_collection = Line3DCollection(p3d, colors='k', linewidths=0.5)
    ax.add_collection3d(face_collection)
    ax.add_collection3d(line_collection)

    # formatting
    set_axes_equal(ax)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    if show==True:
        plt.show()
    return fig, ax


def get_unitcell_center(unit_cell_vertices):
    # return xyz coordinates of the center of the unit cell

    x = unit_cell_vertices[:,0]
    y = unit_cell_vertices[:,1]
    z = unit_cell_vertices[:,2]

    center = np.array([0.5*(np.amax(v) - np.amin(v)) for v in [x, y, z]])

    return center


    


