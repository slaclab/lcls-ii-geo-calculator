import numpy as np
import matplotlib.pyplot as plt
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

    side_labels = ["<-0->", "  side 'a'", "  side 'b'", "  side 'c'"]
    for i in range(4):
        ax.text(x[i], y[i], z[i], side_labels[i])

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


