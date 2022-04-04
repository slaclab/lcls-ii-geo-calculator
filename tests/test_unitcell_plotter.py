import numpy as np
from plotting.plot_lattice import plot_unitcell, get_unitcell_hull_centered
from tools.angle_calc_utils import rotationmat3D

cell_lengths = [2,3,4]
cell_angles = [np.deg2rad(70), np.deg2rad(70), np.deg2rad(90)]

vertices = get_unitcell_hull_centered(cell_lengths, cell_angles)
fig, ax = plot_unitcell(vertices, show=True)

rot_90_about_z = rotationmat3D(np.pi/2, [0, 0, 1])
figrot, axrot = plot_unitcell(vertices, rotation_mat=rot_90_about_z, show=True)

