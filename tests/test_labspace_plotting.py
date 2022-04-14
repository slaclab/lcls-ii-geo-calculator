from DiffracSample import DiffracSample
from plotting.plot_momentum_vectors import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from plotting.plot_lattice import set_axes_equal

# test plot of momentum transfer vectors
# and detector position representation

# spoof inputs
sample = DiffracSample()
sample.k0 = sample.beam_axis
sample.kp = np.array([1, 1, 0.001])
sample.kp = sample.kp / np.linalg.norm(sample.kp)
nu = np.arctan2(sample.kp[1], sample.kp[0])  # radians
delta = np.arcsin(sample.kp[2])
sample.ucell_a, sample.ucell_b, sample.ucell_c = [1,1,1]
center = np.array([0,0,0])

fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)

ax = plot_momentum_transfer(sample, center, ax)
ax = plot_detector_target(sample, center, nu, delta, ax)
set_axes_equal(ax)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()

