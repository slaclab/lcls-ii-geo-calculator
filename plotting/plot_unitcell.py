import matplotlib.pyplot as plt
import numpy as np

def plot_lattice_unitcell(bg, t0=None):

    # plot the unit cell from a set of lattice vectors
    # matlab ref code calls the variable bg (inv space) but seems to also use same code for realspace (at) vectors?
    dim = np.shape(bg)[0]
    if t0 is None:
        t0 = np.zeros((dim,1))

    ibg = np.linalg.inv(bg).transpose()

    for i in range(0, dim):
        ibg[:, i] = ibg[:, i] / np.linalg.norm(ibg[:, i])

    v1 = -np.sum(np.multiply(ibg, bg[..., :]+t0), 0).transpose()
    v2 = np.sum(np.multiply(ibg[..., :], t0), 0).transpose()

    b = np.concatenate((v1, v2), 0)
    A = np.concatenate((-ibg.transpose(), ibg.transpose()), 0)

        




