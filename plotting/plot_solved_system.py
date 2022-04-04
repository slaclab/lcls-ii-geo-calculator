from plotting.plot_lattice import *
from tools.angle_calc_utils import rotationmat3D

def plot_diffrac_system(sample, rotation_matrix, nu, delt, ax=None):

    xax = np.array([1,0,0])
    yax = np.array([0,1,0])
    zax = np.array([0,0,1])

    scale = np.max(sample.at)

    # get plot of realspace lattice
    fig, ax = plot_lattice_rotation(sample.at, rotation_matrix, prev_rotationmat3d=np.eye(3),
                                    show=False, hide_orig=True, ax=ax)

    # turn normal/inplane vectors into the sample face plane
    f1 = sample.inplane_vect.dot(rotation_matrix)
    f2 = -f1
    f3 = np.cross(sample.normal_vect.dot(rotation_matrix), f1)
    f4 = -f3
    face = np.stack([f1, f2, f3, f4], axis=1)* scale
    face_x = face[0,:]
    face_y = face[1,:]
    face_z = face[2,:]

    ax.plot_trisurf(face_x, face_y, face_z, color='k', alpha=.4, linewidth=0)
    ax.plot(0,0,0, color='k', label='Sample face plane')

    detector_zeroline = xax
    detector_rmat = rotationmat3D(delt, -yax).dot(rotationmat3D(nu, zax))
    detector_line = detector_zeroline.dot(detector_rmat)
    detx = np.array([-detector_line[0], detector_line[0]]) * scale
    dety = np.array([-detector_line[1], detector_line[1]])* scale
    detz = np.array([-detector_line[2], detector_line[2]])* scale
    ax.plot(detx, dety, detz, color='magenta', linestyle='dashed', label='detector')
    ax.legend()
    #plt.show()
    return fig, ax

