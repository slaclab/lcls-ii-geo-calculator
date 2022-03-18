from pylab import *
from tools.angle_calc_utils import *
from KappaDiffractometer import KappaDiffractometer
import plotly.graph_objects as go
import plotly.express as px

def plotly_lattice(at, phi, kappa, omega, sample_rot=None,sample_rot_axis=np.array([0,0,1]), kappa_rot_mat=None):

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
    #print(coords)
    #print(np.shape(coords)) #12x3
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # plot surfaces parallel to b1, b2 vectors and perpendicular to b3
    for j in [0, 4, 8]:
        x, y, z = coords[j:j+4, 0], coords[j:j+4, 1], coords[j:j+4, 2]
        xr, yr, zr = rot_coords[j:j+4, 0], rot_coords[j:j+4, 1], rot_coords[j:j+4, 2]
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
        ax.plot(x, y, z, color='r')
        ax.plot(xr, yr, zr, color='b')


    # plot lattice index points
    i = 0
    for vi in v:
        randcolor = '#%06X' % randint(0, 0xFFFFFF)
        ax.scatter(vi[0], vi[1], vi[2], marker = 'o', color=randcolor)
        vr = vrot[i]
        ax.scatter(vr[0], vr[1], vr[2], marker='v', color=randcolor)
        i += 1

    # plot beamline
    beam_x = np.array([b1[0]/2, b1[0]/2])
    beam_y = np.array([b2[1]/2, b2[1]/2])
    beam_z = np.array([0, b3[2]*2])
    ax.plot(beam_x, beam_y, beam_z, color='g', linestyle='dashed')

    # formatting
    set_axes_equal(ax)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.show()