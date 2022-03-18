from pylab import *
from KappaDiffractometer import KappaDiffractometer
from DiffracSample import DiffracSample
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from Arrow3D import *
from patch3d import *

def plot_scattering_geometry(sample, hkl, phi, kappa, omega):
    # Plot the scattering geometry
    # this uses the package `arrow3` from matlabcentral to plot:
    # https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3/

    # sample_inplane = [0,1,0];
    lambda0=sample.lambda0
    bg = sample.bg.transpose()
    Ghkl = sample.Ghkl
    k0 = sample.k0 # vector

    # parameters for plotting arrows
    arrow_width = 3
    arrow_height= 5
    o = np.array([0,0,0])
    X = np.array([1,0,0])
    Y = np.array([0,1,0])
    Z = np.array([0,0,1])

    # get the rotation matrix for the diffraction condition
    # Rtot = find_rotation_grazing_exit(geometry, sample_normal, hkl, beta);
    # rather use the Huber angles in `geometry` to define the scattering geometry:
    diff = KappaDiffractometer()
    Rtot = np.dot(diff.kappa_rotation_matrix(phi, kappa, omega),sample.sample_rmat);

    ### the final rotated vectors are:
    finalRotGhkl    = np.dot(Rtot,Ghkl)
    finalRotNsample = np.dot(Rtot, sample.normal_vect)
    finalkp         = k0 + finalRotGhkl


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot axes (?)
    #ax.arrow3D(o, X, 'r:', arrow_width/3, arrow_height)
    #ax.arrow3D(o, Y, 'g:', arrow_width/3, arrow_height)
    #ax.arrow3D(o, Z, 'b:', arrow_width/3, arrow_height)

    # plot `k0`
    ax.arrow3D(o, k0, mutation_scale=20, arrowstyle='-|>')
    # plot `rotGhkl`
    ax.arrow3D(o, finalRotGhkl, mutation_scale=20, arrowstyle='-|>')
    # plot `kp`
    ax.arrow3D(o, finalkp,  mutation_scale=20, arrowstyle='-|>')

    # plot the rotated normal `rotNsample`
    ax.arrow3D(o, 0.5*finalRotNsample, mutation_scale=20, arrowstyle='-|>')

    # Make a square starting from `sample_normal` and `sample_inplane`
    v2 = cross(sample.normal_vect, sample.inplane_vect)
    vs = np.stack([o, sample.inplane_vect, v2, v2+sample.inplane_vect], axis=-1)
    vs = vs-mean(vs,0) + o # center at `o`
    faces = [1,2,4,3]

    p1 = patch3d('vertices',vs, 'faces',faces);
    p1.FaceAlpha=0.8;
    p1.FaceColor=[.95,.5,.4];

    # find the final rotation axis and angle
    [V, e] = eig(Rtot);
    e = diag(e);
    # rot_axis = V(:,imag(e) == 0);
    #
    # # the trace of Rtot is 1 + 2*cos(th) with th the rotation angle
    # rot_angle= acosd((trace(Rtot) - 1)/2);
    #
    # # to fix the sign of the angle (not well-defined in the trace above) do
    # # this:
    # if norm(rotationmat3D(rot_angle, rot_axis) - Rtot) > 1e-15
    #     rot_angle = -rot_angle;
    # end
    #
    # rotate(p1, rot_axis, rot_angle, o);
    #
    # xlabel('x')
    # ylabel('y')
    # zlabel('z')



