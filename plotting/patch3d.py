import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def patch3d(ax, x, y, z, v, vmin=0, vmax=100, cmap_name='viridis'):
    cmap = mpl.cm.get_cmap(cmap_name)               # Get colormap by name
    c = cmap(mpl.colors.Normalize(vmin, vmax)(v))   # Normalize value and get color
    pc = Poly3DCollection([list(zip(x,y,z))])       # Create PolyCollection from coords
    pc.set_facecolor(c)                             # Set facecolor to mapped value
    pc.set_edgecolor('k')                           # Set edgecolor to black
    ax.add_collection3d(pc)                         # Add PolyCollection to axes
    return pc

def view(ax, code):
    if code == 2: #view(2) sets the default two-dimensional view, az = 0, el = 90.
        ax.view_init(90, 0)     # (args are reversed from MATLAB)

    if code == 3: #view(3) sets the default three-dimensional view, az = –37.5, el = 30.
        ax.view_init(30, -37.5)

