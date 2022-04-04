import _tkinter
from tkinter import *
from tkinter import ttk, filedialog, scrolledtext
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg)
from tools.conversions import *
from KappaDiffractometer import KappaDiffractometer
from DiffracSample import DiffracSample
from plotting.plot_lattice import *
from plotting.plot_momentum_vectors import *


class BeamGeoApp(Tk):
    def __init__(self):
        super().__init__()

        # instantiate sample object, diffractometer object for GUI to modulate
        self.sample = DiffracSample()
        self.diffractometer = KappaDiffractometer()

        self.title("LCLSII Beam Geometry Calculator")

        self.resizable(1, 1)
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=2)

        # some GUI construction components need to be referenced outside of constructor functions
        # for changing associated text etc
        self.pathentry = None
        self.cell_a_entry = None
        self.cell_b_entry = None
        self.cell_c_entry = None
        self.cell_alpha_entry = None
        self.cell_beta_entry = None
        self.cell_gamma_entry = None
        self.fig1_mainaxes = None
        self.canvas1 = None

        # input variables to track
        self.h_miller = DoubleVar()
        self.l_miller = DoubleVar()
        self.k_miller = DoubleVar()
        self.snorm_x = DoubleVar()
        self.snorm_y = DoubleVar()
        self.snorm_z = DoubleVar()
        self.e_beam = DoubleVar()

        self.wavelen = DoubleVar()
        self.wavelen.set(0.0)
        self.incd_vs_exit = StringVar()
        self.angle_spec = DoubleVar()
        self.deltmax = DoubleVar()

        # lattice structure/inputs variables
        self.lattice_input_mode = StringVar()
        self.lattice_config_filepath = StringVar()
        self.cell_a = DoubleVar()
        self.cell_b = DoubleVar()
        self.cell_c = DoubleVar()
        self.cell_alpha = DoubleVar()
        self.cell_beta = DoubleVar()
        self.cell_gamma = DoubleVar()
        self.crystalgroup_list = ['<not set>', 'Cubic', 'Hexagonalic', 'Rhombic', 'Tetragonal',
                                  'Orthorhombic', 'Monoclinic', 'Triclinic']
        self.crystalgroup = StringVar()
        self.crystalgroup.set(self.crystalgroup_list[0])

        # output calculations (initialize to display zero)
        self.phi = DoubleVar()
        self.phi.set(0)
        self.kappa = DoubleVar()
        self.kappa.set(0)
        self.omega = DoubleVar()
        self.omega.set(0)
        self.delta = DoubleVar()
        self.delta.set(0)
        self.nu = DoubleVar()
        self.nu.set(0)
        self.k0 = StringVar()
        self.k0.set("[0,0,0]")
        self.kp = StringVar()
        self.kp.set("[0,0,0]")

        # track current 3d rotation for plotting (initialize at identity rotation)
        self.current_rotation = np.eye(3)
        self.last_rotation = None

        # call constructor functions
        leftcol = ttk.Frame(self)
        run_button = ttk.Button(leftcol, width=20, text="Run", command=lambda: self.run_command())
        rightcol = ttk.Frame(self)

        self.draw_latticeinfo_block(leftcol, 0, 0)
        self.draw_userinput_block(leftcol, 0, 1)

        run_button.grid(column=0, row=2, pady=5)
        self.draw_output_text(leftcol, 0, 3)

        self.draw_status_box(leftcol, 0, 4)

        self.draw_lattice_fig(rightcol, 0, 0)

        leftcol.grid(column=0, row=0)
        rightcol.grid(column=1, row=0)

    # overhead function to draw user inputs: calls functions to draw hkl, sample normal, beam energy etc
    def draw_userinput_block(self, parent, col, row):
        user_in = ttk.Frame(parent)
        user_in['borderwidth'] = 2
        user_in['relief'] = 'sunken'
        inputs_label = ttk.Label(user_in, text="User Inputs")
        inputs_label.grid(column=0, row=0)

        self.draw_hkl_input(user_in, 0, 1)
        self.draw_samplenorm_input(user_in, 0, 2)
        self.draw_beamenergy_input(user_in, 0, 3)
        self.draw_inc_exit_radios(user_in, 0, 4)
        self.draw_inc_ext_angle_entry(user_in, 0, 5)
        self.draw_delta_limit_box(user_in, 0, 6)

        user_in.grid(column=col, row=row, sticky=NW, padx=3, pady=3)

    # draw calculation outputs (text)
    def draw_output_text(self, parent, col, row):
        output_frame = ttk.Frame(parent)

        calcs_out = ttk.Frame(output_frame)
        calcs_out['borderwidth'] = 2
        calcs_out['relief'] = 'sunken'
        diff_angles_label = ttk.Label(calcs_out, text="Diffractometer Angles [deg]")
        diff_angles_label.grid(column=0, row=0, pady=3, padx=5)
        det_angles_label = ttk.Label(calcs_out, text="Detector Angles [deg]")
        det_angles_label.grid(column=0, row=2, pady=3, padx=5)

        # construct diffractometer angle outputs
        diffang_frame = ttk.Frame(calcs_out)
        diffang_frame['borderwidth'] = 1
        diffang_frame['relief'] = 'solid'
        klabel = ttk.Label(diffang_frame, text="Kappa:")
        k = ttk.Label(diffang_frame, textvariable=self.kappa)
        phlabel = ttk.Label(diffang_frame, text="   Phi:")
        ph = ttk.Label(diffang_frame, textvariable=self.phi)
        wlabel = ttk.Label(diffang_frame, text="   Omega:")
        w = ttk.Label(diffang_frame, textvariable=self.omega)
        i = 0
        for item in [klabel, k, phlabel, ph, wlabel, w]:
            item.grid(column=i, row=0)
            i += 1
        diffang_frame.grid(column=0, row=1)

        # construct detector angle outputs
        detang_frame = ttk.Frame(calcs_out)
        detang_frame['borderwidth'] = 1
        detang_frame['relief'] = 'solid'
        dlabel = ttk.Label(detang_frame, text="Delta:")
        d = ttk.Label(detang_frame, textvariable=self.delta)
        nlabel = ttk.Label(detang_frame, text="   Nu:")
        n = ttk.Label(detang_frame, textvariable=self.nu)
        i = 0
        for item in [dlabel, d, nlabel, n]:
            item.grid(column=i, row=0)
            i += 1
        detang_frame.grid(column=0, row=3)

        calcs_out.grid(column=0, row=0, sticky=NW, padx=3, pady=1)

        # additional outputs (momentum vector, bragg angle?)
        extras_frame = ttk.Frame(output_frame)
        extras_frame['borderwidth'] = 1
        extras_frame['relief'] = 'sunken'
        mvec_label = ttk.Label(extras_frame, text="Momentum transfer")
        kframe = ttk.Frame(extras_frame)
        k0_label = ttk.Label(kframe, text="k_0:")
        kp_label = ttk.Label(kframe, text="  k_p:")
        k0 = ttk.Label(kframe, textvariable=self.k0)
        kp = ttk.Label(kframe, textvariable=self.kp)
        k0_label.grid(column=0, row=0, padx=3)
        k0.grid(column=1, row=0)
        kp_label.grid(column=2, row=0, padx=3)
        kp.grid(column=3, row=0)

        mvec_label.grid(column=0, row=0)
        kframe.grid(column=0, row=1)

        extras_frame.grid(column=1, row=0, sticky=NW, padx=3, pady=1)

        output_frame.grid(column=col, row=row, sticky=NW)

    # draw a text box to print status and error messages
    def draw_status_box(self, parent, col, row):
        frame = ttk.Frame(parent)
        textbox = scrolledtext.ScrolledText(frame, width=60, height=6)
        textbox.grid(column=0, row=0, sticky=W)
        frame.grid(column=col, row=row, sticky=NW)

    # overhead function calls commands to draw lattice radios, file inputs, manual entries
    def draw_latticeinfo_block(self, parent, col, row):
        lattice_frame = ttk.Frame(parent)
        lattice_frame['borderwidth'] = 2
        lattice_frame['relief'] = 'sunken'
        lattice_label = ttk.Label(lattice_frame, text="Lattice Structure Inputs")
        lattice_label.grid(column=0, row=0)

        self.draw_lattice_filepath_in(lattice_frame, 0, 1)
        self.draw_man_lattice_inputs(lattice_frame, 0, 2)

        lattice_frame.grid(column=col, row=row, sticky=NW, padx=3, pady=3)

    # draw figure placeholder
    def draw_lattice_fig(self, parent, col, row):
        latfig_frame = ttk.Frame(parent)  # put in frame just in case? for formatting
        fig = plt.figure()
        self.fig1_mainaxes = fig.add_subplot(projection='3d')
        self.canvas1 = FigureCanvasTkAgg(fig, master=latfig_frame)
        self.canvas1.get_tk_widget().grid(column=0, row=0)
        self.canvas1.draw()

        latfig_frame.grid(column=col, row=row)

    # draw input boxes for miller indices HKL
    def draw_hkl_input(self, parent, col, row):
        hkl = ttk.Frame(parent)
        hkl_label = ttk.Label(hkl, text='Miller Indices')
        hlabel = ttk.Label(hkl, text='H:')
        klabel = ttk.Label(hkl, text='K:')
        llabel = ttk.Label(hkl, text='L:')
        h_in = ttk.Entry(hkl, width=4, textvariable=self.h_miller)
        k_in = ttk.Entry(hkl, width=4, textvariable=self.k_miller)
        l_in = ttk.Entry(hkl, width=4, textvariable=self.l_miller)
        hkl_label.grid(column=0, row=0, sticky=W, padx=5)

        hlabel.grid(column=1, row=0)
        h_in.grid(column=2, row=0, padx=2)

        klabel.grid(column=3, row=0)
        k_in.grid(column=4, row=0, padx=2)

        llabel.grid(column=5, row=0)
        l_in.grid(column=6, row=0, padx=2)

        hkl.grid(column=col, row=row, padx=2, pady=2, sticky=W)

    # draw input boxes for sample normal vector
    def draw_samplenorm_input(self, parent, col, row):
        sn_frame = ttk.Frame(parent)
        vec_elements = ttk.Frame(sn_frame)
        sn_label1 = ttk.Label(sn_frame, text='Specify vector normal to sample face plane')
        sn_label2 = ttk.Label(sn_frame, text='using same (sample) basis as unit cell vectors')
        xlabel = ttk.Label(vec_elements, text='x:')
        ylabel = ttk.Label(vec_elements, text='y:')
        zlabel = ttk.Label(vec_elements, text='z:')
        x_in = ttk.Entry(vec_elements, width=4, textvariable=self.snorm_x)
        y_in = ttk.Entry(vec_elements, width=4, textvariable=self.snorm_y)
        z_in = ttk.Entry(vec_elements, width=4, textvariable=self.snorm_z)

        sn_label1.grid(column=0, row=0, sticky=W, padx=5)
        sn_label2.grid(column=0, row=1, sticky=W, padx=5)
        vec_elements.grid(column=0, row=2, sticky=W, padx=5)

        xlabel.grid(column=0, row=0)
        x_in.grid(column=1, row=0, padx=2)

        ylabel.grid(column=2, row=0)
        y_in.grid(column=3, row=0, padx=2)

        zlabel.grid(column=4, row=0)
        z_in.grid(column=5, row=0, padx=2)

        sn_frame.grid(column=col, row=row, padx=2, pady=5, sticky=W)

    # draw beam energy input box
    def draw_beamenergy_input(self, parent, col, row):
        # beam energy input
        e_beam_frame = ttk.Frame(parent)
        elabel = ttk.Label(e_beam_frame, text="Photon Energy [eV]: ")
        lamlabel = ttk.Label(e_beam_frame, text=" -> wavelength [Ang]: ")
        lamdisp = ttk.Label(e_beam_frame, textvariable=self.wavelen)
        e_in = ttk.Entry(e_beam_frame, width=6, textvariable=self.e_beam,
                         validate="focusout", validatecommand=lambda: self.update_wavelen())
        elabel.grid(column=0, row=0, padx=2)
        e_in.grid(column=1, row=0)
        lamlabel.grid(column=2, row=0)
        lamdisp.grid(column=3, row=0)

        e_beam_frame.grid(column=col, row=row, padx=5, pady=5)

    # DEPRECATED: draw radio button selection for lattice input mode
    def draw_latticein_radios(self, parent, col, row):
        # select input style for lattice info
        radioframe = ttk.Frame(parent)
        select_cif = ttk.Radiobutton(radioframe, text='Read CIF/.vasp file', variable=self.lattice_input_mode,
                                     value='file')
        # select_vasp = ttk.Radiobutton(radioframe, text='.vasp file', variable=self.lattice_input_mode, value='vasp')
        select_man = ttk.Radiobutton(radioframe, text='Manual input', variable=self.lattice_input_mode, value='man')

        select_cif.grid(column=0, row=0, padx=5)
        # select_vasp.grid(column=1, row=0, padx=5)
        select_man.grid(column=1, row=0, padx=5)

        radioframe.grid(column=col, row=row, padx=5, pady=2, sticky=W)

    # draw filepath selection tools
    def draw_lattice_filepath_in(self, parent, col, row):
        # input filepath for lattice input file if applicable
        latticepath_frame = ttk.Frame(parent)
        pathlabel = ttk.Label(latticepath_frame, text="Filepath")
        self.pathentry = ttk.Entry(latticepath_frame, width=30, textvariable=self.lattice_config_filepath)
        browser_button = ttk.Button(latticepath_frame, width=5, text="Browse", command=lambda: self.open_file_browser())
        readfile_button = ttk.Button(latticepath_frame, width=15, text="Read File")
        pathlabel.grid(column=0, row=0, sticky=W)
        self.pathentry.grid(column=0, row=1, sticky=W)
        browser_button.grid(column=1, row=1, sticky=W)

        latticepath_frame.grid(column=col, row=row, padx=3, pady=3, sticky=W)

    # draw manual input boxes for lattice inputs + space group dropdown menu
    def draw_man_lattice_inputs(self, parent, col, row):
        # include inputs for lattice vectors, angles etc
        # autofill with results from CIF file or else let user specify
        # also, dropdown with cell crystal type

        spec_lattice_frame = ttk.Frame(parent)
        cellvec_text1 = ttk.Label(spec_lattice_frame, text="Unit cell lengths [A]")

        inputs_frame = ttk.Frame(spec_lattice_frame)
        alabel = ttk.Label(inputs_frame, text="a:")
        blabel = ttk.Label(inputs_frame, text="b:")
        clabel = ttk.Label(inputs_frame, text="c:")
        self.cell_a_entry = ttk.Entry(inputs_frame, width=10, textvariable=self.cell_a,
                                      validate="focusout",
                                      validatecommand=lambda: self.update_cell_a())
        self.cell_b_entry = ttk.Entry(inputs_frame, width=10, textvariable=self.cell_b,
                                      validate="focusout",
                                      validatecommand=lambda: self.update_cell_b())
        self.cell_c_entry = ttk.Entry(inputs_frame, width=10, textvariable=self.cell_c,
                                      validate="focusout",
                                      validatecommand=lambda: self.update_cell_c())
        alabel.grid(column=0, row=0)
        self.cell_a_entry.grid(column=1, row=0)
        blabel.grid(column=0, row=1)
        self.cell_b_entry.grid(column=1, row=1)
        clabel.grid(column=0, row=2)
        self.cell_c_entry.grid(column=1, row=2)

        lattang_title = ttk.Label(spec_lattice_frame, text="Unit cell angles [deg]")
        angles_frame = ttk.Frame(spec_lattice_frame)
        alphlabel = ttk.Label(angles_frame, text="alpha:")
        betalabel = ttk.Label(angles_frame, text="beta:")
        gammlabel = ttk.Label(angles_frame, text="gamma:")
        self.cell_alpha_entry = ttk.Entry(angles_frame, width=10, textvariable=self.cell_alpha,
                                          validate="focusout",
                                          validatecommand=lambda: self.update_cell_alpha())
        self.cell_beta_entry = ttk.Entry(angles_frame, width=10, textvariable=self.cell_beta,
                                         validate="focusout",
                                         validatecommand=lambda: self.update_cell_beta())
        self.cell_gamma_entry = ttk.Entry(angles_frame, width=10, textvariable=self.cell_gamma,
                                          validate="focusout",
                                          validatecommand=lambda: self.update_cell_gamma())
        alphlabel.grid(column=0, row=0)
        self.cell_alpha_entry.grid(column=1, row=0)
        betalabel.grid(column=0, row=1)
        self.cell_beta_entry.grid(column=1, row=1)
        gammlabel.grid(column=0, row=2)
        self.cell_gamma_entry.grid(column=1, row=2)

        crystalgroup_label = ttk.Label(spec_lattice_frame, text="Crystal space group:")
        crystalgroup_ddown = ttk.OptionMenu(spec_lattice_frame, self.crystalgroup, *self.crystalgroup_list)
        crystalgroup_ddown.config(width=15)

        cellvec_text1.grid(column=0, row=0)
        lattang_title.grid(column=1, row=0)

        inputs_frame.grid(column=0, row=2, padx=5)
        angles_frame.grid(column=1, row=2)

        crystalgroup_label.grid(column=0, row=3, sticky=W)
        crystalgroup_ddown.grid(column=1, row=3, sticky=W)

        spec_lattice_frame.grid(column=col, row=row, sticky=W)

    # draw radiobutton selection for angle constraint type
    def draw_inc_exit_radios(self, parent, col, row):
        # angle selection radio buttons (incidence vs exit angle specification)
        angle_selection = ttk.Frame(parent)
        select_incidence = ttk.Radiobutton(angle_selection, text='Set Incidence Angle', variable=self.incd_vs_exit,
                                           value='inc')
        select_exit = ttk.Radiobutton(angle_selection, text='Set Exit Angle', variable=self.incd_vs_exit,
                                      value='ext')
        select_free = ttk.Radiobutton(angle_selection, text='Unconstrained', variable=self.incd_vs_exit,
                                      value='free')
        select_incidence.grid(column=0, row=0, padx=5)
        select_exit.grid(column=1, row=0, padx=5)
        select_free.grid(column=2, row=0, padx=5, sticky=W)

        angle_selection.grid(column=col, row=row, padx=5, sticky=W)

    # draw input box for incidence or exit angle
    def draw_inc_ext_angle_entry(self, parent, col, row):
        # entry boxes to specify entry/exit angle
        anglespec_frame = ttk.Frame(parent)

        alabel = ttk.Label(anglespec_frame, text="Specify Incidence or Exit Angle:")
        unitslabel = ttk.Label(anglespec_frame, text='[deg]')
        angle_val = ttk.Entry(anglespec_frame, width=4, textvariable=self.angle_spec)
        alabel.grid(column=0, row=0, padx=2)
        angle_val.grid(column=1, row=0)
        unitslabel.grid(column=2, row=0, padx=2)

        anglespec_frame.grid(column=col, row=row, padx=5, sticky=W)

    # draw input box for detector angle limits
    # TODO add limits on kappa angle
    def draw_delta_limit_box(self, parent, col, row):
        # entry boxes to specify limits on allowable detector angles
        detector_lim_frame = ttk.Frame(parent)

        deltmax_label = ttk.Label(detector_lim_frame, text="Set max detector angle Delta (azimuthal):")
        unitslabel = ttk.Label(detector_lim_frame, text='[deg]')
        deltmax = ttk.Entry(detector_lim_frame, width=4, textvariable=self.deltmax)

        deltmax_label.grid(column=0, row=0, sticky=W)
        deltmax.grid(column=1, row=0, sticky=W)
        unitslabel.grid(column=2, row=0, sticky=W, padx=2)

        detector_lim_frame.grid(column=col, row=row, sticky=W, padx=5, pady=2)

    # use system file browser to locate cif/vasp filepaths
    def open_file_browser(self):
        file = filedialog.askopenfilename()
        if file is not None:
            self.lattice_config_filepath = file

            # place the selected file into the path entry box for reference
            if self.pathentry is not None:
                self.pathentry.delete(0, END)
                self.pathentry.insert(0, file)
            try:
                self.read_lattice_file()
            except:
                pass
                # TODO exception handling

            # TODO other checks on file etc
        else:
            pass
            # TODO print error outputs

    # if unit cell information provided not-manually, update the text boxes to reflect new info
    def update_unitcell_readout(self):
        # angles associated with the GUI class and its internal var's will be in degrees
        # angles associated with diffractometer or sample classes will be in radians
        lengths = [self.sample.ucell_a, self.sample.ucell_b, self.sample.ucell_c]
        angles = [np.round(np.rad2deg(k), 3) for k in
                  [self.sample.ucell_alpha, self.sample.ucell_beta,
                   self.sample.ucell_gamma]]  # these are internal storage so units radians

        update_vals = lengths + angles
        disp_update_list = [self.cell_a_entry, self.cell_b_entry, self.cell_c_entry,
                            self.cell_alpha_entry, self.cell_beta_entry, self.cell_gamma_entry]
        var_update_list = [self.cell_a, self.cell_b, self.cell_c,
                           self.cell_alpha, self.cell_beta, self.cell_gamma]

        for i in range(6):
            var_update_list[i].set(update_vals[i])
            disp_update_list[i].delete(0, END)
            disp_update_list[i].insert(0, update_vals[i])

    # try reading input file to populate lattice structure variables
    def read_lattice_file(self):
        path = self.lattice_config_filepath
        if path is not None:
            if path.endswith(".cif"):
                print("Found cif file")
                self.sample.readCIF(path)
                self.update_unitcell_readout()
            elif path.endswith(".vasp"):
                print("vasp file reader not yet implemented")
                # TODO vasp file reader
            else:
                pass
        else:
            pass
        self.update_figure()

    ## the following 6 functions are for individual lattice parameter updates
    # TODO wrap things in try-except logic probably
    def update_cell_a(self):
        a = self.cell_a.get()
        self.sample.ucell_a = a
        self.sample.try_update_lattice()
        self.update_figure()
        return False

    def update_cell_b(self):
        b = self.cell_b.get()
        self.sample.ucell_b = b
        self.sample.try_update_lattice()
        self.update_figure()
        return False

    def update_cell_c(self):
        c = self.cell_c.get()
        self.sample.ucell_c = c
        self.sample.try_update_lattice()
        self.update_figure()
        return False

    def update_cell_alpha(self):
        # convert to radians before setting value in sample object
        alpha = np.deg2rad(self.cell_alpha.get())
        self.sample.ucell_alpha = alpha
        self.sample.try_update_lattice()
        self.update_figure()
        return False

    def update_cell_beta(self):
        beta = np.deg2rad(self.cell_beta.get())
        self.sample.ucell_beta = beta
        self.sample.try_update_lattice()
        self.update_figure()
        return False

    def update_cell_gamma(self):
        gamma = np.deg2rad(self.cell_gamma.get())
        self.sample.ucell_gamma = gamma
        self.sample.try_update_lattice()
        self.update_figure()
        return False

    # update wavelength variable internally/on display when energy is changed
    def update_wavelen(self):
        try:
            energy_float = self.e_beam.get()
            wavelen = photonE_to_wavelen_A(energy_float)
            self.wavelen.set(np.round(wavelen, 5))
            self.sample.set_beam_energy(energy_float)
            self.sample.k0 = self.sample.get_k0()
        except:
            self.wavelen.set(0)
        return False

    # collect the HKL indices and send them to DiffracSample object
    def update_hkl_from_inputs(self):
        # replace empty strings with zeros
        for item in [self.h_miller, self.k_miller, self.l_miller]:
            try:
                item.get()
            except _tkinter.TclError:
                item.set(0.0)
        self.sample.bragg_hkl = np.array([self.h_miller.get(), self.k_miller.get(), self.l_miller.get()])

    def update_samplenormal_from_inputs(self):
        # replace empty strings with zeros
        for item in [self.snorm_x, self.snorm_y, self.snorm_z]:
            try:
                item.get()
            except _tkinter.TclError:
                item.set(0.0)

        self.sample.set_sample_normal(np.array([self.snorm_x.get(), self.snorm_y.get(), self.snorm_z.get()]))

    def update_figure(self):

        # if no lattice defined, no update
        if not self.sample.is_lattice_defined:
            return
        else:
            # update figure object
            # if no rotation, will plot unrotated lattice
            self.fig1_mainaxes.clear()
            unit_cell = get_unitcell_hull_centered([self.sample.ucell_a, self.sample.ucell_b, self.sample.ucell_c],
                                          [self.sample.ucell_alpha, self.sample.ucell_beta, self.sample.ucell_gamma])
            fig, self.fig1_mainaxes = plot_unitcell(unit_cell, self.current_rotation, ax=self.fig1_mainaxes)

            self.fig1_mainaxes = plot_momentum_transfer(self.sample, np.array([0,0,0]), self.fig1_mainaxes)

            self.canvas1.draw()


    def run_command(self):
        # try: run geometry calcs with provided inputs

        # collect hkl, sample incidence plane normal inputs and port to backend classes
        self.update_hkl_from_inputs()
        self.update_samplenormal_from_inputs()

        opti_init_condition = -15  # TODO this is arbitrary for now (and therefore maybe not optimal)

        if self.incd_vs_exit.get() == "inc":
            # calculate based on constrained incidence angle
            print("calculating incidence angle rotation")
            incd_angle = np.deg2rad(self.angle_spec.get())
            ref_rotation = self.sample.find_rotation_grazing_incidence(incd_angle, opti_init_condition)

        elif self.incd_vs_exit.get() == "ext":
            print("Calculating exit angle rotation")
            # calculate based on constrained exit angle
            exit_angle = np.deg2rad(self.angle_spec.get())
            ref_rotation = self.sample.find_rotation_grazing_exit(exit_angle, opti_init_condition)
            self.last_rotation = self.current_rotation
            self.current_rotation = ref_rotation
        else:
            # TODO implement unconstrained solution
            pass
        print("Reference rotation:")
        print(ref_rotation)

        # calculate instrument angles
        self.diffractometer.set_kappa_from_eulerian(self.current_rotation)
        self.diffractometer.calc_detector_angles(self.sample.kp, self.sample.lambda0)

        # update display
        angle_decimal_places = 4
        self.kappa.set(np.round(np.rad2deg(self.diffractometer.kappa), angle_decimal_places))
        self.phi.set(np.round(np.rad2deg(self.diffractometer.phi), angle_decimal_places))
        self.omega.set(np.round(np.rad2deg(self.diffractometer.omega), angle_decimal_places))
        self.nu.set(np.round(np.rad2deg(self.diffractometer.nu), angle_decimal_places))
        self.delta.set(np.round(np.rad2deg(self.diffractometer.delta), angle_decimal_places))

        # output momentum vector
        k0_string = np.array2string(self.sample.k0, formatter={'float_kind': lambda x: "%.2f" % x}, separator=',')
        kp_string = np.array2string(self.sample.kp, formatter={'float_kind': lambda x: "%.2f" % x}, separator=',')
        self.k0.set(k0_string)
        self.kp.set(kp_string)

        self.update_figure()


app = BeamGeoApp()
app.mainloop()
