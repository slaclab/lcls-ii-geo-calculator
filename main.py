import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from KappaDiffractometer import KappaDiffractometer
from DiffracSample import DiffracSample
from cmd import Cmd
import os
from datetime import datetime
from plotting.plot_lattice import get_unitcell_hull_centered, plot_unitcell
from plotting.plot_momentum_vectors import plot_momentum_transfer, plot_detector_target


## ALL CALCULATIONS IN RADIANS,         ##
## DEGREES ONLY FOR DISPLAY/UI PURPOSES ##

# Provide a command line interface for setting up and running geometry calculations
# typing '?' in the running interface brings up a list of available commands
class GeoCalculatorShell(Cmd):
    prompt = 'geo>> '
    intro = "Welcome! Type ? to list commands"
    def __init__(self):
        super().__init__()
        self.sample = DiffracSample()
        self.diffractometer = KappaDiffractometer()
        self.current_rotation = np.eye(3)
        self.last_calc_type = None
        self.fixed_angle = None

    def do_exit(self, inp):
        '''Type exit to quit app.'''
        print("Quitting application...")
        return True

    def do_search_CIF_files(self, directory):
        '''Search a given directory for CIF filepaths and list in the terminal [NOT IMPLEMENTED]'''

    def do_readCIF(self, file_input):
        '''Read a CIF file by name (if in searchable location) or filepath.'''
        #TODO better filepath completion or figure out a browse functionality
        if not file_input.endswith('.cif'):
            print("Missing .cif extension or file is wrong type!")
        else:
            if os.path.exists(file_input):
                print("Trying to read CIF file...")
                self.sample.readCIF(file_input)
                print("Found data:")
                print("------------------------------")
                print("Unit cell dimensions")
                print("Cell lengths in Angstrom:")
                print("a: ", str(self.sample.ucell_a), " | ",
                      "b: ", str(self.sample.ucell_b), " | ",
                      "c: ", str(self.sample.ucell_c))
                print("Cell angles in degrees:")
                print("alpha: ", str(np.rad2deg(self.sample.ucell_alpha)), " | ",
                      "beta: ", str(np.rad2deg(self.sample.ucell_beta)), " | ",
                      "gamma: ", str(np.rad2deg(self.sample.ucell_gamma)))
                print("------------------------------")
            else:
                print("No item at <" + file_input + "> found!")

    def do_set_cell_length_a(self, a_in):
        '''edit unit cell side length a (input in Angstrom as a single numerical value)'''
        try:
            self.sample.ucell_a = float(a_in.strip())
            self.sample.try_update_lattice()
        except TypeError:
            print("Please provide a numerical input")

    def do_set_cell_length_b(self, b_in):
        '''edit unit cell side length b (input in Angstrom as a single numerical value)'''
        try:
            self.sample.ucell_b = float(b_in.strip())
            self.sample.try_update_lattice()
        except TypeError:
            print("Please provide a numerical input")

    def do_set_cell_length_c(self, c_in):
        '''edit unit cell side length c (input in Angstrom as a single numerical value)'''
        try:
            self.sample.ucell_c = float(c_in.strip())
            self.sample.try_update_lattice()
        except TypeError:
            print("Please provide a numerical input")

    def do_set_all_cell_lengths(self, input):
        '''Specify unit cell lengths as comma-separated triplet a,b,c or single value for cubic cells.'''
        try:
            input = input.strip().split(',')
            if np.size(input) > 1:
                [self.sample.ucell_a, self.sample.ucell_b, self.sample.ucell_c] = [float(value.strip()) for value in input]
            else:
                [self.sample.ucell_a, self.sample.ucell_b, self.sample.ucell_c] = [float(input[0])] * 3
            self.sample.try_update_lattice()
        except TypeError as e:
            print("Check input formatting!")
            print("sys error statements:")
            print(e)

    def do_set_cell_angle_alpha(self, alpha_in):
        '''edit unit cell angle alpha (input in degrees as single numerical value)'''
        try:
            alpha_deg = float(alpha_in.strip())
            self.sample.ucell_alpha = np.deg2rad(alpha_deg)
            self.sample.try_update_lattice()
        except TypeError:
            print("Please provide a numerical input")

    def do_set_cell_angle_beta(self, beta_in):
        '''edit unit cell angle beta (input in degrees as single numerical value)'''
        try:
            beta_deg = float(beta_in.strip())
            self.sample.ucell_beta = np.deg2rad(beta_deg)
            self.sample.try_update_lattice()
        except TypeError:
            print("Please provide a numerical input")

    def do_set_cell_angle_gamma(self, gamma_in):
        '''edit unit cell angle gamma (input in degrees as single numerical value)'''
        try:
            gamma_deg = float(gamma_in.strip())
            self.sample.ucell_gamma = np.deg2rad(gamma_deg)
            self.sample.try_update_lattice()
        except TypeError:
            print("Please provide a numerical input")

    def do_set_all_cell_angles(self, input):
        '''Specify unit cell angles as comma-separated triplet alpha,beta,gamma or single value for identical angles.'''
        try:
            input = input.strip().split(',')
            if np.size(input) > 1:
                [self.sample.ucell_alpha, self.sample.ucell_beta, self.sample.ucell_gamma] = [np.deg2rad(float(value.strip())) for value in input]
            else:
                [self.sample.ucell_alpha, self.sample.ucell_beta, self.sample.ucell_gamma] = [np.deg2rad(float(input[0]))] * 3
            self.sample.try_update_lattice()
        except TypeError as e:
            print("Check input formatting!")
            print("sys error statements:")
            print(e)

    def do_set_photon_energy(self, E_eV):
        '''Provide incident photon energy in units of eV.'''
        try:
            E_eV = float(E_eV.strip())
            self.sample.set_beam_energy(E_eV)
            self.sample.get_k0()
        except TypeError as e:
            print("Please provide a numerical input")

    def do_HKL_input(self, hkl_string):
        '''Provide Miller indices as comma-separated float values h,k,l.'''
        try:
            hkl_list = hkl_string.strip().split(',')
            [h, k, l] = [float(item.strip()) for item in hkl_list]
            self.sample.bragg_hkl = np.array([h, k, l])
        except ValueError or TypeError:
            print("Please check input formatting; values should be comma-separated floats with no other characters.")

    def do_select_incident_face(self, samplenorm_string):
        '''Provide normal vector for sample face incident to incoming radiation. Example input 1,0,0
            sets incident face to be parallel to unit cell "a" direction.'''
        try:
            samplenorm_list = samplenorm_string.strip().split(',')
            [x,y,z] = [float(item.strip()) for item in samplenorm_list]
            self.sample.set_sample_normal(np.array([x,y,z]))
        except ValueError or TypeError:
            print("Check formatting: values should be comma-separated numerical")

    def do_print_working_sample_info(self, inp):
        '''Display information about sample currently loaded into program.'''

        print("------------------------------")
        print("Unit cell dimensions")
        print("Cell lengths in Angstrom:")
        print("a: ", str(self.sample.ucell_a), " | ",
              "b: ", str(self.sample.ucell_b), " | ",
              "c: ", str(self.sample.ucell_c))
        print("Cell angles in degrees:")
        if not None in [self.sample.ucell_alpha, self.sample.ucell_beta, self.sample.ucell_gamma]:
            print("alpha: ", str(np.rad2deg(self.sample.ucell_alpha)), " | ",
                  "beta: ", str(np.rad2deg(self.sample.ucell_beta)), " | ",
                  "gamma: ", str(np.rad2deg(self.sample.ucell_gamma)))
        else:
            angles = [self.sample.ucell_alpha, self.sample.ucell_beta, self.sample.ucell_gamma]
            for i in range(3):
                if angles[i] is None:
                    angles[i] = 0
            angles_str = "alpha: {0} | beta: {1} | gamma: {2}".format(*[str(np.rad2deg(ang)) for ang in angles])
            print(angles_str)
        print("------------------------------")
        print("Photon energy: ", str(self.sample.beam_energy), " [eV]")
        print("Photon wavelength: ", str(self.sample.lambda0), " [Angstrom]")
        print("------------------------------")
        print("Miller indices")
        print("H: ", str(self.sample.bragg_hkl[0]), " | ",
              "K: ", str(self.sample.bragg_hkl[1]), " | ",
              "L: ", str(self.sample.bragg_hkl[2]))
        print("------------------------------")

    def do_calc_rotation_for_incidence_angle(self, inc_angle):
        '''Calculate sample rotation and instrument positions for input incidence angle in degrees.'''
        try:
            alpha = np.deg2rad(float(inc_angle.strip()))
        except TypeError or ValueError:
            print("Check that incidence angle input is float or int")
            return False # TODO make sure this doesn't exit app

        self.fixed_angle = inc_angle
        rotation_3d = self.sample.find_rotation_grazing_incidence(alpha, -15)
        self.diffractometer.set_kappa_from_eulerian(rotation_3d)
        self.diffractometer.calc_detector_angles(self.sample.kp, self.sample.lambda0)
        self.last_calc_type = "Fixed incidence angle"
        print("Rotation transformation matrix:")
        print(np.array2string(rotation_3d))
        print("Instrument angles:")
        print("Phi: " + str(np.rad2deg(self.diffractometer.phi)))
        print("Kappa: " + str(np.rad2deg(self.diffractometer.kappa)))
        print("Omega: " + str(np.rad2deg(self.diffractometer.omega)))
        print("Detector nu: " + str(np.rad2deg(self.diffractometer.nu)))
        print("Detector delta: " + str(np.rad2deg(self.diffractometer.delta)))
        print("Momentum Transfer")
        print("k0: " + np.array2string(self.sample.k0))
        print("kp: " + np.array2string(self.sample.kp))

    def do_calc_rotation_for_exit_angle(self, ext_angle):
        '''Calculate sample rotation and instrument positions for input exit angle in degrees.'''
        try:
            beta = np.deg2rad(float(ext_angle.strip()))
        except TypeError or ValueError:
            print("Check that incidence angle input is float or int")
            return False  # TODO make sure this doesn't exit app

        self.fixed_angle = ext_angle
        rotation_3d = self.sample.find_rotation_grazing_exit(beta, -15)
        self.diffractometer.set_kappa_from_eulerian(rotation_3d)
        self.diffractometer.calc_detector_angles(self.sample.kp, self.sample.lambda0)
        self.last_calc_type = "Fixed exit angle"
        print("Rotation transformation matrix:")
        print(np.array2string(rotation_3d))
        print("Instrument angles:")
        print("Phi: " + str(np.rad2deg(self.diffractometer.phi)))
        print("Kappa: " + str(np.rad2deg(self.diffractometer.kappa)))
        print("Omega: " + str(np.rad2deg(self.diffractometer.omega)))
        print("Detector nu: " + str(np.rad2deg(self.diffractometer.nu)))
        print("Detector delta: " + str(np.rad2deg(self.diffractometer.delta)))
        print("Momentum Transfer")
        print("k0: " + np.array2string(self.sample.k0))
        print("kp: " + np.array2string(self.sample.kp))

    def do_make_unitcell_plot(self, inp):
        '''Produce plot showing unit cell shape/orientation with momentum transfer if applicable'''
        fig = plt.figure()
        ax = Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(ax)

        if not self.sample.is_lattice_defined:
            print("Error: not enough lattice information specified to plot unit cell!")
            return
        else:
            unit_cell = get_unitcell_hull_centered([self.sample.ucell_a, self.sample.ucell_b, self.sample.ucell_c],
                                          [self.sample.ucell_alpha, self.sample.ucell_beta, self.sample.ucell_gamma])
            fig, ax = plot_unitcell(unit_cell, self.current_rotation, ax)
            if not None in [self.sample.kp, self.sample.k0]:
                ax = plot_momentum_transfer(self.sample, np.array([0,0,0]), ax)
                ax = plot_detector_target(self.sample, np.array([0,0,0]), self.diffractometer.nu,
                                      self.diffractometer.delta, ax)

            plt.show(block=False)

    def do_save_summary_to_txt(self, filename):
        '''Summarize current sample/calculations and save as .txt file with given filename.'''
        current_time = datetime.now()
        dt_string = current_time.strftime("%d/%m/%Y %H:%M:%S")

        # if no filename supplied, create one using datetime info
        if filename is None or "":
            filename = current_time.strftime("geocalc_outputs_%d-%m-%Y.txt")
            if os.path.exists(filename):
                filename = filename.removesuffix(".txt")
                filename = filename + current_time.strftime("_%Hh%M.txt")
        elif not filename.endswith(".txt"):
            filename = filename + ".txt"

        # make and write to txt file
        file = open(filename, "w+")

        file.write("LCLSii geometry calculator outputs "+dt_string+"\n")
        file.write("-----------------------------------\r\n")
        if self.sample.lattice_filepath is not None:
            file.write("Sample data referenced from: "+self.sample.lattice_filepath+"\r\n")
        file.write("-----------------------------------\n")
        file.write("UNIT CELL GEOMETRY\n")
        file.write("'a' side in line with lab 'x' vector before rotation.\n")
        file.write("alpha: angle between (b,c); beta: (a,c); gamma: (a,b)\n")
        file.write("Cell lengths in Angstrom:\n")
        file.write("a: ", str(self.sample.ucell_a), " | ",
              "b: ", str(self.sample.ucell_b), " | ",
              "c: ", str(self.sample.ucell_c), "\n")
        file.write("Cell angles in degrees:\n")
        file.write("alpha: ", str(np.rad2deg(self.sample.ucell_alpha)), " | ",
              "beta: ", str(np.rad2deg(self.sample.ucell_beta)), " | ",
              "gamma: ", str(np.rad2deg(self.sample.ucell_gamma)),"\r\n")
        file.write("------------------------------\n")
        file.write("Photon energy: ", str(self.sample.beam_energy), " eV\n")
        file.write("Photon wavelength: ", str(self.sample.lambda0), " Angstrom\r\n")
        file.write("------------------------------\n")
        file.write("Miller indices\n")
        file.write("H: ", str(self.sample.bragg_hkl[0]), " | ",
              "K: ", str(self.sample.bragg_hkl[1]), " | ",
              "L: ", str(self.sample.bragg_hkl[2]),"\r\n")
        file.write("------------------------------\n")
        if self.last_calc_type is not None:
            file.write("Calculation type: {} at angle {} deg.\n".format(self.last_calc_type, self.fixed_angle))
            file.write("Phi: {}\n".format(str(np.rad2deg(self.diffractometer.phi))))
            file.write("Kappa: {}\n".format(str(np.rad2deg(self.diffractometer.kappa))))
            file.write("Omega: {}\n".format(str(np.rad2deg(self.diffractometer.omega))))
            file.write("Detector nu: {}\n".format(str(np.rad2deg(self.diffractometer.nu))))
            file.write("Detector delta: {}\r\n".format(str(np.rad2deg(self.diffractometer.delta))))
            file.write("Momentum Transfer\n")
            file.write("k0: {}\n".format(np.array2string(self.sample.k0)))
            file.write("kp: {}\n".format(np.array2string(self.sample.kp)))
            file.write("------------------------------\n")

        file.close()

GeoCalculatorShell().cmdloop()
