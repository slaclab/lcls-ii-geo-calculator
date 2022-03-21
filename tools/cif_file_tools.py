import numpy as np
from numpy import sin, cos, sqrt
from angle_calc_utils import *


def cif_symm_to_bravais_lattice(cif_file, use_gemmi=False):
    # get realspace(?) lattice vectors from cif symmetry information
    if (use_gemmi):
        pass
    else:
        # reading from manually generated dict
        try:
            a = float(cif_file['_cell_length_a'])
            b = float(cif_file['_cell_length_b'])
            c = float(cif_file['_cell_length_c'])
            alpha = np.radians(float(cif_file['_cell_angle_alpha']))
            beta = np.radians(float(cif_file['_cell_angle_beta']))
            gamma = np.radians(float(cif_file['_cell_angle_gamma']))
        except KeyError as e:
            print("The CIF file seems to be missing some cell length or angle data.")
            print(e)

        bravais = bravais_matrix(a, b, c, alpha, beta, gamma)

    return bravais, [a, b, c], [alpha, beta, gamma]


def man_read_cif_to_dict (relpath):
    # (partially) reads cif info to a dictionary using a relative filepath
    # manual file read (old, switching to gemmi)
    # todo cleanup multiline data/headers
    cif_dict = {}
    with open(relpath, 'rb') as f:
        for bline in f:
            line = bline.decode('ascii')
            if line[0] != '#': # disregard header
                #print(line)
                linearray = line.split()
                if len(linearray) > 1:
                    cif_dict[linearray[0]] = ' '.join(linearray[1:])
    f.close()

    return cif_dict

def cif_to_symm_type(cif_dict):
    # get the symmetry type/group classification from the cif dict
    try:
        itnum = int(cif_dict['_space_group_IT_number'])
    except KeyError as e:
        itnum = None
        print("can't read group IT number from cif dictionary")
        print(e)

    # group number classifications taken from https://en.wikipedia.org/wiki/List_of_space_groups
    if itnum <= 2:
        symm_type = "triclinic"
    elif itnum <= 15:
        symm_type = "monoclinic"
    elif itnum <= 74:
        symm_type = "orthorhombic"
    elif itnum <= 142:
        symm_type = "tetragonal"
    elif itnum <= 167:
        symm_type = "trigonal"
    elif itnum <= 194:
        symm_type = "hexagonal"
    elif itnum <= 230:
        symm_type = "cubic"
    else:
        print("group IT number", itnum, " out of bounds!")
        symm_type = "n/a"

    return symm_type


# functions below are old: trying to use gemmi package to read cif file
# gemmi could theoretically be a more robust cif reader than the manual one (which relies on reading raw text)
# implemented here, but was struggling to get it to function
def read_cif(pathname):

    try:
        pass
        # f = cif.read(pathname)
    except [KeyError, ValueError] as e:
        print("Could not read CIF file from ", pathname)
        print(e)
        f = None

    return f


def cif_get_spacegroup(cif_dict):
    # returns a gemmi "spacegroup" object using the HM name from the cif file
    try:
        hm_name = cif_dict['_symmetry_space_group_name_H-M']
        print(hm_name)
    except KeyError as e:
        print('Could not find Hermann Mauguin name in cif file')

    try:
        pass
        alpha = (float(cif_dict['_cell_angle_alpha']))
        #gamma = float(cif_dict['_cell_angle_gamma'])
        #sg = gemmi.SpaceGroup(hm_name)#, alpha = alpha, gamma = gamma)

        return None #sg
    except ValueError as e:
        print(e)
        print("trying to look up via number...")

    try:
        pass
        group_no = int(cif_dict["_space_group_IT_number"])
        #sg = gemmi.find_spacegroup_by_number(group_no)
        return None #sg
    except Exception as e:
        print(e)
        print("could not look up by number")

    return None





