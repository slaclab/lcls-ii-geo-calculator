import sys

from tools.cif_file_tools import *

testpath1 = "../CIFfiles/9002806.cif"
testpath2 = "../CIFfiles/AMS_DATA.cif"

dict1 = man_read_cif_to_dict(testpath1)
brav1 = cif_symm_to_bravais_lattice(dict1)
print(brav1)
sg1 = cif_get_spacegroup(dict1)
print(sg1)

dict2 = man_read_cif_to_dict(testpath2)
brav2 = cif_symm_to_bravais_lattice(dict2)
print(brav2)
sg2 = cif_get_spacegroup(dict2)
print(sg2)

# todo figure out how to check bravais lattices

