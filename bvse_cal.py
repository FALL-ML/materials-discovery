'''
__author__ = "Shuting Chi"
__copyright__ = "Copyright 2019,"
__version__ = "1.0"
__maintainer__ = "Shuting Chi"
__email__ = "shu_cst@163.com"
__date__ = "Mar 19, 2019"
'''
import ase.io
import numpy as np
import os, sys

import Structure as struct
import BVAnalysis

def bv_calculation(dict, cif_file):
    # atoms = ase.io.read(cif_file, store_tags=True, onduplicates='keep')
    atoms = ase.io.read(cif_file, store_tags=True)
    struc = struct.Structure()
    struc.GetAseStructure(atoms)
    bvs = BVAnalysis.BVAnalysis()
    bvs.SetStructure(struc)
    bvs.SetMoveIon(dict['transport_ion'])
    bvs.ValenceOfMoveIon = dict['valence']
    bvs.SetLengthResolution(dict['resolution'])
    bvs.CaluBVSE(None)

    bvse_file = ''
    for k in dict:
        bvse_file = bvse_file + str(dict[k])+ '_'
    bvse_file=bvse_file.strip('_')
    # bvse_file = os.path.join(os.path.dirname(cif_file), bvse_file)
    bvse_file = cif_file.split(".")[0] + bvse_file

    bvs.SaveBVSData(filename=bvse_file)
    bvs.SaveBVSEData(filename=bvse_file)
    bvs.SaveBVELData(filename=bvse_file)
    return bvs.Ea, bvse_file


