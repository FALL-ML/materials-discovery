import os, sys, re
import copy
from bvse_cal import bv_calculation
from pymatgen.io.cif import CifWriter

from monty.io import zopen

from pymatgen.io.cif import CifParser
from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation
from pymatgen.transformations.standard_transformations import DiscretizeOccupanciesTransformation
from pymatgen.transformations.standard_transformations import ConventionalCellTransformation
from pymatgen.transformations.standard_transformations import OxidationStateDecorationTransformation

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen import symmetry
from pymatgen import Structure
import pymatgen.io.ase as ase_io
import pymatgen.io.cif as cif_io

from ase.build import bulk
from ase.io import espresso
import ase.io

import Structure as customstruc
import BVAnalysis

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import itertools
from collections import OrderedDict
import pwscf_input

import pandas as pd

import json

import pymatgen.core.periodic_table as pt

from multiprocessing import Pool
from functools import partial
import worker_neb_bvse

from cavd.netio import *
from cavd.channel import Channel
from cavd.netstorage import AtomNetwork, connection_values_list
from cavd.local_environment import CifParser_new, LocalEnvirCom
from cavd.get_Symmetry import get_symnum_sites, get_equivalent_vornet,get_labeled_vornet
from cavd.recovery import rediscovery, rediscovery_kdTree, rediscovery_byRad_kdTree, rediscovery_byRad_kdTree_onlyVertex
from cavd import bmd_com
from cavd_bvse.mergecluster import load_voids_channels_from_file, load_struc, load_bvse

from cavd_bvse.non_equivalent_paths import non_equivalent_paths

import tkinter as tk
from tkinter import filedialog

from shutil import copyfile


# Default settings - manually set the first two: conducting_ion and cif_filename
mobile_ion = "Li"
cif_filename = 'mp-1222492.cif'
root_name_of_file = cif_filename.replace(".cif", "")

# Automatically generate root directories
root_dictionary = dict(
    cif_input_root=os.path.join(os.getcwd(), '{}_cifs\\'.format(mobile_ion.lower())),
    cif_output_root=os.path.join(os.getcwd(), '{}_cifs_sanitized\\'.format(mobile_ion.lower())),
    cavd_output_root=os.path.join(os.getcwd(), '{}_cavd_outputs\\'.format(mobile_ion.lower())),
    bvse_output_root=os.path.join(os.getcwd(), '{}_bvse_outputs\\'.format(mobile_ion.lower())),
    scf_input_root=os.path.join(os.getcwd(), '{}_scf_inputs\\'.format(mobile_ion.lower())),
    scf_output_root=os.path.join(os.getcwd(), '{}_scf_outputs\\'.format(mobile_ion.lower())),
    vcr_input_root=os.path.join(os.getcwd(), '{}_vcr_inputs\\'.format(mobile_ion.lower())),
    vcr_output_root=os.path.join(os.getcwd(), '{}_vcr_outputs\\'.format(mobile_ion.lower())),
    nep_output_root=os.path.join(os.getcwd(), '{}_nep_outputs\\'.format(mobile_ion.lower())),
    neb_root=os.path.join(os.getcwd(), '{}_neb_cifs\\'.format(mobile_ion.lower())))

for path in root_dictionary.values():
    try:
        os.mkdir(path)
    except:
        pass

# Automatically generate structure-specific folders
paths_dictionary = dict(
    cif_input_dir=os.path.join(os.getcwd(), '{}_cifs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    cif_output_dir=os.path.join(os.getcwd(), '{}_cifs_sanitized\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    cavd_output_dir=os.path.join(os.getcwd(), '{}_cavd_outputs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    bvse_output_dir=os.path.join(os.getcwd(), '{}_bvse_outputs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    scf_input_dir=os.path.join(os.getcwd(), '{}_scf_inputs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    scf_output_dir=os.path.join(os.getcwd(), '{}_scf_outputs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    vcr_input_dir=os.path.join(os.getcwd(), '{}_vcr_inputs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    vcr_output_dir=os.path.join(os.getcwd(), '{}_vcr_outputs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    nep_output_dir=os.path.join(os.getcwd(), '{}_nep_outputs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)),
    neb_dir=os.path.join(os.getcwd(), '{}_neb_cifs\\{}\\'.format(mobile_ion.lower(), root_name_of_file)))


for path in paths_dictionary.values():
    try:
        os.mkdir(path)
    except:
        pass

def ordering_check(filename=None, cif_input_dir=paths_dictionary['cif_input_dir'],
                   cif_output_dir=paths_dictionary['cif_output_dir']):
    # look for the file in the Li Cifs/ folder
    filepath = cif_input_dir + filename

    #     stru = Structure.from(filepath)

    # pull in structure
    with zopen(filepath, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    print("Spacegroup of input: {}".format(stru.get_space_group_info()))

    # if the structure is ordered then just write the .cif to a "sanitized" folder
    if (stru.is_ordered):
        # save raw version
        copyfile((cif_input_dir + filename), (cif_output_dir + "raw_" + filename))
        #         w = CifWriter(stru, symprec=None)
        #         w.write_file(cif_output_dir + "raw_" + filename)

        # save symprec version
        w = CifWriter(stru, symprec=True)
        w.write_file(cif_output_dir + "symprec_" + filename)
        print("Input is ordered")


    # if the structure is not ordered, attempt to order it
    else:
        try:
            # convert input into a conventional cell
            trans = ConventionalCellTransformation()
            stru = trans.apply_transformation(stru)

            # convert the conventional cell into an ordered cell
            trans = OrderDisorderedStructureTransformation()
            stru = trans.apply_transformation(stru, return_ranked_list=100)

            # use the structureMatcher to knock out equivalent sites
            matcher = StructureMatcher()
            groups = matcher.group_structures([d["structure"] for d in stru])

            # write the lowest energy structure to the "sanitized" folder
            groups[0][0].to(filename=cif_output_dir + "raw_" + filename)

            # save symprec version
            w = CifWriter(groups[0][0], symprec=True)
            w.write_file(cif_output_dir + "symprec_" + filename)

            print("Input sucessfully ordered")

        except Exception as e:
            # if conversion of an unordered structure fails:
            print("ERROR: Disordered and {}".format(e))


def outvesta_fcs(filename, migrant, ext="", ntol=0.02, lower=0, upper=10, symprec=0.01,
                 cif_output_dir=paths_dictionary['cif_output_dir'], \
                 cavd_output_dir=paths_dictionary['cavd_output_dir'],
                 bvse_output_dir=paths_dictionary['bvse_output_dir'], scaling=1):
    filepath = cif_output_dir + ext + filename

    with zopen(filepath, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]

    # attempt on a supercell
    if scaling != 1:
        stru.make_supercell([scaling, scaling, scaling])

    sitesym = parser.get_sym_opt()
    species = [str(sp).replace("Specie ", "") for sp in stru.species]
    elements = [re.sub('[^a-zA-Z]', '', sp) for sp in species]

    if migrant not in elements:
        raise ValueError("The input migrant ion not in the input structure! Please check it.")

    effec_radii, migrant_radius, migrant_alpha, nei_dises, coordination_list = LocalEnvirCom(stru, migrant)

    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filepath, migrant, effec_radii, True)
    vornet, edge_centers, fcs, faces = atmnet.perform_voronoi_decomposition(True, ntol)

    add_fcs_vornet = vornet.add_facecenters(faces)
    sym_vornet, voids = get_labeled_vornet(add_fcs_vornet, sitesym, symprec)
    voids_abs = []

    for void in sym_vornet.nodes:
        voids_abs.append(void[2])

    bottlenecks = []
    for bt in sym_vornet.edges:
        bottlenecks.append(bt[2])

    fcens = []
    for fc in fcs:
        fcens.append(fc[0])

    vorosites = [voids_abs, bottlenecks, fcens]
    recover_rate, recover_state, migrate_mindis = rediscovery_kdTree(stru, migrant, vorosites)
    prefixname = filename.replace(".cif", "")
    newpath = cavd_output_dir + ext + prefixname

    # compute the R_T
    conn_val = connection_values_list(newpath + ".resex", sym_vornet)
    channels = Channel.findChannels2(sym_vornet, atmnet, lower, upper, newpath + ".net")
    prefixname = filename.replace(cif_output_dir, "")
    prefixname = prefixname.replace(".cif", "")

    # output vesta file for visualization
    newpath = bvse_output_dir + ext + prefixname
    Channel.writeToVESTA(channels, atmnet, newpath)

    return conn_val, recover_rate, stru


def simple_bvse(filename, cif_output_dir=paths_dictionary['cif_output_dir'], ext="", \
                bvse_condition_key={'origin': 0, 'origin_id': 0, 'transport_ion': 'Li', 'valence': 1,
                                    'resolution': 0.1}, \
                config_dict={'origin': 0, 'origin_id': 0, 'transport_ion': 'Li', 'valence': 1, 'resolution': 0.1}):
    cif_filename = cif_output_dir + ext + filename

    bvse_condition = {}
    for k in bvse_condition_key:
        bvse_condition[k] = config_dict[k]
    Ea, bvse_file = bv_calculation(bvse_condition, cif_filename)

    bvse = Ea['BVSE']
    bvel = Ea['BVEL']
    bvse_dict = {'bvse_1d': bvse[0], 'bvse_2d': bvse[1], 'bvse_3d': bvse[2], 'bvel_1d': bvel[0], 'bvel_2d': bvel[1],
                 'bvel_3d': bvel[2]}
    bvse_filedir = bvse_file + '_BVSE.npy'

    return Ea, bvse_file, bvse_filedir, bvse


def cavd_bvse_stabilizer(filename, migrant, ntol, lower, upper, \
                         cif_output_dir=paths_dictionary['cif_output_dir'], \
                         cavd_output_dir=paths_dictionary['cavd_output_dir'], \
                         bvse_output_dir=paths_dictionary['bvse_output_dir']):
    # try to evaluate the raw file
    extension = "raw_"
    try:
        conn_val, recover_rate, structure = outvesta_fcs(filename, migrant, ext=extension, ntol=ntol, lower=lower,
                                                         upper=upper, \
                                                         cif_output_dir=cif_output_dir, \
                                                         cavd_output_dir=cavd_output_dir, \
                                                         bvse_output_dir=bvse_output_dir)
        Ea, bvse_file, bvse_filedir, bvse = simple_bvse(filename, ext=extension, cif_output_dir=cif_output_dir)
        return Ea, conn_val, recover_rate, extension, bvse_file, bvse_filedir, bvse

    except Exception as e:
        print("Raw structure incompatible with BVSE: {}".format(e))

    # if the raw file didn't work, try to evaluate the symprec file
    extension = "symprec_"
    try:
        recover_rates = []
        structures = []
        # Check three cells: [1x, 1x, 1x],  [2x, 2x, 2x], and [3x, 3x, 3x]
        # Recovery rates sometimes improve at larger supercells
        for i in [0, 1, 2]:
            # Below vales for lower and upper are specific to Li. Values are taken from CAVD reference paper.
            conn_val, recover_rate, structure = outvesta_fcs(filename, migrant, ext=extension, ntol=ntol, lower=lower,
                                                             upper=upper, \
                                                             cif_output_dir=cif_output_dir, \
                                                             cavd_output_dir=cavd_output_dir, \
                                                             bvse_output_dir=bvse_output_dir, scaling=(i + 1))
            recover_rates.append(recover_rate)
            structures.append(structure)
            if recover_rate == 1:
                break

        # find the id that had the maximum recovery rate and return that
        max_rate = max(recover_rates)
        max_idx = recover_rates.index(max_rate)
        w = CifWriter(structures[max_idx], symprec=True)
        w.write_file(cif_output_dir + extension + filename)
        conn_val, recover_rate, structure = outvesta_fcs(filename, migrant, ext=extension, ntol=ntol, lower=lower,
                                                         upper=upper, \
                                                         cif_output_dir=cif_output_dir, \
                                                         cavd_output_dir=cavd_output_dir, \
                                                         bvse_output_dir=bvse_output_dir)
        Ea, bvse_file, bvse_filedir, bvse = simpleBVSE(filename, ext=extension, cif_output_dir=cif_output_dir)
        return Ea, conn_val, recover_rate, extension, bvse_file, bvse_filedir, bvse

    except Exception as e:
        print("Symprec structure incompatible with BVSE: {}".format(e))

    # Return None for all values if no representations were compatible with simpleBVSE
    raise ValueError("No representations were compatible with simpleBVSE")



if __name__ == '__main__':


    # run orderingCheck
    # creates two potentially viable structures named: raw_+cif_filename.cif, and symprec_+cif_filename.cif
    ordering_check(cif_filename)

    # run cavd_bvse_stabilizer() on the input structure
    # the functon runs CAVD+BVSE on the two potentially viable structures, in order
    # for the symprec_ structure it creates supercells to try to optimize the recover_rate
    # returns results for whichever structure works first: raw_ > symprec_ (symprec is the last resort)
    # the extension that works is returned
    Ea, conn_val, recover_rate, extension, bvse_file, bvse_filedir, bvse = cavd_bvse_stabilizer(cif_filename, 'Li',
                                                                                                ntol=0.02, lower=0,
                                                                                                upper=10)
    print("CAVD ConnVal's are: {}".format(conn_val))
    print("BVSE Ea values are: {}".format(Ea))

    sanitized_cif_location = paths_dictionary['cif_output_dir'] + extension + cif_filename
    cavd_file_location = paths_dictionary['cavd_output_dir'] + extension + cif_filename.split('.')[0] + '.net'
    print(sanitized_cif_location)
    print(cavd_file_location)

    mep = non_equivalent_paths(sanitized_cif_location, bvse_filedir, cavd_file_location, bvse, mobile_ion)
    print("{} paths identified.".format(len(mep.paths)))