import os
import re
import ase.io
from Structure import Structure as customstruc
import BVAnalysis

def neb_bvse(filename, transport_ion="Li", valence=1, mobile_ion_interact=False):
    bvse_condition_key = {'origin': 0, 'origin_id': 0, 'transport_ion': transport_ion, 'valence': valence,
                          'resolution': 0.1}
    config_dict = {'origin': 0, 'origin_id': 0, 'transport_ion': transport_ion, 'valence': valence, 'resolution': 0.1}
    bvse_condition = {}
    for k in bvse_condition_key:
        bvse_condition[k] = config_dict[k]
    atoms = ase.io.read(os.path.join('./NEB_BVSE_Comparison', filename), store_tags=True)
    if mobile_ion_interact:
        for atom in range(atoms.positions.shape[0]):
            if atom == atoms.positions.shape[0]-1:
                continue
            else:
                if atoms.info['_atom_site_type_symbol'][atom] == "Li+":
                    atoms.info['_atom_site_type_symbol'][atom] = "A+"
                    atoms.info['_atom_site_label'][atom] = 'A' + \
                                                           re.findall(r'(\w+?)(\d+)', atoms.info['_atom_site_label'][atom])[
                                                               0][-1]
                    atoms.numbers[atom] = 0
        atoms.info['_atom_type_symbol'].append('A+')
        atoms.info['_atom_type_oxidation_number'].append(
            atoms.info['_atom_type_oxidation_number'][atoms.info['_atom_type_symbol'].index('Li+')])
        interact_label = '_' + transport_ion + '_interaction'
    else:
        interact_label = ''
    struc = customstruc()
    struc.GetAseStructure(atoms)
    bvs = BVAnalysis.BVAnalysis()
    bvs.SetStructure(struc)
    bvs.SetMoveIon(bvse_condition['transport_ion'])
    bvs.ValenceOfMoveIon = bvse_condition['valence']
    bvs.SetLengthResolution(bvse_condition['resolution'])
    bvs.CaluBVSE(None)
    bvse_file = ''
    for k in bvse_condition:
        bvse_file = bvse_file + str(bvse_condition[k]) + '_'
    bvse_file = bvse_file.strip('_')
    bvse_file = '_'.join(filename.split('\\')[-1].split('.')[0].split('_')[:5]) + '_' + \
                filename.split('\\')[-1].split('.')[0].split('_')[-1] + interact_label + '_' + bvse_file
    bvs.SaveBVSData(filename="./NEB_BVSE_Comparison/" + bvse_file)
    bvs.SaveBVSEData(filename="./NEB_BVSE_Comparison/" + bvse_file)
    bvs.SaveBVELData(filename="./NEB_BVSE_Comparison/" + bvse_file)

