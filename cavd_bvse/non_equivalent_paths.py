'''
__author__ = "Penghui Mi"
__copyright__ = "Copyright 2019,"
__version__ = "1.0"
__maintainer__ = "Penghui Mi"
__email__ = "phmim@shu.edu.cn"
__date__ = "Mar 19, 2019"
'''
from cavd_bvse.selectvoidbyBVSE import SelectVoidByBVSE
from cavd_bvse.mergecluster import MergeCluster
from cavd_bvse.cal_all_paths import MigrationPaths
from pymatgen.core.structure import Structure
from cavd_bvse.mergecluster import load_voids_channels_from_file, load_struc, load_bvse


def non_equivalent_paths(filename_CIF, filename_BVSE, filename_CAVD, bv_Ea, moveion='Li'):
    voids, channels = load_voids_channels_from_file(filename_CAVD)
    struc = load_struc(filename_CIF)
    energy = load_bvse(filename_BVSE)

    mc = MergeCluster(voids, channels, struc, energy, filename_CIF)
    mc.to_net(filename_CAVD)
    svbb = SelectVoidByBVSE(struc, energy)
    svbb.setpara(mc.mergedvoids, mc.mergedchannels, bv_Ea)
    svbb.selectvoid()
    svbb.calchannellabel()
    svbb.to_net(filename_CIF)

    mp = MigrationPaths(struc, energy, moveion)
    mp.setvoids(svbb.voids)
    mp.setchannels(svbb.channels)
    mp.bulid_migrationnet()
    mp.read_cif_nodles()
    mp.find_endpoints()
    mp.cal_equal_path()
    mp.calpathsenergy()
    # mp.showenergy(filename_CIF)
    mp.savedata(filename_CIF)
    return mp

def showpaths(paths,filename):
    struc = Structure.from_str(open(filename).read(), fmt="cif")
    struc.to(fmt='POSCAR', filename=filename.split(".")[0] + '_POSCAR')
    struc1 = Structure.from_file(filename.split(".")[0] + '_POSCAR')
    for i in range(len(paths)):
        for j in range(1, len(paths[i]) - 1):
            struc1.insert(0, 'He', paths[i][j])
    struc1.to(filename=filename.split(".")[0] + '_POSCARpath')

