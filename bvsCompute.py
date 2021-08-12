'''
__author__ = "Shuting Chi"
__copyright__ = "Copyright 2019,"
__version__ = "1.0"
__maintainer__ = "Shuting Chi"
__email__ = "shu_cst@163.com"
__date__ = "Mar 19, 2019"
'''
import bvse_cal, datetime, cif_preprocess
import ase.io, os, sys, traceback, json
from interact import *

#mkdir launch_dir according to the config_file
def launch_dir(config_file):
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    os.chdir(BASE_DIR)

    if not os.path.exists(os.path.join(BASE_DIR, 'launches')):
        os.mkdir(os.path.join(BASE_DIR, 'launches'))
    new_dir = os.path.splitext(os.path.basename(config_file))[0]
    new_dir = os.path.join(os.path.join(BASE_DIR, 'launches'), new_dir)
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    return new_dir

def bvsCompute(config_file):
    try:
        with open(config_file, "r") as f:
            config_dict=json.load(f)
            f.close()

        new_dir = launch_dir(config_file)

        inte = interact_mongodb()
        inte.set_identifier(config_dict['origin'], config_dict['origin_id'])
        cif_id=config_dict['structid']
        
        cif_filename = os.path.join(new_dir, config_dict['origin']+'_'+str(config_dict['origin_id'])+'.cif')
        if inte.download(cif_id, 'cif', cif_filename):
            try:
                cif_preprocess.check(cif_filename, config_dict['transport_ion'], config_dict['isorder'])
            except Exception as err:
                inte.insert_to_bvse(None, {'error': str(err)}, {'structid': config_dict['structid'], 'origin':config_dict['origin'], 'origin_id': config_dict['origin_id']})
            else:
                Ea, bvse_file = bvse_cal.bv_calculation(config_dict, cif_filename)
                BVSE = Ea['BVSE']
                BVEL = Ea['BVEL']
                bvse_dict={'bvse_1d': BVSE[0], 'bvse_2d': BVSE[1], 'bvse_3d': BVSE[2], 'bvel_1d': BVEL[0], 'bvel_2d': BVEL[1], 'bvel_3d': BVEL[2]}

                config_dict.pop('isorder')
                inte.insert_to_bvse(bvse_file, bvse_dict, config_dict)
    except:
        raise Exception('BVSE calculation failed.')

if __name__ == '__main__':
    bvsCompute(sys.argv[1])
