"""
PWscfInput class was design to speed up
input-generation processes for QUANTUM ESPRESSO
runs.

***Required***
Atomic Simulation Environment (ASE), numpy

***How to use***
If you know ASE, you know how to use this class.
See the provided example (make_pwscf_input.py) for details.

***Note***
This class is distributed in the hope
that it will benefit QE and/or ASE users.
There is NO WARRANTY of any form. Users are supposed
to carefully check the class and its obtained results.

***History***
+++ The python class was started back in 2010.
+++ Its author did not keep track of what have been changed.
+++ The code was last modified on 06/19/2013 for cleaning
up before being distributed. It was renamed to PWscfInput.


****************************************
Author: Duy Le
Department of Physics
University of Central Florida

Email: duy.le@ucf.edu
Website: http://www.physics.ufc.edu/~dle
****************************************
"""

import numpy as np
from ase.constraints import FixAtoms, FixScaled


class PWscfInput:
    def __init__(self, atoms):
        """
        Used multiple sub-classes

        """
        self.atoms = atoms
        self.control = Control()
        self.system = System(atoms)
        self.starting_magnetization \
            = StartingMagnetization \
            (self.system.structure.ntyp)
        self.electrons = Electrons()
        self.ions = Ions()
        self.cell = Cell()
        self.atomic_species = AtomicSpecies(self.atoms)
        self.kpoints = Kpoints()

    def write_input(self, filename):
        write_pwscf_input(self, filename)


class Control:
    def __init__(self):
        """
        &CONTROL section
        """
        self.settings = ControlSettings()
        self.io = ControlIO()
        self.ion_relax = ControlIonRelax()


class ControlSettings:
    def __init__(self):
        self.calculation = 'relax'
        self.restart_mode = 'from_scratch'
        self.prefix = 'pwscf'
        self.pseudo_dir = 'PATH_TO_PSEUDO_DIR'


class ControlIO:
    def __init__(self):
        self.outdir = '.'
        self.verbosity = 'low'


class ControlIonRelax:
    def __init__(self):
        self.tprnfor = True
        self.tstress = True
        self.forc_conv_thr = 1.0e-3
        self.etot_conv_thr = 1.0e-4
        self.nstep = 100


class System:
    def __init__(self, atoms):
        self.structure = SystemStructure(atoms)
        self.ecut = SystemEcut()
        self.occupations = SystemOccupations()
        self.spin_pol = SystemSpinPol()


class SystemStructure:
    def __init__(self, atoms):
        self.ibrav = 0
        self.a = 1.0
        self.nat = len(atoms)
        mol = get_reduce_atom_list(atoms)
        self.ntyp = len(mol)
        # incorrect if atoms are not grouped.
        # self.nbnd= 100


class SystemEcut:
    def __init__(self):
        self.ecutwfc = 50
        self.ecutrho = 200


class SystemOccupations:
    def __init__(self):
        self.occupations = 'smearing'
        self.smearing = 'fd'
        self.degauss = 0.007


class SystemSpinPol:
    def __init__(self):
        self.nspin = 1


class StartingMagnetization:
    def __init__(self, ntyp):
        starting_magnetization \
            = np.zeros([ntyp])
        self.starting_magnetization \
            = starting_magnetization


class Electrons:
    def __init__(self):
        self.diagonalization \
            = 'david'
        self.conv_thr = 1.0e-8
        self.electron_maxstep \
            = 100
        self.mixing_mode = 'plain'
        self.mixing_beta = 0.4


class Ions:
    def __init__(self):
        self.ion_dynamics = 'bfgs'
        self.wfc_extrapolation \
            = 'none'
        self.pot_extrapolation \
            = 'none'


class Cell:
    def __init__(self):
        self.cell_dynamics = 'bfgs'
        self.press = '0.0'
        self.press_conv_thr = '0.5'


class AtomicSpecies:
    def __init__(self, atoms):
        mol = get_reduce_atom_list(atoms)
        ntyp = len(mol)
        self.ntyp = ntyp
        self.symbol = mol
        mass = np.zeros([ntyp], dtype=np.float)
        self.mass = mass
        pseudo_potential = np.zeros([ntyp], dtype='S30')
        pseudo_potential[:] = 'Please_set_pseudo_file'
        self.pseudo_potential = pseudo_potential


class Kpoints:
    def __init__(self):
        self.type = 'automatic'
        self.nk = [15, 15, 1]
        self.sk = [0, 0, 0]


def get_reduce_atom_list(atoms):
    """
    Get list of atomic symbol then reduce it.
    New version of ASE should have this option
    already.
    """
    mol = atoms.get_chemical_symbols()
    if len(mol) > 1:
        for i in range(len(mol) - 1, 0, -1):
            for j in range(i):
                if mol[j] == mol[i]:
                    del mol[i]
                    break
    return mol


def write_k_points(kpoints, f):
    f.write('K_POINTS ' + kpoints.type + '\n')
    f.write("%4i %4i %4i %3i %3i %3i" % (
        kpoints.nk[0], kpoints.nk[1], kpoints.nk[2],
        kpoints.sk[0], kpoints.sk[1], kpoints.sk[2]) + '\n')


def write_key(item, dict):
    value = vars(dict)[item]
    if type(value) == str:
        str_value = '\'' + value + '\','
    if type(value) == float:
        str_value = str(value) + ','
    if type(value) == int:
        str_value = str(value) + ','
    if type(value) == bool:
        if (value):
            str_value = '.true.,'
        else:
            str_value = '.false.,'
    item_len = item.__len__()
    default_len = 25
    add_str = ''
    for i in range(item_len, default_len):
        add_str += ' '
    string = item + add_str + ' = ' + str_value
    return string


def write_array_key(item, dict, f):
    array_value = vars(dict)[item]
    for i in range(len(array_value)):
        value = array_value[i]
        item_len = item.__len__()
        default_len = 25
        add_str = ''
        for j in range(item_len + 3, default_len):
            add_str += ' '

        string = item + '(' + str(i + 1) + ')' + add_str + ' = '
        f.write(string + str(value) + ',\n')


def write_atomic_species(atomic_species, f):
    for i in range(atomic_species.ntyp):
        f.write("%5s %8.4f %s" % ( \
            atomic_species.symbol[i], \
            atomic_species.mass[i], \
            atomic_species.pseudo_potential[i]
        )+'\n')


def write_structure(atoms, f):
    f.write('ATOMIC_POSITIONS  crystal \n')
    sflags = np.zeros((len(atoms), 3), dtype=bool)
    newsflags = np.ones((len(atoms), 3), dtype=np.int)
    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixScaled):
                sflags[constr.a] = constr.mask
            elif isinstance(constr, FixAtoms):
                sflags[constr.index] = [True, True, True]
            for i in range(len(atoms)):
                for j in range(3):
                    if sflags[i, j]:
                        newsflags[i, j] = 0

    for i in range(len(atoms)):
        f.write('%3s %20.14f %20.14f %20.14f %3i %3i %3i' % ( \
            atoms.get_chemical_symbols()[i], \
            atoms.get_scaled_positions()[i, 0], \
            atoms.get_scaled_positions()[i, 1], \
            atoms.get_scaled_positions()[i, 2], \
            newsflags[i, 0], newsflags[i, 1], newsflags[i, 2])+'\n')

    f.write('CELL_PARAMETERS\n')
    for i in range(3):
        f.write('%20.14f %20.14f %20.14f' % ( \
            atoms.cell[i, 0], \
            atoms.cell[i, 1], \
            atoms.cell[i, 2]) + '\n')


def write_pwscf_input(object, filename):
    f = open(filename, 'w')
    """ Write CONTROL section """
    f.write('&CONTROL\n')
    f.write('! .control.settings.\n')
    dict = object.control.settings
    for item in vars(dict):
        f.write(write_key(item, dict)+'\n')

    f.write('\n')
    f.write('! .control.io.\n')
    dict = object.control.io
    for item in vars(dict):
        f.write(write_key(item, dict)+'\n')

    f.write('\n')
    f.write('! .control.ion_relax.\n')
    dict = object.control.ion_relax
    for item in vars(dict):
        f.write(write_key(item, dict) + '\n')

    f.write('/\n')

    f.write('\n')

    """ &SYSTEM section """
    f.write('&SYSTEM\n')
    f.write('! .system.structure.\n')
    dict = object.system.structure
    for item in vars(dict):
        f.write(write_key(item, dict)+'\n')

    f.write('\n')
    f.write('! .system.ecut.\n')
    dict = object.system.ecut
    for item in vars(dict):
        f.write(write_key(item, dict)+'\n')

    f.write('\n')
    f.write('! .system.occupations.\n')
    dict = object.system.occupations
    for item in vars(dict):
        f.write(write_key(item, dict) + '\n')

    if object.system.spin_pol.nspin == 2:
        f.write('\n')
        f.write('! .system.spin_pol.\n')
        dict = object.system.spin_pol
        for item in vars(dict):
            f.write(write_key(item, dict) + '\n')

        f.write('! .system.starting_magnetization.\n')
        dict = object.starting_magnetization
        for item in vars(dict):
            write_array_key(item, dict, f)

    f.write('/\n')
    f.write('\n')

    """ &ELECTRONS section """
    f.write('&ELECTRONS\n')
    f.write('! .electrons.\n')
    dict = object.electrons
    for item in vars(dict):
        f.write(write_key(item, dict) + '\n')

    f.write('/\n')
    f.write('\n')

    """ &IONS section """
    f.write('&IONS\n')
    f.write('! .ions.\n')
    dict = object.ions
    for item in vars(dict):
        f.write(write_key(item, dict) + '\n')

    f.write('/\n')
    f.write('\n')

    """ &Cell section """
    f.write('&CEll\n')
    f.write('! .cell.\n')
    dict = object.cell
    for item in vars(dict):
        f.write(write_key(item, dict) + '\n')

    f.write('/\n')
    f.write('\n')

    """ ATOMIC_SPECIES section """
    f.write('ATOMIC_SPECIES\n')
    write_atomic_species(object.atomic_species, f)
    f.write('\n')
    write_structure(object.atoms, f)
    f.write('\n')
    write_k_points(object.kpoints, f)