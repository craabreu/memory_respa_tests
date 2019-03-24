import argparse
import atomsmm
import numpy as np
import pandas as pd

import math
import time

from sys import stdout
from simtk import openmm
from simtk import unit
from simtk.openmm import app

parser = argparse.ArgumentParser()
#parser.add_argument('--timestep', dest='timestep', help='time step size', type=int, required=True)
parser.add_argument('--nsteps', dest='nsteps', help='number of steps', type=int, required=True)
parser.add_argument('--device', dest='device', help='the GPU device', default='None')
parser.add_argument('--seed', dest='seed', help='the RNG seed', type=int, default=0)
parser.add_argument('--platform', dest='platform', help='the computation platform', default='CUDA')
parser.add_argument('--part', dest='part', help='the case part', type=int, default=0)
args = parser.parse_args()

seed = int(1000*time.time()) % 16384 if args.seed == 0 else args.seed
print(f'Employed RNG seed is {seed}')

solute = 'phenol'
solvent = 'water'
base = '{}-in-{}'.format(solute, solvent)
platform_name = args.platform
methods = {1:'Langevin', 2:'SIN-R'}
method = methods[1]

steps_per_state = args.nsteps

dt = 0.5*unit.femtoseconds
temp = 298.15*unit.kelvin
rcut = 12*unit.angstroms
rswitch = 11*unit.angstroms
rcutIn = 8*unit.angstroms
rswitchIn = 5*unit.angstroms
tau = 10*unit.femtoseconds
gamma = 0.1/unit.femtoseconds
reportInterval = 40

platform = openmm.Platform.getPlatformByName(platform_name)
properties = dict(Precision='mixed') if platform_name == 'CUDA' else dict()
if args.device != 'None':
    properties['DeviceIndex'] = args.device

pdb = app.PDBFile(f'{base}.pdb')
residues = [atom.residue.name for atom in pdb.topology.atoms()]
solute_atoms = set(i for (i, name) in enumerate(residues) if name == 'aaa')

forcefield = app.ForceField(f'{base}.xml')
openmm_system = forcefield.createSystem(pdb.topology,
                                        nonbondedMethod=openmm.app.PME,
                                        nonbondedCutoff=rcut,
                                        rigidWater=False,
                                        removeCMMotion=False)

nbforce = openmm_system.getForce(atomsmm.findNonbondedForce(openmm_system))
nbforce.setUseSwitchingFunction(True)
nbforce.setSwitchingDistance(rswitch)
nbforce.setUseDispersionCorrection(True)

solvation_system = atomsmm.SolvationSystem(openmm_system, solute_atoms,
                                           use_softcore=True,
                                           softcore_group=0,
                                           split_exceptions=True)

respa_system = solvation_system
loops = [1]

if method == 'Langevin':
    integrator = atomsmm.integrators.Langevin_R_Integrator(dt, loops, temp, gamma, has_memory=True)
else:
    integrator = atomsmm.SIN_R_Integrator(dt, loops, temp, tau, gamma)
integrator.setRandomNumberSeed(seed)

simulation = openmm.app.Simulation(pdb.topology, respa_system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temp, seed)

states_file = 'alchemical_states' + ('.inp' if args.part == 0 else f'_{args.part}.inp')
states_data = pd.read_csv(states_file, sep='\s+', comment='#')
parameterStates = states_data[['lambda_vdw', 'lambda_coul']]
simulate = states_data['weight'] != -np.inf
for state in reversed(states_data.index):
    if simulate[state]:
        for name, value in parameterStates.iloc[state].items():
            simulation.context.setParameter(name, value)
            print(f'{name} = {value}')
        dataReporter = atomsmm.ExtendedStateDataReporter(stdout, reportInterval, separator=',',
            step=True, potentialEnergy=True, temperature=True,
            speed=True, extraFile=f'{base}_data-{state:02d}.csv')
        multistateReporter = atomsmm.ExtendedStateDataReporter(f'{base}_energy-{state:02d}.csv',
            reportInterval, separator=',', step=True, globalParameterStates=parameterStates)
        simulation.reporters = [dataReporter, multistateReporter]
        simulation.step(steps_per_state)
