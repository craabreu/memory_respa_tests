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
parser.add_argument('--timestep', dest='timestep', help='time step size', type=int, required=True)
parser.add_argument('--nsteps', dest='nsteps', help='number of steps', type=int, required=True)
parser.add_argument('--device', dest='device', help='the GPU device', default='None')
parser.add_argument('--seed', dest='seed', help='the RNG seed', type=int, default=0)
parser.add_argument('--platform', dest='platform', help='the computation platform', default='CUDA')
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

dt = args.timestep*unit.femtoseconds
temp = 298.15*unit.kelvin
rcut = 12*unit.angstroms
rswitch = 11*unit.angstroms
rcutIn = 8*unit.angstroms
rswitchIn = 5*unit.angstroms
tau = 10*unit.femtoseconds
gamma = 0.1/unit.femtoseconds
reportInterval = 90//args.timestep

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
                                           softcore_group=1,
                                           split_exceptions=True)

if args.timestep <= 3:
    respa_system = solvation_system
    for force in solvation_system.getForces():
        if isinstance(force, openmm.NonbondedForce):
            force.setReciprocalSpaceForceGroup(1)
            force.setForceGroup(1)
    loops = [2*args.timestep, 1]
else:
    respa_system = atomsmm.RESPASystem(solvation_system, rcutIn, rswitchIn, fastExceptions=True)
    loops = [6, args.timestep//3, 1]

if method == 'Langevin':
    integrator = atomsmm.integrators.Langevin_R_Integrator(dt, loops, temp, gamma, has_memory=True)
else:
    integrator = atomsmm.SIN_R_Integrator(dt, loops, temp, tau, gamma)
integrator.setRandomNumberSeed(seed)

simulation = openmm.app.Simulation(pdb.topology, respa_system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temp, seed)

states_data = pd.read_csv('alchemical_states.inp', sep='\s+', comment='#')
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
