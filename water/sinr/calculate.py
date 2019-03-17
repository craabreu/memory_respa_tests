import argparse
import atomsmm
import pandas as pd

import math
import time

from sys import stdout
from simtk import openmm
from simtk import unit
from simtk.openmm import app

parser = argparse.ArgumentParser()
parser.add_argument('--device', dest='device', help='the GPU device', default='None')
parser.add_argument('--secdev', dest='secdev', help='the secondary GPU device', default='None')
parser.add_argument('--seed', dest='seed', help='the RNG seed', type=int, default=0)
parser.add_argument('--platform', dest='platform', help='the computation platform', default='CUDA')
args = parser.parse_args()

seed = int(1000*time.time()) % 16384 if args.seed == 0 else args.seed
print(f'Employed RNG seed is {seed}')

# base = f'sinr-{args.timestep:02d}fs'
platform_name = args.platform

rcut = 12*unit.angstroms
rswitch = 11*unit.angstroms

platform = openmm.Platform.getPlatformByName(platform_name)
properties = dict(Precision='mixed') if platform_name == 'CUDA' else dict()
if args.device != 'None':
    properties['DeviceIndex'] = args.device

pdb = app.PDBFile('water.pdb')
forcefield = app.ForceField('water.xml')
openmm_system = forcefield.createSystem(pdb.topology,
                                        nonbondedMethod=openmm.app.PME,
                                        nonbondedCutoff=rcut,
                                        rigidWater=False,
                                        removeCMMotion=False)

nbforce = openmm_system.getForce(atomsmm.findNonbondedForce(openmm_system))
nbforce.setUseSwitchingFunction(True)
nbforce.setSwitchingDistance(rswitch)
nbforce.setUseDispersionCorrection(True)

integrator = openmm.CustomIntegrator(0)

simulation = openmm.app.Simulation(pdb.topology, openmm_system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

dataReporter = atomsmm.ExtendedStateDataReporter(stdout, 1, separator=',',
        step=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        speed=True)

simulation.reporters += [dataReporter]
simulation.step(1)
