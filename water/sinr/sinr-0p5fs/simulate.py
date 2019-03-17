import argparse
import atomsmm
import pandas as pd

import math
import time

from sys import stdout
from simtk import openmm
from simtk import unit
from simtk.openmm import app

import os
numThreads = 8
os.putenv('OPENMM_CPU_THREADS', str(numThreads))

parser = argparse.ArgumentParser()
parser.add_argument('--nsteps', dest='nsteps', help='number of steps', type=int, required=True)
parser.add_argument('--device', dest='device', help='the GPU device', default='None')
parser.add_argument('--secdev', dest='secdev', help='the secondary GPU device', default='None')
parser.add_argument('--seed', dest='seed', help='the RNG seed', type=int, default=0)
parser.add_argument('--platform', dest='platform', help='the computation platform', default='CUDA')
args = parser.parse_args()

seed = int(1000*time.time()) % 16384 if args.seed == 0 else args.seed
print(f'Employed RNG seed is {seed}')

base = f'sinr-0p5fs'
platform_name = args.platform

dt = 0.5*unit.femtoseconds
temp = 298.15*unit.kelvin
rcut = 12*unit.angstroms
rswitch = 11*unit.angstroms
rcutIn = 8*unit.angstroms
rswitchIn = 5*unit.angstroms
tau = 10*unit.femtoseconds
gamma = 0.1/unit.femtoseconds
reportInterval = 180

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
nbforce.setUseDispersionCorrection(False)

integrator = openmm.LangevinIntegrator(temp, 1.0/unit.picoseconds, dt)

simulation = openmm.app.Simulation(pdb.topology, openmm_system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temp, seed)

computer = atomsmm.PressureComputer(openmm_system,
        pdb.topology,
#        openmm.Platform.getPlatformByName('CPU'),
        openmm.Platform.getPlatformByName('CUDA'),
        dict(Precision='mixed', DeviceIndex=args.secdev),
        temperature=temp)

dataReporter = atomsmm.ExtendedStateDataReporter(stdout, reportInterval, separator=',',
        step=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        atomicVirial=True,
        atomicPressure=True,
        nonbondedVirial=True,
        molecularVirial=True,
        molecularPressure=True,
        molecularKineticEnergy=True,
        coulombEnergy=True,
        pressure_computer=computer,
        speed=True,
        extraFile=f'{base}.csv')

pdbReporter = openmm.app.PDBReporter(f'{base}.pdb', 4*reportInterval)

simulation.reporters += [dataReporter, pdbReporter]
simulation.step(args.nsteps)
