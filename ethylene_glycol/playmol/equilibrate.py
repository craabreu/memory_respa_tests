import argparse
import atomsmm
import numpy as np

from sys import stdout
from simtk import openmm
from simtk import unit
from simtk.openmm import app

parser = argparse.ArgumentParser()
parser.add_argument('-file', dest='file', required = True, help='Base name for PDB and XML files')
parser.add_argument('-steps', dest='nsteps', required = True, help='Number of MD steps')
parser.add_argument('-platform', dest='platform', required = True, choices=['CPU', 'CUDA', 'OpenCL'], help='Simulation platform')
args = parser.parse_args()
basename = args.file
nsteps = int(args.nsteps)

dt = 1*unit.femtoseconds
temp = 298.15*unit.kelvin
rcut = 12*unit.angstroms
rswitch = 11*unit.angstroms
rcutIn = 8*unit.angstroms
rswitchIn = 5*unit.angstroms
seed = 98745
platform_name = args.platform
reportInterval = 100

platform = openmm.Platform.getPlatformByName(platform_name)
properties = dict(Precision='mixed') if platform_name == 'CUDA' else dict()

pdb = app.PDBFile('{}_raw.pdb'.format(basename))

residues = [atom.residue.name for atom in pdb.topology.atoms()]
solute_atoms = set(i for (i, name) in enumerate(residues) if name == 'aaa')
solvent_atoms = set(range(pdb.topology.getNumAtoms())) - solute_atoms

forcefield = app.ForceField('{}.xml'.format(basename))

openmm_system = forcefield.createSystem(pdb.topology,
                                        nonbondedMethod=openmm.app.PME,
                                        nonbondedCutoff=rcut,
                                        rigidWater=False,
                                        constraints=None,
                                        removeCMMotion=False)
nbforce = openmm_system.getForce(atomsmm.findNonbondedForce(openmm_system))
nbforce.setUseSwitchingFunction(True)
nbforce.setSwitchingDistance(rswitch)
nbforce.setUseDispersionCorrection(True)
openmm_system.addForce(openmm.MonteCarloBarostat(1*unit.atmospheres, temp, 25))

integrator = openmm.LangevinIntegrator(temp, 1.0/unit.picoseconds, 1.0*unit.femtosecond)
simulation = openmm.app.Simulation(pdb.topology, openmm_system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temp, seed)

reporter = atomsmm.ExtendedStateDataReporter(stdout, reportInterval, separator=',', step=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True,
    volume=True, density=True, speed=True, extraFile='properties.csv')
simulation.reporters.append(reporter)
simulation.step(nsteps)

state = simulation.context.getState(getPositions=True)
coords = state.getPositions(asNumpy=True).value_in_unit(unit.angstroms)
volume = state.getPeriodicBoxVolume().value_in_unit(unit.angstroms**3)
Lbox = [volume**(1/3)]*3

out = open('box.temp', 'w')
print('box lengths {} {} {}'.format(*Lbox), file=out)
print('reset xyz', file=out)
print('build', file=out)
print(pdb.topology.getNumAtoms(), file=out)
for atom, coord in zip(pdb.topology.atoms(), coords):
    print(atom.name, *coord, file=out)
out.close()
