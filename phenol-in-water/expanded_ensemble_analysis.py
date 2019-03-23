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
parser.add_argument('--nrows', dest='nrows', help='the number of rows to read', type=int, default=None)
parser.add_argument('--interval', dest='reportInterval', help='the report interval', type=int, default=90)
args = parser.parse_args()

solute = 'phenol'
solvent = 'water'
#base = '{}-in-{}_energy'.format(solute, solvent)
base = '{}-in-{}_energy'.format(solute, solvent)
temp = 298.15*unit.kelvin
states_data = pd.read_csv('expanded_ensemble_states.inp', sep='\s+', comment='#')
expandedEnsembleReporter = atomsmm.reporters.ExpandedEnsembleReporter(stdout,
    args.reportInterval, separator=',', states=states_data, temperature=temp)
expandedEnsembleReporter.read_csv(f'{base}.csv', nrows=args.nrows)
expandedEnsembleReporter.state_sampling_analysis(staging_variable='lambda_vdw')
expandedEnsembleReporter.walking_time_analysis(history=True)
