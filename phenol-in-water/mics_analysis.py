# -*- coding: utf-8 -*-
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mics as mx

from simtk import unit

T = 298.15*unit.kelvin
solute = 'phenol'
solvent = 'water'
kB = unit.BOLTZMANN_CONSTANT_kB*unit.AVOGADRO_CONSTANT_NA
states = pd.read_csv('../../alchemical_states.inp', sep='\s+', comment='#')
nstates = len(states.index)
states = states[states.weight != -np.inf].drop('weight', axis=1)
files = [f'{solute}-in-{solvent}_energy-{i:02d}.csv' for i in states.index]
renamer = lambda x: x.replace('Energy[','E').replace('] (kJ/mole)', '')
kT = (kB*T).value_in_unit(unit.kilojoules_per_mole)

mx.verbose = True
samples = mx.pooledsample()
for i, file in zip(states.index, files):
    print(f"Reading file {file}")
    df = pd.read_csv(file)
    df.drop(index=range(3000), inplace=True)
    df.rename(renamer, axis='columns', inplace=True)
    prev = max(i-1, 0)
    next = min(i+1, nstates-1)
    samples += mx.sample(df, f'beta*(E{i}-E0)', acfun=f'E{next}-E{prev}', beta=1/kT,
                         lambda_vdw=states.lambda_vdw[i],
                         lambda_coul=states.lambda_coul[i])

mixture = mx.mixture(samples, engine=mx.MICS(tol=1.0E-11))
# mixture = mx.mixture(samples.subsampling(), engine=mx.MBAR())
f = mixture.free_energies()
kT = (kB*T).value_in_unit(unit.kilocalories_per_mole)
#f['lambda[vdw+coul]'] = f['lambda_vdw'] + f['lambda_coul']
f['delta_g'] = kT*f['f']
f['error(delta_g)'] = kT*f['df']
print(f)
f.plot(x='lambda_vdw', y='delta_g', yerr='error(delta_g)')
mixture.histograms('potential', bins=200).plot(x='potential')
mixture.histograms('u0', bins=200).plot(x='u0')
# plt.figure()
# plt.imshow(mixture.Overlap)
plt.show()

cwd = os.path.split(os.getcwd())[1]
plt.savefig(f'../{cwd}.png')
f.to_csv(f'../{cwd}.csv')
