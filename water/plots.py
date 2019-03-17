import itertools
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import figstyle
import os
from scipy.fftpack import fft

dt_scaling = 'log'
methods = dict(sinr='SIN(R)', langevin='Memory Langevin')
tools = methods.keys()
step_size = {'0.5': 0.5}
label = {'0.5': '0.5 fs (*)'}
for case in ['01', '03', '06', '09', '15', '30', '45', '90']:
    step_size[case] = float(case)
    label[case] = '{} fs'.format(int(case))

def data(type, timesteps):
    repo = {tool:[] for tool in tools}
    for dt, tool in itertools.product(timesteps, tools):
        file = '{}/results/dt{}fs_{}.csv'.format(tool, dt.replace('.', 'p'), type)
        if os.path.isfile(file):
            repo[tool].append(pd.read_csv(file, skipinitialspace=True))
    return repo

# Intermolecular radial distribution functions:
def plot_radial_distribution_functions(timesteps):
    rdf = data('rdf', timesteps)
    for tool, method in methods.items():
        fig, ax = plt.subplots(3, 1, figsize=(3.37,5.6), sharex=True)
        fig.suptitle(f'{method}: Radial distribution functions')
        fig.subplots_adjust(hspace=0.1)
        ax[-1].set_xlabel('Distance (\\AA)')
        for dt, gr in zip(timesteps, rdf[tool]):
            distance = gr['Distance [pm]']/100
            for i, pair in enumerate(['O-O', 'O-H', 'H-H']):
                ax[i].plot(distance, gr[f'g({pair})'], label=label[dt])
                ax[i].set_ylabel(f'g({pair})')
        ax[2].legend(loc='lower right', ncol=2)
        fig.savefig(f'{tool}_rdf.png')

# Bond length distributions:
def plot_bond_length_distributions(timesteps):
    bond = data('bond', timesteps)
    for tool, method in methods.items():
        fig, ax = plt.subplots(1, 1, figsize=(3.37,2.3), sharex=True)
        fig.suptitle(f'{method}: Bond length distributions')
        ax.set_xlabel('Distance (\\AA)')
        for dt, gr in zip(timesteps, bond[tool]):
            distance = gr['# Distance [pm]']/100
            ax.plot(distance, gr['g(r)'], label=label[dt])
            ax.set_ylabel(f'frequency')
        ax.legend()
        fig.savefig(f'{tool}_bond.png')

# Angle distributions:
def plot_angle_distributions(timesteps):
    angle = data('angle', timesteps)
    for tool, method in methods.items():
        fig, ax = plt.subplots(1, 1, figsize=(3.37,2.3), sharex=True)
        fig.suptitle(f'{method}: Angle distributions')
        ax.set_xlabel('Angle (\\textdegree)')
        for dt, gr in zip(timesteps, angle[tool]):
            distance = gr['# Angle (degree)']
            ax.plot(distance, gr['Occurrence'], label=label[dt])
            ax.set_ylabel(f'frequency')
        ax.legend()
        fig.savefig(f'{tool}_angle.png')

# Average bonds and angles:
def plot_bond_and_angle_averages():
    fig, ax = plt.subplots(2, 1, figsize=(3.37,4.6), sharex=True)
    fig.suptitle('Average bond lengths and angles')
    ax[-1].set_xlabel('Time step size (fs)')
    for i, type in enumerate(['bond', 'angle']):
        for tool, method in methods.items():
            df = pd.read_csv(f'{tool}/results/{type}_stats.csv')
            dt, mean = df['dt'], df['mean']
            if type == 'bond':
                mean /= 100
            ax[i].plot(dt, mean, marker='o', label=method)
            ax[i].set_xscale(dt_scaling)
    ax[0].set_ylabel('Bond length (\\AA)')
    ax[0].legend(loc='best', title='$r_0 = 1.012$ \\AA')
    ax[1].set_ylabel('Angle (\\textdegree)')
    ax[1].legend(loc='best', title='$\\theta_0 = 113.24$\\textdegree')
    fig.savefig('average_bonds_and_angles.png')

# Thermo Properties:
def plot_properties():
    n = 3
    fig, ax = plt.subplots(n, 1, figsize=(3.37,n*2+0.3), sharex=True)
    fig.suptitle('Average Properties')
    ax[-1].set_xlabel('Time step size (fs)')
    properties = dict(PotEng='Potential Energy (kJ/mol)',
                      Press='Atomic Pressure (atm)',
                      MolPress='Molecular Pressure (atm)')
    for tool, method in methods.items():
        df = pd.read_csv(f'{tool}/results/properties.csv')
        for i, type in enumerate(properties.keys()):
            dt, mean = df['dt'], df[type]
            ax[i].plot(dt, mean, marker='o', label=method)
            ax[i].legend()
            ax[i].set_xscale(dt_scaling)
            ax[i].set_ylabel(properties[type])
    fig.savefig('average_properties.png')

all = ['0.5', '01', '03', '06', '09', '15', '30', '45', '90']
plot_radial_distribution_functions(all)
plot_bond_length_distributions(all)
plot_angle_distributions(all)
plot_bond_and_angle_averages()
plot_properties()
plt.show()
