import matplotlib.pyplot as plt
params = {
	'figure.subplot.hspace': 0.2,
	'figure.titlesize': 12,
	'font.size': 8,
	'text.usetex': True,
	'lines.linewidth': 0.5,
	'lines.markersize': 4,
	'errorbar.capsize' : 2,
	'legend.frameon': False,
	'savefig.format': 'png',
	'savefig.dpi': 600,
	'savefig.bbox': 'tight',
}
plt.rcParams.update(params)
