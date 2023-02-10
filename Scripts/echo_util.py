import matplotlib.pyplot as plt 
import numpy as np
import matplotlib 
from cycler import cycler

# output file frequency 
dt_plot_sim = 5.0

# Arbitrary renormalisation for skymaps to make colourmaps have less 000
skymap_renorm = 1000.0

ndt_max = 87
plot_times = [24, 42, 64, 68]

def set_mod_defaults():
	'''
	set some default plot parameters
	uses Computer Modern font 
	(use set_times function to use Times instead)
	'''
	plt.rcParams["text.usetex"] = "True"
	plt.rcParams['font.serif']=['cm']
	plt.rcParams['font.family']='serif'	
	plt.rcParams['text.latex.preamble']=r'\usepackage{amsmath}'
	plt.rcParams['font.size']=18
	plt.rcParams['xtick.labelsize']=15
	plt.rcParams['ytick.labelsize']=15
	plt.rcParams['legend.fontsize']=14
	plt.rcParams['axes.titlesize']=16
	plt.rcParams['axes.labelsize']=16
	plt.rcParams['axes.linewidth']=2
	plt.rcParams["lines.linewidth"] = 2.2
	plt.rcParams['xtick.top']='True'
	plt.rcParams['xtick.bottom']='True'
	plt.rcParams['xtick.minor.visible']='True'
	plt.rcParams['xtick.direction']='out'
	plt.rcParams['ytick.left']='True'
	plt.rcParams['ytick.right']='True'
	plt.rcParams['ytick.minor.visible']='True'
	plt.rcParams['ytick.direction']='out'
	plt.rcParams['xtick.major.width']=1.5
	plt.rcParams['xtick.minor.width']=1
	plt.rcParams['xtick.major.size']=4
	plt.rcParams['xtick.minor.size']=3
	plt.rcParams['ytick.major.width']=1.5
	plt.rcParams['ytick.minor.width']=1
	plt.rcParams['ytick.major.size']=4
	plt.rcParams['ytick.minor.size']=3

def get_mappable(N, vmin=0, vmax=1, cmap_name = "Spectral", return_func = False):
	my_cmap = matplotlib.cm.get_cmap(cmap_name)
	colors = my_cmap(np.linspace(0,1,num=N))

	norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
	mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_name)
	if return_func:
		fcol = colour_func(norm, cmap_name)
		return (mappable, colors, mappable.to_rgba)
	else:
		return (mappable, colors)

def set_times():
	'''
	set fonts to times.
	'''
	from matplotlib import rc
	rc('font',**{'family':'serif','serif':['Times']})
	rc('text', usetex=True)


def set_cycler(cmap_name = "viridis", N = None):
	'''
	set the cycler to use a colormap
	'''
	if cmap_name == "default" or (N is None and cmap_name != "colorblind"):
		my_cycler = cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
	elif cmap_name == "colorblind":
		my_cycler = cycler('color', ['#0072B2', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9'])
	else:
		_, colors = get_mappable(N, cmap_name = cmap_name)
		#if type(style) == str:
		my_cycler = (cycler(color=colors)) 

	plt.rc('axes', prop_cycle=my_cycler)


def read_skymap(ntime, directory = "./", skip_header=3):
    '''
    get skymap data run from cr-reverb and stored in a directory
    '''
    fname = "skymap_{:08d}.out".format(ntime)
    data = np.genfromtxt("{}/output_000/{}".format(directory, fname), unpack=True, skip_header=skip_header)
    return (data)
