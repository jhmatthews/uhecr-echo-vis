import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
import os
import echo_util as util
import pandas as pd
from tqdm import tqdm
from constants import *

xsc = [-1.1,  3.2,  3.9,  3.7,  2.8,  2.3, -2.4, -3.2, -1.7]
ysc = [4.8,  3.8,  2.3, -0.3, -1.9, -2.4, -2.5,  2.8,  2.8]
res = 0.03
delta_t = 0.03 * 1e6 * PARSEC / C / YR / 1e6
dt_plot = delta_t * 10.0

CenAColor = "C6"
def draw_scatterers(radius=10.0):
	for i in range(len(xsc)):
		circle1 = plt.Circle((xsc[i],ysc[i]), 10.0 * res, fill = False, color = "k", ls="-")
		#print (i)
		plt.gca().add_artist(circle1)


	#plt.scatter(xsc, ysc)

	plt.scatter([0],[0], marker="+", color="C0", zorder=3, s=100)
	plt.scatter([-1.5],[3.2], marker="o", color = CenAColor, zorder=3, s=100)

def set_lims():
	lim = 5.2
	plt.xlim(-lim,lim)
	plt.ylim(-lim,lim)

def draw_ellipse(f1 = [0,0], f2 = [-1.5,3.2], time_delay = 31.3, color="C0", fill=True, lw=3, ls=(0,(8,1))):
	f1 = np.array(f1)
	f2 = np.array(f2)
	time_as_distance = (time_delay * 1e6 * YR * C) / PARSEC / 1e6
	a = time_as_distance / 2.0 # semi major axis
	#print (a)
	c = np.linalg.norm(f2-f1) / 2.0
	centre = 0.5 * (f2 - f1)
	diff = f2 - f1
	angle = np.arctan(diff[1]/diff[0])
	#print (angle)
	eccentricity = c / a 
	b = np.sqrt(1-(eccentricity*eccentricity)) * a 
	ellipse1 = matplotlib.patches.Ellipse((centre), 2 * a, 2 * b, np.degrees(angle), fill = fill, lw=lw, color=color, alpha=1.0, ls=ls)
	plt.gca().add_artist(ellipse1)


def run(cmap = "viridis", fmt="pdf"):
	shape = (601,601) # number of bins outputted from program
	util.set_mod_defaults()
	util.set_times()

	plt.figure(figsize=(6.35,5))
	n1 = 6
	dn = 6
	t0 = 11.7
	dt = 0.9784*1.0
	timestamps = np.arange(0,43,6.0 * dt)
	tcrits = [20.6,33.3]
	timestamps = np.sort(np.concatenate( (np.array(tcrits)-t0,timestamps)))
	mappable, colors = util.get_mappable(N=len(timestamps), vmin=t0+timestamps[0], 
		                                    vmax=t0+timestamps[-1], cmap_name=cmap)

	#tcrits = [20.6,33.3]
	times_to_use = timestamps[::-1]
	print ("Drawing ellipses...")
	for j in tqdm(range(len(times_to_use))):
		ntime = times_to_use[j]
		t_plot = ntime
		time_string = r"$t={:.1f}~{{\rm Myr}}$".format(t_plot)
		#print (j, time_string)
		xbin_summed = np.zeros(shape )
		#print 
		ls = (0,(8,1))
		#ls = "-"
		lw = 2
		for t in tcrits:
			if np.fabs(t_plot + t0 - t) < 0.1:
				ls = "-"
				lw = 4.5

		#if ntim
		draw_ellipse(f1 = [0,0], f2 = [-1.5,3.2], time_delay = t_plot + t0, color=colors[::-1][j], fill=False, ls=ls, lw = lw)


	draw_scatterers()
	#plt.colorbar()
	set_lims()

	plt.ylabel(r"$y~({\rm Mpc})$", labelpad=-2, fontsize=18)
	plt.xlabel(r"$x~({\rm Mpc})$", fontsize=18, labelpad=-2)
	cbar = plt.colorbar(mappable=mappable)
	cbar.set_label("Ballistic Arrival Time (Myr)")
			#plt.tight_layout()

	plt.text(0.2,-0.3,"$f_2$", color="C0", fontsize=18)
	plt.text(-1.2,3.1,r"$f_1$", color=CenAColor, fontsize=18)
	plt.subplots_adjust(right=0.96, left=0.11, bottom=0.11, top=0.98)
	plt.grid(ls=":")

	figure_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Figures'))
	plt.savefig("{}/fig5.{}".format(figure_dir, fmt), format=fmt)

if __name__ == "__main__":
	run()
