import numpy as np 
import matplotlib.pyplot as plt 
import echo_util as util
from constants import *
import os 
import astropy.coordinates as coords
import astropy.units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord

def flip_data2(l):
	if np.isscalar(l):
		l = np.array([l])
	l[l<np.pi] *= -1
	l[l>np.pi] = (2.0*np.pi) - l[l>np.pi]

	if len(l) == 1:
		return l[0]
	else:
		return l

def draw_planes(colors=["k", "C1"]):
	labels = ["Galactic", "Supergalactic"]
	frames = ["galactic", "supergalactic"]
	ls = [":", ":"]
	x = np.arange(0,360.1,0.1)
	y = np.zeros_like(x)

	for i, f in enumerate(frames[1:]):
		coords = SkyCoord(x*u.degree, y*u.degree, frame=f)
		ldata = flip_data2(coords.galactic.l.radian)
		args = np.argsort(ldata)
		plt.plot( ldata[args][:-1], coords.galactic.b.radian[args][:-1], c=colors[i], label=labels[i+1], lw=2, ls=ls[i])

def set_labels():
	ticks = np.array([-180,-135,-90,-45,0,45,90,135,180])
	ticks = np.array([-180,-120,-60,0,60,120,180])

	#ticks = np.array([180,120,60,0,360,300,240])
	ticks_rad = np.radians(np.array(ticks))
	plt.gca().set_xticks(ticks_rad)

	#yticks = [-89.999,-75,-60,-45,-30,-15,0,15,30,45,60,75,89.999]
	#plt.gca().set_yticks(np.radians(np.array(yticks)))
	plt.gca().set_longitude_grid_ends(90)
	yticks = np.linspace(-np.pi/2,np.pi/2,7)
	plt.gca().set_yticks(yticks)
	plt.gca().set_yticklabels([r"${}^\circ$".format(i) for i in [-90,-60,-30,0,30,60,90]], fontsize=18)

	tick_labels = ticks
	tick_labels[tick_labels>0] = 360-tick_labels[tick_labels>0]
	tick_labels[tick_labels<0] *= -1
	tick_labels = [r"${}^\circ$".format(i) for i in tick_labels]
	plt.gca().set_xticklabels(tick_labels, fontsize=18)


def get_galaxies():
	coords = dict()
	coords["M83"] = SkyCoord("13 37 00.919 -29 51 56.74", frame='icrs', unit=(u.hourangle, u.deg))
	coords["CenA"] = SkyCoord("13 25 27.61509104 -43 01 08.8056025 ", frame='icrs', unit=(u.hourangle, u.deg))
	coords["M82"] = SkyCoord("09 55 52.430 +69 40 46.93", frame='icrs', unit=(u.hourangle, u.deg))
	coords["NGC253"] = SkyCoord("00 47 33.134 -25 17 19.68", frame='icrs', unit=(u.hourangle, u.deg))
	coords["IC342"] = SkyCoord("03 46 48.514 +68 05 45.98", frame='icrs', unit=(u.hourangle, u.deg))
	coords["Circinus"] = SkyCoord("14 13 09.906 -65 20 20.47", frame='icrs', unit=(u.hourangle, u.deg))
	coords["NGC4945"] = SkyCoord("13 05 27.279 -49 28 04.44", frame='icrs', unit=(u.hourangle, u.deg))
	coords["M94"] = SkyCoord("12 50 53.0737971432 +41 07 12.900884628", frame='icrs', unit=(u.hourangle, u.deg))
	coords["M64"] = SkyCoord("12 56 43.696 +21 40 57.57", frame='icrs', unit=(u.hourangle, u.deg))
	coords["Maffei1"] = SkyCoord("02 36 35.4662282832 +59 39 17.507008584", frame='icrs', unit=(u.hourangle, u.deg))

	labels = dict()
	labels["M83"] = "M83"
	labels["CenA"] = "Cen A"
	labels["M82"] = r"M81 \& M82"
	labels["NGC253"] = "NGC 253"
	labels["IC342"] = "IC 342"
	labels["Circinus"] = "Circinus"
	labels["NGC4945"] = "NGC 4945"
	labels["M94"]="M94"
	labels["M64"]="M64"
	labels["Maffei1"] = r"Maffei 1 \& 2"
	return (coords, labels)

def init_figure(subplot, projection="aitoff"):
	fig = plt.figure(figsize=(9, 5))
	ax = fig.add_subplot(subplot, projection=projection)
	ax.grid(True)

	return (fig, ax)



def add_Galaxies(sbg_color="k", cena_color="C6", text=True):
	i = 0
	j = 0

	# offsets for labels of galaxies (done by hand)
	offsets = dict()
	offsets["M83"] = (20,5)
	offsets["CenA"] = (-60,0)
	offsets["M82"] = (-75,5)
	offsets["NGC253"] = (30,13)
	offsets["IC342"] = (45,5)
	offsets["Circinus"] = (15,-30)
	offsets["NGC4945"] = (30,10)
	offsets["M94"] = (20,-9)
	offsets["M64"] = (-30,-12)
	offsets["Maffei1"] = (0,-7)

	# get the coordinates and labels for each galaxy
	coords, labels = get_galaxies()

	fs = 20
	ms = 250
	doff = [1,0.05,0.05,1.0,0.08,0.1,0.08,1,1,1] # fraction by which to offset point of the arrow

	for keys, values in coords.items():
		ra_rad = coords[keys].galactic.l.radian

		ra_rad = flip_data2(ra_rad)
		ycoord_to_use = coords[keys].galactic.b.radian
		if keys == "NGC4945":
			ha = "left"
			va = "bottom"
		else:
			ha = "center"
			va = "center"

		if keys == "CenA":
			plt.scatter(ra_rad, ycoord_to_use, s=ms, c=cena_color, marker="o", label=None)
		else:
			if i == 0:
				label ="SBGs"
			else:
				label = None
			plt.scatter(ra_rad, ycoord_to_use, s=ms, c=sbg_color, marker="o", label=None)
			i+=1

		if text:
			if keys == "CenA":
				text_color = cena_color
			else:
				text_color = sbg_color

			xytext = (ra_rad + np.radians(offsets[keys][0]), ycoord_to_use + np.radians(offsets[keys][1]))
			xypoint = (ra_rad + np.radians(doff[j]*offsets[keys][0]), ycoord_to_use+ np.radians(doff[j]*offsets[keys][1]))
			# if arrow[j]:
			plt.annotate(labels[keys], 
						 xypoint, 
						 xytext=xytext,
						 arrowprops=dict(arrowstyle="->",connectionstyle="arc3", color=text_color), 
						 color=text_color, fontsize=fs, fontweight='bold',
						 horizontalalignment=ha, verticalalignment=va)
		j+=1


def run():
	print ("Making figure 2...")

	util.set_mod_defaults()
	util.set_times()

	init_figure(111, projection="hammer")
	draw_planes()
	add_Galaxies()
	set_labels()

	plt.subplots_adjust(top=0.92, bottom=0.08, right=0.96,left=0.04)
	figure_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Figures'))
	plt.savefig("{}/fig2.pdf".format(figure_dir))
	print ("Done.")

if __name__ == "__main__":
	run()


