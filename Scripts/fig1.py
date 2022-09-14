import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib
from constants import *
import echo_util as util 
import os 
util.set_mod_defaults()
util.set_times()

xsc = [-1.1,  3.2,  3.9,  3.7,  2.8,  2.3, -2.4, -3.2, -1.7]
ysc = [4.8,  3.8,  2.3, -0.3, -1.9, -2.4, -2.5,  2.8,  2.8]
names = ["M83", "M64", "M94", "M81", "IC 342", "Maffei 1", "NGC 253", "Circinus", "NGC 4945"]
res = 0.03
delta_t = 0.03 * 1e6 * PARSEC / C / YR / 1e6
dt_plot = delta_t * 10.0

def draw_scatterers(radius=10.0):
    fontsize=15
    yoffs = [-0.2,0.2,0.2,0.2,0.2,-0.7,0.2,+0.4,-0.7]
    xoffs = [0.3,0.2,0.2,0.2,0.2,-0.7,0.2,-1.4,0]
    for i in range(len(xsc)):
        circle1 = plt.Circle((xsc[i],ysc[i]), 10.0 * res, fill = True, color = "k")
        #print (i)
        plt.gca().add_artist(circle1)
        plt.text(xsc[i]+xoffs[i], ysc[i]+yoffs[i], names[i], fontweight="bold", fontsize=fontsize)

    plt.scatter([0],[0], marker="+", color="C0", zorder=3, s=100)
    plt.scatter([-1.5],[3.2], color = "C6", zorder=3, s=100)
    #plt.text(xsc[i]+xoffs[i], ysc[i]+yoffs[i], names[i], fontname="Times New Roman", fontweight="bold", fontsize=16)
    plt.text(-1.4,3.3, "Cen A", color="C6", fontsize=fontsize)
    plt.text(-0.3,-0.6, "MW", color="C0", fontsize=fontsize)
    #plt.text(xsc[i]+xoffs[i], ysc[i]+yoffs[i], "MW", color="C0", fontsize=16)

    circle_centre = (0.362, 0.718)
    circle2 = plt.Circle(circle_centre, 3.746, fill = False, color = "k")
    #print (i)
    plt.gca().add_artist(circle2)
    plt.scatter(circle_centre[0], circle_centre[1], marker="x", color = "k", zorder=3, s=60)


def set_lims():
    lim = 5.2
    plt.xlim(-lim,lim)
    plt.ylim(-lim,lim)

def run():
    plt.figure(figsize=(5.62,5.48/0.998))
    draw_scatterers()
    set_lims()
    plt.ylabel(r"$y~({\rm Mpc})$", labelpad=-2, fontsize=18)
    plt.xlabel(r"$x~({\rm Mpc})$", fontsize=18, labelpad=-1)
    plt.subplots_adjust(right=0.96, left=0.11, bottom=0.11, top=0.98)
    plt.grid(ls=":")
    figure_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Figures'))
    plt.savefig("{}/fig1.pdf".format(figure_dir))

if __name__ == "__main__":
    run()
