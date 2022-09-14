import numpy as np
import matplotlib.pyplot as plt
import echo_util as util
import pandas as pd
from constants import *
from tqdm import tqdm
import h5py, os
util.set_mod_defaults()
util.set_times()

res = 0.03
delta_t = 0.03 * 1e6 * PARSEC / C / YR / 1e6
dt_plot = delta_t * 10.0


def draw_scatterers(ax, radius=10.0):
    xsc = [-1.1,  3.2,  3.9,  3.7,  2.8,  2.3, -2.4, -3.2, -1.7]
    ysc = [4.8,  3.8,  2.3, -0.3, -1.9, -2.4, -2.5,  2.8,  2.8]
    for i in range(len(xsc)):
        circle1 = plt.Circle((xsc[i], ysc[i]), 10.0 *
                             res, fill=False, color="k")
        #print (i)
        ax.add_artist(circle1)
    ax.scatter([0], [0], marker="+", color="k")
    ax.scatter([-1.5], [3.2], color="C6")


def create_hdf5_dataset(Models, timestamps, hf_fname="Data/xy_data.h5"):
    """
    Create hfd5 database from cr-reverb raw outputs
    """

    shape = (601, 601)  # number of bins outputted from program
    compression_kwargs = {"compression": "gzip", "compression_opts": 9}

    hf = h5py.File(hf_fname, 'w')
    print("Creating hdf5 dataset...")
    for irun, Model in enumerate(Models):
        print("Reading data for {}".format(Model))

        g = hf.create_group(Model)

        for j in tqdm(range(len(timestamps))):
            ntime = timestamps[j]

            data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Data'))
            directory = "{}/{}/output_000/".format(data_dir, Model)
            fname = "xbin_{:08d}.out".format(ntime)
            x, y, xbin_xyt = np.genfromtxt(
                "{}/{}".format(directory, fname), usecols=np.arange(0, 3), unpack=True, skip_header=3)

            x = np.reshape(x, (shape))
            y = np.reshape(y, (shape))
            xbin = np.reshape(xbin_xyt, (shape))

            # store data in hdf5
            if j == 0:
                g.create_dataset('x', data=x, **compression_kwargs)
                g.create_dataset('y', data=y, **compression_kwargs)
            g.create_dataset("cr{:03d}".format(ntime),
                             data=xbin, **compression_kwargs)

    hf.close()


def figure_from_hdf5(hf, timestamps, Models):
    """
    Make spatial density figures from and hfd5 database
    """

    # labels and save names
    model_name = ["Model A (Pulse)", "Model B (Decline)",
                  "Model C (Delayed Escape)"]
    figname = ["fig3", "figb1", "figb2"]

    for irun, Model in enumerate(Models):

        fig, axes = plt.subplots(figsize=(15.84, 4), nrows=1, ncols=4)
        
        # get data for this specific model
        group = hf.get(Model)

        # loop over each timestamp and plot results
        print("Plotting for {}".format(Model))
        for j in tqdm(range(len(timestamps))):
            ntime = timestamps[j]
            t_plot = dt_plot * ntime
            time_string = r"$t={:.1f}~{{\rm Myr}}$".format(t_plot)

            # get data from HDF5 group
            x = group.get("x")
            y = group.get("y")
            xbin = group.get("cr{:03d}".format(ntime))

            # plot and adjust cosmetics
            ax = axes[j]
            mappable1 = ax.pcolormesh(x, y, np.log10(
                xbin), cmap="Blues", vmin=-2.99, vmax=2)

            draw_scatterers(ax)

            # limits, labels  and ticks 
            ax.set_xlim(-9, 9)
            ax.set_ylim(-9, 9)
            ax.grid(ls=":", alpha=0.5)
            ax.set_xticks([-5, 0, 5])
            ax.set_yticks([-5, 0, 5])
            if j > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel(r"$y~({\rm Mpc})$", labelpad=-2, fontsize=18)
                ax.text(-7, 7, model_name[irun], fontsize=14)

            ax.set_xlabel(r"$x~({\rm Mpc})$", labelpad=-2, fontsize=18)
            ax.text(-7, -8, time_string)

        # adjust subplots and add colorbar
        top = 0.98
        bottom = 0.13
        right = 0.94
        left = 0.049
        wspace = 0.05
        fig_height = ((top-bottom))
        ax1 = plt.gcf().add_axes([right + 0.005, bottom, 0.012, fig_height])
        cbar1 = plt.colorbar(cax=ax1, mappable=mappable1,
                             orientation="vertical")
        cbar1.set_label(r"$\log_{10}[{\rm UHECR~Density~(Arb.)}]$")

        plt.subplots_adjust(hspace=wspace, wspace=wspace,
                            right=right, left=left, bottom=bottom, top=top)
        
        figure_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Figures'))
        plt.savefig(
            "{}/{}.png".format(figure_dir, figname[irun]), format="png", dpi=300)


def run(timestamps=[4, 12, 21, 34], load=True, Models=["ModelA", "ModelB", "ModelC"]):

    # data for Models A,B,C should be stored in this filename
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Data'))
    hf_fname = "{}/xy_data.h5".format(data_dir)

    # create database from raw cr-reverb output
    if load == False:
        create_hdf5_dataset(Models, timestamps, hf_fname=hf_fname)

    print("Loading hdf5 dataset...")
    hf = h5py.File(hf_fname, "r")
    figure_from_hdf5(hf, timestamps, Models)
    hf.close()

    print("All done.")


if __name__ == "__main__":
    run(load = True)
