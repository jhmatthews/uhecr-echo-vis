import numpy as np
import matplotlib.pyplot as plt
import echo_util as util
import pandas as pd
from constants import *
from tqdm import tqdm
import h5py, os
util.set_mod_defaults()
util.set_times()
import fig3
import argparse
import Make_Figures

res = 0.03
delta_t = 0.03 * 1e6 * PARSEC / C / YR / 1e6
dt_plot = delta_t * util.dt_plot_sim



def figure_from_hdf5(hf, timestamps, Models, 
    model_name = ["Model A (Pulse)", "Model B (Decline)",
                  "Model C (Delayed Escape)"]):
    """
    Make spatial density figures from and hfd5 database
    """

    # labels and save names
    # model_name = ["MW Halo", "Model A (Pulse)", "Model B (Decline)",
    #               "Model C (Delayed Escape)"]
    figname = Models

    for irun, Model in enumerate(Models):

        
        # get data for this specific model
        group = hf.get(Model)

        # loop over each timestamp and plot results
        print("Plotting for {}".format(Model))

        for j in tqdm(range(len(timestamps))):
            fig, ax = plt.subplots(figsize=(5, 4), nrows=1, ncols=1)
            ntime = timestamps[j]
            t_plot = dt_plot * ntime
            time_string = r"$t={:.1f}~{{\rm Myr}}$".format(t_plot)

            # get data from HDF5 group
            x = group.get("x")
            y = group.get("y")
            xbin = group.get("cr{:03d}".format(ntime))

            # plot and adjust cosmetics
        
            mappable1 = ax.pcolormesh(x, y, np.log10(
                xbin), cmap="Blues", vmin=-4.99, vmax=0)

            fig3.draw_scatterers(ax)

            # limits, labels  and ticks 
            ax.set_xlim(-9, 9)
            ax.set_ylim(-9, 9)
            ax.grid(ls=":", alpha=0.5)
            ax.set_xticks([-5, 0, 5])
            ax.set_yticks([-5, 0, 5])
      
            ax.set_ylabel(r"$y~({\rm Mpc})$", labelpad=-2, fontsize=18)
            ax.text(-7, 7, model_name[irun], fontsize=14)

            ax.set_xlabel(r"$x~({\rm Mpc})$", labelpad=-2, fontsize=18)
            ax.text(-7, -8, time_string)

            # adjust subplots and add colorbar
            top = 0.98
            bottom = 0.13
            right = 0.8
            left = 0.12
            wspace = 0.05
            fig_height = ((top-bottom))
            ax1 = plt.gcf().add_axes([right + 0.01, bottom, 0.05, fig_height])
            cbar1 = plt.colorbar(cax=ax1, mappable=mappable1,
                                 orientation="vertical")
            cbar1.set_label(r"$\log_{10}[{\rm UHECR~Density~(Arb.)}]$")

            plt.subplots_adjust(hspace=wspace, wspace=wspace,
                                right=right, left=left, bottom=bottom, top=top)
            
            figure_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Movies/Frames'))
            #print (jm_util.get_aspect(plt.gcf().get_axes()[0]))
            plt.savefig(
                "{}/xy_{}_{:03d}.png".format(figure_dir, figname[irun], ntime), format="png", dpi=300)


def run(timestamps=np.arange(0,util.ndt_max,1), load=True, Models=["ModelA", "ModelB", "ModelC"]):

    # data for Models A,B,C should be stored in this filename
    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Data'))
    hf_fname = "{}/xy_data_all.h5".format(data_dir)

    # create database from raw cr-reverb output
    if load == False:
        fig3.create_hdf5_dataset(Models, timestamps, hf_fname=hf_fname)

    print("Loading hdf5 dataset...")
    hf = h5py.File(hf_fname, "r")
    figure_from_hdf5(hf, timestamps, Models)
    hf.close()

    print("All done.")


if __name__ == "__main__":
    args = Make_Figures.parse_args()
    load = (args.raw == False)
    run(load = load)
