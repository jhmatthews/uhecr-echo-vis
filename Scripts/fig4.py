import numpy as np
import matplotlib.pyplot as plt
from constants import *
import sys
import echo_util as util
from tqdm import tqdm
import os


def total_flux_from_skymaps(directory, timestamps):
    '''
    read skymap data at each timestamp and sum the data
    '''
    # array to store fluxes
    flux = np.zeros(len(timestamps), dtype=float)

    # read each skymap and store in array
    for j in tqdm(range(len(timestamps))):
        ntime = timestamps[j]
        data = util.read_skymap(int(ntime), directory=directory)
        flux[j] = np.sum(data[2])

    return (flux)


def save_all_fluxes(savename, t, fluxA, fluxB, fluxC):
    to_save = np.column_stack((t, fluxA, fluxB, fluxC))
    np.savetxt(savename, to_save)


def prepare_data(load=True):

    delta_t = 0.03 * 1e6 * PARSEC / C / YR / 1e6
    dt_plot = delta_t * 10.0

    shape = (601, 601)  # number of bins outputted from program

    timestamps = np.arange(0, 45, 1.0)
    runs = ["ModelA", "ModelB", "ModelC"]
    flux = []

    header = "t (Myr)"
    for run in runs:
        header+="\t{}".format(run)

    data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Data'))
    savename = "{}/local_flux.dat".format(data_dir)
    if load:
        data = np.genfromtxt(savename, unpack=True)
        t = data[0, :]
        flux = data[1:, :]
    else:
        print("Recalculating local flux from outputs.")
        print("Warning, reading this data can take a while...")
        for imod, run in enumerate(runs):
            print(run)
            directory = "{}/{}".format(data_dir, run)
            t = timestamps * dt_plot

            flux.append(total_flux_from_skymaps(directory, timestamps))

        save_all_fluxes(savename, t, flux[0], flux[1], flux[2])

    return (t, flux)


def make_plots(t, flux, labels=["A", "B", "C"], colors=["C0", "C2", "C1"]):
    ''' 
    Make figure 4 from Taylor et al. CoG echoes paper
    '''
    plt.figure(figsize=(6, 4))

    # plot each model output
    for imod, folder in enumerate(flux):
        n = flux[imod]
        n0 = flux[0]
        plt.plot(t, n / np.max(n0), lw=3,
                 label="Model {}".format(labels[imod]), c=colors[imod])

    # plot exponential with 3 Myr decay
    plt.plot(t, np.exp(-t/3.0), c="C3", ls="-.",
             lw=2.5, label=r"$e^{-t/(3~{\rm Myr})}$")

    # cosmetics
    plt.legend(frameon=True)
    plt.ylim(3e-5, 2)
    plt.xlim(0, 42)
    ax = plt.gca()
    for t in [11.7, 20.6, 33.3]:
        ax.axvline(t, lw=2, ls="--", alpha=0.8, color="k")
        trans = ax.get_xaxis_transform()
        time_string = r"${:.1f}~{{\rm Myr}}$".format(t)
        plt.text(t+0.01, .05, time_string, transform=trans,
                 rotation=270, fontsize=14)

    plt.xlabel(r"$t~({\rm Myr})$", fontsize=18, labelpad=-1)
    plt.ylabel("Local UHECR Density (Arb.)", fontsize=16)
    plt.subplots_adjust(right=0.98, top=0.98, left=0.13, bottom=0.14)
    plt.semilogy()
    figure_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Figures'))
    plt.savefig("{}/fig4.pdf".format(figure_dir))

def run(load=True, reset_cycler = True):
    util.set_mod_defaults()
    util.set_times()
    util.set_cycler("colorblind")

    t, flux = prepare_data(load=load)
    make_plots(t, flux)

    if reset_cycler:
        # revert to default colors
        util.set_cycler('default')

if __name__ == "__main__":
    run(load = True)

