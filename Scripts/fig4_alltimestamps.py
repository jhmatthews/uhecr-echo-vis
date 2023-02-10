import numpy as np
import matplotlib.pyplot as plt
import echo_util as util
from tqdm import tqdm
import os
import fig4
import Make_Figures

def make_plots(t, flux, labels=["A", "B", "C"], colors=["C0", "C2", "C1"]):
    ''' 
    Make figure 4 from Taylor et al. CoG echoes paper
    '''

    # plot each model output
    runs = ["ModelA", "ModelB", "ModelC","ModelB_smooth2","ModelB_no4945"]
    for imod, folder in enumerate(runs):
        print ("Making plot for Model A...")
        timestamps=np.arange(0,util.ndt_max,1)
        for ntime in tqdm(timestamps):
            plt.figure(figsize=(6, 4))
            n = flux[imod]
            n0 = flux[0]
            if ntime > 0:
                plt.plot(t[:ntime+1], n[:ntime+1] / np.max(n0), lw=3,
                        label="Model {}".format(labels[imod]), c="C0")
                if n[ntime] / np.max(n0) < 3e-5:
                    nplot = 3e-5
                else:
                    nplot = n[ntime] / np.max(n0)
                plt.scatter(t[ntime], nplot, c="C0")
            #else:
            plt.plot(t, n / np.max(n0), lw=3,
                    label="Model {}".format(labels[imod]), c="C0", alpha=0.3)

            # cosmetics
            plt.ylim(3e-5, 2)
            plt.xlim(0, 42)
            ax = plt.gca()
            for t0 in [11.7, 20.6, 33.3]:
                ax.axvline(t0, lw=2, ls="--", alpha=0.8, color="k")
                trans = ax.get_xaxis_transform()
                time_string = r"${:.1f}~{{\rm Myr}}$".format(t0)
                plt.text(t0+0.01, .05, time_string, transform=trans,
                         rotation=270, fontsize=14)

            plt.xlabel(r"$t~({\rm Myr})$", fontsize=18, labelpad=-1)
            plt.ylabel("Local UHECR Density (Arb.)", fontsize=16)
            plt.subplots_adjust(right=0.98, top=0.98, left=0.13, bottom=0.14)
            plt.semilogy()
            figure_dir = os.path.abspath(os.path.join(os.path.dirname(__file__ ), '..', 'Movies/Frames'))
                    #print (jm_util.get_aspect(plt.gcf().get_axes()[0]))
            plt.savefig(
                "{}/timeseries_{}_{:03d}.png".format(figure_dir, folder, ntime), format="png", dpi=300)


def run(load=True, reset_cycler = True):
    util.set_mod_defaults()
    util.set_times()
    util.set_cycler("colorblind")

    t, flux = fig4.prepare_data(load=load)
    make_plots(t, flux)

    if reset_cycler:
        # revert to default colors
        util.set_cycler('default')

if __name__ == "__main__":
    args = Make_Figures.parse_args()
    load = (args.raw == False)
    run(load = load)

