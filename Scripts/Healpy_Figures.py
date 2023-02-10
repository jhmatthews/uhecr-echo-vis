import echo_util as util
import reverb_hp_subroutines as reverb
import warnings
import os
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy import constants as const
from astropy.coordinates import SkyCoord
from healpy.newvisufunc import projview, newprojplot
from tqdm import tqdm
plt.rcParams["figure.dpi"] = 300
print("Healpy version {}".format(hp.__version__))
util.set_mod_defaults()
util.set_times()
warnings.filterwarnings('ignore')


data_dir = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'Data'))
figure_dir = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'Figures'))
frame_dir = os.path.abspath(os.path.join(
    os.path.dirname(__file__), '..', 'Movies/Frames'))


def make_figs_9_10(root, ntimes=[34], folder=".", savefolder="."):

    print("Making composition figures for {}...".format(root))
    for ntime in tqdm(ntimes):
        #print("ntime:", ntime)
        directory = "{}/{}".format(folder, root)
        hpx_map, hpx_mapA, groups, lnA_bins, summed = reverb.get_maps_from_data2(
            ntime, nside=64, smooth=20, directory=directory, use_lnA_bins=True, summed=True)
        hpx_mapA_masked = hp.ma(np.log(hpx_mapA))

        yr = 365.25 * 24.0 * 3600.0
        delta_t = 0.03 * 1e6 * const.pc.cgs.value / const.c.cgs.value / yr / 1e6
        dt_plot = delta_t * util.dt_plot_sim

        #hpx_mapA_masked.mask = (hpx_map < 1e-2)
        to_plot = groups["7"]+groups["8"]
        to_plot = groups["0"]
        #print (groups.keys())
        override_plot_properties = {
            "figure_width": 9, "figure_size_ratio": 5./9.},
        fontsize_dict = {"xtick_label": 18,
                         "ytick_label": 18, "title": 18, "cbar_label": 18}

        all_groups = np.zeros_like(groups["0"])
        for i in groups.keys():
            all_groups += groups[i]

        to_plots = [groups["7"]+groups["8"], groups["2"]]

        #print (np.sum(summed), np.sum(hpx_map))

        cmaps = ["Oranges", "Purples"]
        labels = ["\ln A > 3.5", "1 < \ln A < 1.5"]
        savename = ["heavymap_{}_{}".format(root, ntime),
                    "lightmap_{}_{}".format(root, ntime)]

        fig = plt.figure()
        for i, to_plot in enumerate(to_plots):

            # plt.subplot(4,2,i+1)
            if np.max(to_plot) < 5e-7:
                vmin = 5e-7
                vmax = 1
            else:
                vmin = 0
                vmax = None

            projview(to_plot * util.skymap_renorm,
                     unit="UHECR Flux",
                     cb_orientation="vertical",
                     graticule=True,
                     projection_type="hammer",
                     graticule_labels=True,
                     xtick_label_color="k",
                     title=r"$t={:.1f}~{{\rm Myr}},~{}$".format(
                         ntime * dt_plot + 0.05, labels[i]),
                     cmap=cmaps[i], min=vmin, max=vmax,
                     override_plot_properties=override_plot_properties,
                     fontsize=fontsize_dict)

            plt.savefig("{}/{}.png".format(savefolder,
                        savename[i]), dpi=300, transparent=True)
            # plt.savefig("reverb_maps/{}.pdf".format(savename[i]), format="pdf")
            # plt.close("all")


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(
        prog='Healpy_Figures',
        description='make Healpy figures for Taylor et al. UHECR echoes')

    parser.add_argument('-m', '--movies', action='store_true',
                        default=False, help="Make movie frames")

    parser.add_argument('--no-figures', action='store_true', default=False,
                        help="Don't make figures")

    parser.add_argument('-f', '--files', nargs='+', type=str,
                        help="which files to make figures from",
                        default=["ModelA", "ModelB", "ModelC"])

    args = parser.parse_args()
    return args


def run_paper_figures(roots=["ModelA", "ModelB", "ModelC"]):
    times = util.plot_times

    yr = 365.25 * 24.0 * 3600.0
    delta_t = 0.03 * 1e6 * const.pc.cgs.value / const.c.cgs.value / yr / 1e6
    dt_plot = delta_t * util.dt_plot_sim
    maxes = [None for t in times]

    for iroot, root in enumerate(tqdm(roots)):

        make_figs_9_10(root, ntimes=[times[-1]], folder=data_dir,
                       savefolder=figure_dir)

        print ("Making skymap figures for {}...".format(root))
        reverb.make_skymap_plots(root, times, dt_plot, smooths=20,
                                 maxes=maxes, skip_header=1, folder=data_dir,
                                 savefolder=figure_dir, renorm = util.skymap_renorm)


def run_movie_frames(roots=["ModelA", "ModelB", "ModelC"]):
    times = np.arange(0, util.ndt_max, 1)

    yr = 365.25 * 24.0 * 3600.0
    delta_t = 0.03 * 1e6 * const.pc.cgs.value / const.c.cgs.value / yr / 1e6
    dt_plot = delta_t * util.dt_plot_sim
    maxes = [None for t in times]

    for iroot, root in enumerate(roots):

        make_figs_9_10(root, ntimes=times, folder=data_dir,
                       savefolder=frame_dir)

        reverb.make_skymap_plots(root, times, dt_plot, smooths=20,
                                 maxes=maxes, skip_header=1, folder=data_dir,
                                 savefolder=frame_dir, renorm = util.skymap_renorm)


if __name__ == "__main__":
    # parse arguments
    args = parse_args()

    # make movies and figures as required
    if args.no_figures == False:
        run_paper_figures(roots=args.files)
    if args.movies:
        run_movie_frames(roots=args.files)
