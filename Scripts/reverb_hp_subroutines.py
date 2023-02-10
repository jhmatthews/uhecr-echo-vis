import numpy as np
import matplotlib.pyplot as plt
#import jm_util 
import astropy.coordinates as coords
import astropy.units as u
import healpy as hp
from astropy import constants as const
from astropy.coordinates import SkyCoord
from healpy.newvisufunc import projview, newprojplot
from tqdm import tqdm

def get_data(ntime, directory = "./", skip_header=0):
    '''
    get skymap data run from cr-reverb and stored in a directory
    '''    
    fname = "skymap_{:08d}.out".format(ntime)
    l,b, skyt,skyA, H, He, N, Fe = np.genfromtxt("{}/output_000/{}".format(directory, fname), usecols=np.arange(0,8), unpack=True, skip_header=skip_header)
    return (l, b, skyt, skyA, H, He, N, Fe)

def get_all_data(ntime, directory = "./", skip_header=3, subdir = "output_000"):
    '''
    get skymap data run from cr-reverb and stored in a directory
    '''    
    fname = "skymap_{:08d}.out".format(ntime)
    data = np.genfromtxt("{}/{}/{}".format(directory, subdir, fname), unpack=True, skip_header=skip_header)
    return (data)

def noisy_map(N=1000000, nside = 64):

    costheta = 2.0 * np.random.random(size=N) - 1.0
    theta = np.arccos(costheta)
    phi = np.random.random(size=N) * 2.0 * np.pi
    hpx_map = hpx_map_from_thetaphi2(theta, phi, nside=nside)
        
    hpx_map[np.isnan(hpx_map)] = 1e-20
    hpx_map[hpx_map <= 0] = 1e-20
    
    return (hpx_map)


def hpx_map_from_thetaphi2(theta, phi, nside=16):
    npix = hp.nside2npix(nside)
    indices = hp.ang2pix(nside, theta, phi)

    idx, counts = np.unique(indices, return_counts=True)

    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts
    
    return (hpx_map)

def hpx_one_map_from_thetaphi(theta, phi, f, nside=16, noise = False, N_noise = 10000000, noise_weight = 0.1):
    npix = hp.nside2npix(nside)
    indices = hp.ang2pix(nside, theta, phi)
    
    hpx_map = np.zeros(npix, dtype=float)
    for i in range(len(f)):   
        hpx_map[indices[i]] += f[i] 
        
    hpx_map[(hpx_map <= 0)] = 1e-20

    return (hpx_map)

def hpx_map_from_thetaphi(theta, phi, f, A, nside=16, noise = False, N_noise = 10000000, noise_weight = 0.1):
    npix = hp.nside2npix(nside)
    indices = hp.ang2pix(nside, theta, phi)
    
    hpx_map = np.zeros(npix, dtype=float)
    hpx_mapA = np.zeros(npix, dtype=float)
    counts = np.zeros_like(f, dtype = int)
    
    # costh = np.cos(0.5 * np.pi - theta)
    
    for i in range(len(f)):   
        hpx_map[indices[i]] += f[i] 
        hpx_mapA[indices[i]] += A[i]
        counts[i] += 1
        
    hpx_map[(hpx_map <= 0)] = 1e-20
    hpx_mapA = hpx_mapA / hpx_map
    hpx_mapA[(hpx_mapA <= 1)] = 1
    
    if noise: 
        hpx_map_noisy = noisy_map(N = N_noise, nside=nside)
        norm = np.sum(hpx_map) / np.sum(hpx_map_noisy)
        print (np.sum(hpx_map), np.sum(hpx_map_noisy), np.sum(hpx_map_noisy * norm * noise_weight))
        hpx_map += (hpx_map_noisy * norm * noise_weight)

#     for i in range(len(hpx_map)):
#         if hpx_map[i] <= 0:
#             hpx_map[i] = -99
#         else:
#             hpx_map[i] = np.log10(hpx_map[i])
            #print (hpx_map[i])
    #idx, counts = np.unique(indices, return_counts=True)
#     for i in range(len(hpx_mapA)):
#         if hpx_mapA[i] >= 1:
#             hpx_mapA[i] = np.log(hpx_mapA[i])
#         else:
#             hpx_mapA[i] = -1

    return (hpx_map, hpx_mapA)


def get_maps_from_data(ntime, nside = 64, directory = ".", smooth = None, N_noise=10000, noise = False, noise_weight = 0.2, skip_header=0):
   
    groups = dict()
    grp_map = dict()
    l, b, skyt, skyA, groups["H"], groups["He"], groups["N"], groups["Fe"] = get_data(ntime, directory=directory, skip_header=skip_header)
    theta = 0.5 * np.pi - np.radians(b.flatten())
    phi = np.radians(l.flatten())
    flux = skyt.flatten()
    mass = skyA.flatten()
    ratio1 = groups["He"]/groups["Fe"]
    ratio2 = groups["He"]/(groups["Fe"]+groups["He"])
    ratio1 = ratio1.flatten()
    ratio2 = ratio2.flatten()
    
    # flatten arrays for groups
    for key in groups.keys():
        groups[key] = groups[key].flatten()
        grp_map[key] = hpx_one_map_from_thetaphi(theta, phi, groups[key], nside=nside)

    ratio1 = grp_map["He"]/grp_map["Fe"]
    ratio2 = grp_map["He"]/(grp_map["Fe"]+grp_map["He"])
    #ratio2 = hpx_one_map_from_thetaphi(theta, phi, ratio2, nside=nside)
    hpx_map, hpx_mapA = hpx_map_from_thetaphi(theta, phi, flux, mass, nside=nside, noise=noise, N_noise=N_noise, noise_weight=noise_weight)
    
    if smooth is not None:
        if smooth > 0:
            hpx_map = hp.smoothing(hpx_map, fwhm=np.radians(smooth))
            hpx_map[hpx_map<=0] = 1e-20
            hpx_mapA = hp.smoothing(hpx_mapA, fwhm=np.radians(smooth))
            
            for key in groups.keys():
                grp_map[key] = hp.smoothing(grp_map[key], fwhm=np.radians(smooth))
                grp_map[key][grp_map[key]<=0] = 1e-20
                
            ratio1 = hp.smoothing(ratio1, fwhm=np.radians(smooth))
            ratio2 = hp.smoothing(ratio2, fwhm=np.radians(smooth))
            #ratio1[ratio1<=0] = ratio1
            #ratio2[ratio2<=0] = ratio2

        
    hpx_mapA[np.isnan(hpx_mapA)] = -1
    hpx_map[np.isnan(hpx_map)] = 1e-20
    hpx_map[hpx_map <= 0] = 1e-20
    
    return (hpx_map, hpx_mapA, grp_map, ratio1, ratio2)


# def bin_groups(lnA_bins = np.arange(0,4.5,0.5)

def get_maps_from_data_MWProp(ntime, nside = 64, directory = ".", smooth = None, subdir = "output_000"):
   
    groups = dict()
    grp_map = dict()
    all_data = get_all_data(ntime, directory=directory, subdir = subdir)
    l = all_data[0]
    b = all_data[1]
    skyt = all_data[2]
    sky_lowlnA = all_data[3] 
    sky_hilnA = all_data[4] 
               
    theta = 0.5 * np.pi - np.radians(b.flatten())
    phi = np.radians(l.flatten())
    flux = skyt.flatten()
    low = sky_lowlnA.flatten()
    high = sky_hilnA.flatten()
    
    # flatten arrays for groups
    hpx_map = hpx_one_map_from_thetaphi(theta, phi, flux, nside=nside)
    high_map = hpx_one_map_from_thetaphi(theta, phi, high, nside=nside)
    low_map = hpx_one_map_from_thetaphi(theta, phi, low, nside=nside)
        
    if smooth is not None:
        if smooth > 0:
            hpx_map = hp.smoothing(hpx_map, fwhm=np.radians(smooth))
            hpx_map[hpx_map<=0] = 1e-20
            
            high_map = hp.smoothing(high_map, fwhm=np.radians(smooth))
            high_map[high_map<=0] = 1e-20
            
            low_map = hp.smoothing(low_map, fwhm=np.radians(smooth))
            low_map[low_map<=0] = 1e-20
            
                    
    hpx_map[np.isnan(hpx_map)] = 1e-20
    low_map[np.isnan(low_map)] = 1e-20
    high_map[np.isnan(high_map)] = 1e-20
   
    return (hpx_map, low_map, high_map)


def get_maps_from_data2(ntime, nside = 64, directory = ".", smooth = None, N_noise=10000, noise = False, noise_weight = 0.2, use_lnA_bins = True, dlnA = 0.5, summed=False):
   
    groups = dict()
    grp_map = dict()
    all_data = get_all_data(ntime, directory=directory)
    l = all_data[0]
    b = all_data[1]
    skyt = all_data[2]
    skyA = all_data[3] 
    if use_lnA_bins:
        lnA_bins = np.arange(0, 4.1 + dlnA, dlnA)
        for i in range(len(lnA_bins)):
            groups[str(i)] = np.zeros_like(skyt)
    else: 
        lnA_bins = None
               
     
    for i in range(4, len(all_data)):
        if use_lnA_bins:
            lnA = np.log(float(i-3))
            ibin = int(lnA / dlnA)
            groups[str(ibin)] += all_data[i,]
        else:
            groups[str(i-3)] = all_data[i,]
    
    
    theta = 0.5 * np.pi - np.radians(b.flatten())
    phi = np.radians(l.flatten())
    flux = skyt.flatten()
    mass = skyA.flatten()
    
    # flatten arrays for groups
    for key in groups.keys():
        groups[key] = groups[key].flatten()
        grp_map[key] = hpx_one_map_from_thetaphi(theta, phi, groups[key], nside=nside)
        
    summed_groups = np.zeros_like(groups["1"])
    for key in groups.keys():
        summed_groups += groups[key]
        
    summed_groups = hpx_one_map_from_thetaphi(theta, phi, summed_groups, nside=nside)
        
    #ratio2 = hpx_one_map_from_thetaphi(theta, phi, ratio2, nside=nside)
    hpx_map, hpx_mapA = hpx_map_from_thetaphi(theta, phi, flux, mass, nside=nside, noise=noise, N_noise=N_noise, noise_weight=noise_weight)
    
    if smooth is not None:
        if smooth > 0:
            hpx_map = hp.smoothing(hpx_map, fwhm=np.radians(smooth))
            hpx_map[hpx_map<=0] = 1e-20
            hpx_mapA = hp.smoothing(hpx_mapA, fwhm=np.radians(smooth))
            
            summed_groups = hp.smoothing(summed_groups, fwhm=np.radians(smooth))
            
            for key in groups.keys():
                grp_map[key] = hp.smoothing(grp_map[key], fwhm=np.radians(smooth))
               
    for key in groups.keys():
        grp_map[key][grp_map[key]<=0] = 1e-20
        grp_map[key][np.isnan(grp_map[key])] = 1e-20
        
    hpx_mapA[np.isnan(hpx_mapA)] = -1
    hpx_map[np.isnan(hpx_map)] = 1e-20
   
    hpx_map[hpx_map <= 0] = 1e-20
    summed_groups[summed_groups <= 0] = 1e-20
        
    if summed:
        return (hpx_map, hpx_mapA, grp_map, lnA_bins, summed_groups)
    else:
        return (hpx_map, hpx_mapA, grp_map, lnA_bins)



def intensity_map(hpx_map, log = False, **kwargs):
    if log:
        projview(np.log10(hpx_map), 
                     unit = "UHECR Flux", 
                     cb_orientation = "vertical", 
                     graticule = True, 
                     projection_type="hammer", 
                     graticule_labels= True,
                     xtick_label_color= "w", 
                     title = title, min=-3, cmap="cividis")
    else:
        projview(hpx_map,
                     unit = "UHECR Flux (Arb.)", 
                     cb_orientation = "vertical", 
                     graticule = True,
                     projection_type="hammer", 
                     graticule_labels= True,
                     xtick_label_color= "w", **kwargs)
        
    return 0

def get_dipoles(roots, ntimes):
    lon = np.zeros((2,len(ntimes)))
    lat = np.zeros((2,len(ntimes)))
    strength = np.zeros((2,len(ntimes)))
    for j, root in enumerate(roots):
        for i, ntime in enumerate(ntimes):
            print (ntime, end=" ")
            hpx_map, hpx_mapA = get_maps_from_data(int(ntime), nside = 64, smooth=smooth, directory=root)
            strength[j,i], vec = hp.pixelfunc.fit_dipole(hpx_map)
            lon[j,i], lat[j,i] = hp.pixelfunc.vec2ang(vec, lonlat = True)
    return (lon, lat, strength)

def make_skymap_plots(root, times, dt_plot, smooths = 20, 
                      maxes = [87.2, 1.01, 0.17, 0.17, 0.17], mins = None,
                      override_plot_properties = {"figure_width": 9, "figure_size_ratio": 5./9.},
                      fontsize_dict = {"xtick_label": 18, "ytick_label": 18, "title": 18, "cbar_label": 18},
                      skip_header=0, savefolder=".", folder=".", renorm = 1.0):
    
    # convert smooth to iterable if needed
    if np.isscalar(smooths):
        smooths = [smooths]
    
    if mins == None:
        mins = np.zeros_like(times)
        
    for smooth in smooths:
        for i,ntime in enumerate(tqdm(times)):
            fig = plt.figure()
            sub = "31{:d}".format(i+1)
            t_plot = ntime * dt_plot

            directory = "{}/{}".format(folder, root)
            
            # get the maps from data 
            hpx_map, hpx_mapA, _, _, _ = get_maps_from_data(ntime, nside = 64, skip_header=skip_header, smooth=smooth, noise=False, noise_weight = 10.0, directory = directory)
            hpx_map[hpx_map < 3e-8] = 1e-10
            if mins[i] == 0.0 and np.max(hpx_map) == 1e-10:
                maxes[i] = 1.0

            time_string = r"$t={:.1f}~{{\rm Myr}}$".format(t_plot)
            intensity_map(hpx_map * renorm, title = time_string, override_plot_properties=override_plot_properties, 
                         fontsize=fontsize_dict, max=maxes[i], min=mins[i])
            plt.xticks(fontsize=16)
            plt.subplots_adjust(top=0.92, bottom=0.08, right=0.96,left=0.04)
            plt.savefig("{}/map_{}_{}_sm{}.png".format(savefolder, root, ntime, smooth), dpi=300, transparent=True)
            
    return 
