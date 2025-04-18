import cogsworth
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import gala.dynamics as gd

import plotting

plt.rc('font', family='serif')
plt.rcParams['text.usetex'] = False
fs = 24

# update various fontsizes to match
params = {'figure.figsize': (12, 8),
          'legend.fontsize': fs,
          'axes.labelsize': fs,
          'xtick.labelsize': 0.9 * fs,
          'ytick.labelsize': 0.9 * fs,
          'axes.linewidth': 1.1,
          'xtick.major.size': 7,
          'xtick.minor.size': 4,
          'ytick.major.size': 7,
          'ytick.minor.size': 4}
plt.rcParams.update(params)


# star_particles = pd.read_hdf("/mnt/home/twagg/supernova-feedback/data/FIRE_star_particles.h5")
# star_particles["ind"] = np.arange(len(star_particles))
# particle_orbits = np.load("/mnt/home/twagg/supernova-feedback/data/particle_orbits.npy", allow_pickle=True)

star_particles = pd.read_hdf("/mnt/home/twagg/supernova-feedback/data/r442_init_recent_stars.h5")
star_particles["ind"] = np.arange(len(star_particles))
particle_orbits = np.load("/mnt/home/twagg/supernova-feedback/data/r442-particle_orbits.npy", allow_pickle=True)

folder = "feedback-variations"
# file_names = ["beta-0.0", "beta-0.5", "beta-1.0", "ce-0.1", "ce-10.0"]
# file_names = ["qcritB-0.0", "qcritB-1000.0", "ccsn-20", "ecsn-265", "no-fallback"]
# file_names = ["imf-1.9", "imf-2.7", "porb-0", "porb-minus1"]
# file_names = ["q-minus1", "q-plus1", "v-disp-0.5", "v-disp-5"]
# file_names = ['singles', 'Z-0.5', 'Z-0.2', 'Z-0.1', 'Z-0.05']
# file_names = ["porb-max-3"]
file_names = ["r442-Zx2"]
    
for file_name in file_names:
    print(f"Doing distances for {file_name}")
    p = cogsworth.pop.load(f"/mnt/home/twagg/ceph/pops/{folder}/{file_name}")
    
    plotting.set_sn_subpop_masks(p)
    
    p._orbits_file = f"/mnt/home/twagg/ceph/pops/{folder}/{file_name}.h5"
    
    bad_binaries = np.isin(p.bpp["bin_num"].values, p.duplicate_sn)
    sn_rows = [p.bpp[(p.bpp["evol_type"] == 15) & ~bad_binaries], p.bpp[(p.bpp["evol_type"] == 16) & ~bad_binaries]]
    kicked_nums = [sn_rows[i]["bin_num"].values for i in range(2)]
    kicked_mask = [np.isin(p.bin_nums, kicked_nums[0]), np.isin(p.bin_nums, kicked_nums[1])]
    
    sn_distances = [np.zeros(len(kicked_nums[0])) * u.kpc, np.zeros(len(kicked_nums[1])) * u.kpc]
    sn_locations = [np.zeros((len(kicked_nums[0]), 3)) * u.kpc, np.zeros((len(kicked_nums[1]), 3)) * u.kpc]
    for i in [0, 1]:
        child_orbits = p.primary_orbits[kicked_mask[i]] if i == 0 else p.secondary_orbits[kicked_mask[i]]
        parent_orbits = particle_orbits[star_particles.loc[p.initC.loc[kicked_nums[i]]["particle_id"].values]["ind"].values]
        for j in range(len(kicked_nums[i])):
            parent_orbit = parent_orbits[j]
            child_orbit = child_orbits[j]
            sn_time = sn_rows[i]["tphys"].iloc[j]

            parent_pos = parent_orbit.pos[(parent_orbit.t - parent_orbit.t[0]).to(u.Myr).value < sn_time][-1]
            child_pos = child_orbit.pos[(child_orbit.t - child_orbit.t[0]).to(u.Myr).value < sn_time][-1]

            sn_distances[i][j] = sum((parent_pos - child_pos).xyz**2)**(0.5)
            sn_locations[i][j] = child_pos.xyz
            
    np.savez(f"/mnt/home/twagg/ceph/pops/{folder}/sn_positions-{file_name}",
             sn_distances[0].to(u.kpc).value, sn_distances[1].to(u.kpc).value,
             sn_locations[0].to(u.kpc).value, sn_locations[1].to(u.kpc).value)

    p.primary_sn_distances, p.secondary_sn_distances = sn_distances[0], sn_distances[1]
    
    fig, axes = plotting.sandpile(p, show=False)
    plt.savefig(f"/mnt/home/twagg/supernova-feedback/plots/sn_distance_hists/{file_name}.pdf", format="pdf", dpi=300, bbox_inches="tight")
