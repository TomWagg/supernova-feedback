import argparse
import numpy as np
import gala.potential as gp
import astropy.units as u
import pandas as pd

import cogsworth

def run_sim(alpha_vir, alpha_ce, mt_eff=-1, ecsn_kick=-20, bhflag=1, subset=None, processes=128, extra_time=200 * u.Myr):
    pot = gp.load("/mnt/home/twagg/supernova-feedback/data/r442_new_potential.yml")
    star_particles = pd.read_hdf("/mnt/home/twagg/supernova-feedback/data/r442_init_recent_stars.h5", key="df")
    
    # increase r442 metallicity to make FIRE m11h
    star_particles["Z"] = star_particles["Z"] * 2

    p = cogsworth.hydro.pop.HydroPopulation(star_particles=star_particles,
                                            max_ev_time=13800.820457297743 * u.Myr + extra_time,
                                            galactic_potential=pot,
                                            m1_cutoff=4,
                                            virial_parameter=alpha_vir,
                                            subset=subset,
                                            cluster_radius=3 * u.pc,
                                            cluster_mass=1e4 * u.Msun,
                                            processes=processes,
                                            BSE_settings={"binfrac": 1.0, "alpha1": alpha_ce, "acc_lim": mt_eff, 'bhflag': bhflag, 'sigmadiv': ecsn_kick})
    p.sample_initial_binaries()
    p.sample_initial_galaxy()

    print("Initial binaries and galaxy sampled, performing stellar evolution")

    p.perform_stellar_evolution()
    
    print("Stellar evolution complete, performing galactic evolution")
    exploding_nums = p.bpp[p.bpp["evol_type"].isin([15, 16])]["bin_num"].unique()
    p = p[exploding_nums]
    
    print("Masked out systems that don't reach core-collapse")

    p.perform_galactic_evolution()
    
    if p._initC is not None and "particle_id" not in p._initC.columns:
        p._initC["particle_id"] = p._initial_binaries["particle_id"]

    # save the results
    p.save(f"/mnt/home/twagg/ceph/pops/feedback-variations/r442-Zx2", overwrite=True)


def main():
    parser = argparse.ArgumentParser(description='Simulation runner')
    parser.add_argument('-a', '--alpha_vir', default=1.0, type=float,
                        help='Star particle virial parameter')
    parser.add_argument('-c', '--alpha_ce', default=1.0, type=float,
                        help='Common-envelope efficiency')
    parser.add_argument('-b', '--beta', default=-1, type=float,
                        help='Mass transfer efficiency')
    parser.add_argument('-k', '--bhflag', default=1, type=int,
                        help='BH kick flag')
    parser.add_argument('-e', '--ecsn-kick', default=-20, type=int,
                        help='ECSN kick strength')
    parser.add_argument('-s', '--subset', default=None, type=int,
                        help='Size of subset of star particles to use')
    parser.add_argument('-p', '--processes', default=128, type=int,
                        help='Number of processes to use')
    parser.add_argument('-t', '--extra_time', default=200, type=int,
                        help='Extra time to evolve for (in Myr)')
    args = parser.parse_args()

    run_sim(alpha_vir=args.alpha_vir, alpha_ce=args.alpha_ce, mt_eff=args.beta, bhflag=args.bhflag, ecsn_kick=args.ecsn_kick, subset=args.subset, processes=args.processes, extra_time=args.extra_time * u.Myr)

if __name__ == "__main__":
    main()
