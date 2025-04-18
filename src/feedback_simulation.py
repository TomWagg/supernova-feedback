import argparse
import gala.potential as gp
import astropy.units as u
import pandas as pd
import numpy as np
import cogsworth

def run_sim(high_mass_slope=None,
            porb_model="sana12",
            q_power_law=0,
            binfrac=1.0,
            processes=128,
            extra_time=200 * u.Myr,
            file=None):
    pot = gp.load("/mnt/home/twagg/supernova-feedback/data/m11h_new_potential.yml")
    star_particles = pd.read_hdf("/mnt/home/twagg/supernova-feedback/data/FIRE_star_particles.h5", key="df")
    # subset = np.load("/mnt/home/twagg/supernova-feedback/data/subset.npy")

    sampling_params = {"porb_model": porb_model, "q_power_law": q_power_law}

    if high_mass_slope is not None:
        sampling_params["primary_model"] = "custom"
        sampling_params["alphas"] = [-1.3, -2.3, high_mass_slope]
        sampling_params["mcuts"] = [0.08, 0.5, 1.0, 150.]
    
    BSE_settings = {}
    if binfrac is not None:
        sampling_params["binfrac"] = binfrac
        sampling_params["keep_singles"] = binfrac == 0.0
        BSE_settings["binfrac"] = binfrac

    p = cogsworth.hydro.pop.HydroPopulation(star_particles=star_particles,
                                            max_ev_time=13736.52127883025 * u.Myr + extra_time,
                                            galactic_potential=pot,
                                            m1_cutoff=4,
                                            virial_parameter=1.0,
                                            subset=None,
                                            cluster_radius=3 * u.pc,
                                            cluster_mass=1e4 * u.Msun,
                                            processes=processes,
                                            sampling_params=sampling_params,
                                            BSE_settings=BSE_settings)
    p.sample_initial_binaries()
    p.sample_initial_galaxy()

    # if this is the fiducial model, save a template
    if high_mass_slope is None and porb_model == "sana12" and q_power_law == 0 and binfrac == 1.0 and extra_time==200 * u.Myr:
        p.save(f"/mnt/home/twagg/ceph/pops/feedback-variations/variation-template", overwrite=True)

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
    if file is None:
        file = f"feedback-sim-porb-{porb_model}-q-{q_power_law}"

    p.save(f"/mnt/home/twagg/ceph/pops/feedback-variations/{file}", overwrite=True)


def main():
    parser = argparse.ArgumentParser(description='Supernova feedback simulation runner')
    parser.add_argument('-p', '--processes', default=128, type=int,
                        help='Number of processes to use')
    parser.add_argument('-t', '--extra_time', default=200, type=int,
                        help='Extra time to evolve for (in Myr)')
    parser.add_argument('-f', '--file', default=None, type=str,
                        help='Output file name')
    parser.add_argument('-P', '--porb_model', default="sana12", type=str,
                        help='Binary orbital period model')
    parser.add_argument('-q', '--q_power_law', default=0, type=float,
                        help='Binary mass ratio power law')
    parser.add_argument('-b', '--binfrac', default=1.0, type=float,
                        help='Binary fraction')
    parser.add_argument('-m', '--high-mass-slope', default=None, type=float,
                        help='High mass slope')
    parser.add_argument('-M', '--porb-max', default=None, type=float,
                        help='Maximum log10 porb')
    args = parser.parse_args()

    # check if args.porb_model is a number and convert to dict if so
    try:
        args.porb_model = {
            "min": 0.15,
            "max": 5.5 if args.porb_max is None else args.porb_max,
            "slope": float(args.porb_model)
        }
    except ValueError:
        pass

    run_sim(high_mass_slope=args.high_mass_slope,
            porb_model=args.porb_model,
            q_power_law=args.q_power_law,
            binfrac=args.binfrac,
            file=args.file,
            processes=args.processes,
            extra_time=args.extra_time * u.Myr)

if __name__ == "__main__":
    main()
