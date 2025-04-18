import argparse
import astropy.units as u
import cogsworth
import pandas as pd


def run_variation(file_name=None, alphavir=None, radius=None):
    print("Loading in the template")
    p = cogsworth.pop.load("/mnt/home/twagg/ceph/pops/feedback-variations/fiducial.h5",
                           parts=["initial_binaries", "initial_galaxy", "stellar_evolution"])
    
    print("Adjusting settings")
    
    particle_ids = p.initC["particle_id"]

    # excuse some minor hacking around the function nothing to see here
    setattr(p, "sample_initial_galaxy", cogsworth.hydro.pop.HydroPopulation.sample_initial_galaxy)
    p._initial_binaries = p.initC
    p.star_particles = pd.read_hdf("/mnt/home/twagg/supernova-feedback/data/FIRE_star_particles.h5")
    p._subset_inds = p.star_particles.index.values

    p._file = None

    if alphavir is not None:
        p.virial_parameter = alphavir

    if radius is not None:
        p.cluster_radius = radius * u.pc
    else:
        p.cluster_radius = 3 * u.pc
    p.cluster_mass = 1e4 * u.Msun

    print("Resampling initial galaxy")
    p.sample_initial_galaxy(p)
    
    print("Starting galactic evolution")
    
    p.perform_galactic_evolution()
    
    p.initC["particle_id"] = particle_ids.loc[p.bin_nums]

    if file_name is None or file_name == "":
        file_name = f"alphavir-{alphavir}-r-{radius}.h5"

    p.save(f"/mnt/home/twagg/ceph/pops/feedback-variations/{file_name}", overwrite=True)


def main():
    parser = argparse.ArgumentParser(description='Feedback galaxy variation runner')
    parser.add_argument('-f', '--file', default="", type=str,
                        help='File name to use')
    parser.add_argument('-a', '--alphavir', default=None, type=float)
    parser.add_argument('-r', '--radius', default=None, type=float)
    args = parser.parse_args()

    run_variation(file_name=args.file, alphavir=args.alphavir, radius=args.radius)

if __name__ == "__main__":
    main()
