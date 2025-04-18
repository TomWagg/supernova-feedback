import argparse
import cogsworth

def reset_sampled_kicks(p):
    cols = ["natal_kick_1", "phi_1", "theta_1", "natal_kick_2", "phi_2", "theta_2"]
    for col in cols:
        p._initC[col] = -100.0

def run_variation(file_name=None, alpha_ce=None, mt_eff=None, ecsn_kick=None, bhflag=None, qcritB=None,
                  Z_factor=None, gamma=None, ccsn_kick=None):
    print("Loading in the template")
    p = cogsworth.pop.load("/mnt/home/twagg/ceph/pops/feedback-variations/variation-template.h5",
                           parts=["initial_binaries", "initial_galaxy"])
    
    print("Adjusting settings")
    
    particle_ids = p.initial_binaries["particle_id"]
    
    if alpha_ce is not None:
        p.BSE_settings["alpha1"] = alpha_ce
        # p.initC["alpha1"] = alpha_ce

    if mt_eff is not None:
        p.BSE_settings["acc_lim"] = mt_eff
        # p.initC["acc_lim"] = mt_eff
        
    if ecsn_kick is not None:
        # reset_sampled_kicks(p)
        p.BSE_settings["sigmadiv"] = ecsn_kick
        # p.initC["sigmadiv"] = ecsn_kick
        
    if ccsn_kick is not None:
        # reset_sampled_kicks(p)
        p.BSE_settings["sigma"] = ccsn_kick
        
    if bhflag is not None:
        # reset_sampled_kicks(p)
        p.BSE_settings["bhflag"] = bhflag
        # p.initC["bhflag"] = bhflag

    if qcritB is not None:
        # set qcrit for kstar = 2,3,4
        qcrit_array = [0.0, 0.0, qcritB, qcritB, qcritB,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        p.BSE_settings["qcrit_array"] = qcrit_array
        # p.initC["qcrit_2"] = qcritB
        # p.initC["qcrit_3"] = qcritB
        # p.initC["qcrit_4"] = qcritB

    if Z_factor is not None:
        if p._initial_binaries is not None:
            p._initial_binaries["metallicity"] *= Z_factor
        # p.initC["metallicity"] *= Z_factor
        p._initial_galaxy._Z *= Z_factor

    if gamma is not None:
        p.BSE_settings["gamma"] = gamma
        # p.initC["gamma"] = gamma

    p._file = None

    print("Starting stellar evolution")
        
    p.perform_stellar_evolution()

    exploding_nums = p.bpp[p.bpp["evol_type"].isin([15, 16])]["bin_num"].unique()
    p = p[exploding_nums]
    
    print("Starting galactic evolution")
    
    p.perform_galactic_evolution()
    
    p.initC["particle_id"] = particle_ids.loc[p.bin_nums]

    if file_name is None or file_name == "":
        file_name = f"ce-{alpha_ce}-beta-{mt_eff}-ecsn-{ecsn_kick}-bhflag-{bhflag}.h5"

    print("Saving file now")

    p.save(f"/mnt/home/twagg/ceph/pops/feedback-variations/{file_name}", overwrite=True)


def main():
    parser = argparse.ArgumentParser(description='Feedback variation runner')
    parser.add_argument('-f', '--file', default="", type=str,
                        help='File name to use')
    parser.add_argument('-c', '--alpha_ce', default=None, type=float,
                        help='Common-envelope efficiency')
    parser.add_argument('-b', '--beta', default=None, type=float,
                        help='Mass transfer efficiency')
    parser.add_argument('-k', '--bhflag', default=None, type=int,
                        help='BH kick flag')
    parser.add_argument('-e', '--ecsn-kick', default=None, type=int,
                        help='ECSN kick strength')
    parser.add_argument('-C', '--ccsn-kick', default=None, type=int,
                        help='CCSN kick strength')
    parser.add_argument('-q', '--qcritB', default=None, type=float,
                        help='Critical mass ratio for kstar = 2,3,4')
    parser.add_argument('-Z', '--Z_factor', default=None, type=float,
                        help='Metallicity factor')
    parser.add_argument('-g', '--gamma', default=None, type=float,
                        help='Angular momentum loss for mass transfer')
    args = parser.parse_args()

    run_variation(file_name=args.file,
                  alpha_ce=args.alpha_ce, mt_eff=args.beta, ecsn_kick=args.ecsn_kick,
                  ccsn_kick=args.ccsn_kick, bhflag=args.bhflag, qcritB=args.qcritB,
                  Z_factor=args.Z_factor, gamma=args.gamma)

if __name__ == "__main__":
    main()
