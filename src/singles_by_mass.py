import cogsworth

print("Loading in the template")
p = cogsworth.pop.load("/mnt/home/twagg/ceph/pops/feedback-variations/variation-template.h5",
                       parts=["initial_binaries", "initial_galaxy"])

print("Adjusting settings")

particle_ids = p.initial_binaries["particle_id"]

p.initial_binaries["porb"] = 1e20
p.initial_binaries["ecc"] = 0.0
print("Starting stellar evolution")
    
p.perform_stellar_evolution()

exploding_nums = p.bpp[p.bpp["evol_type"].isin([15, 16])]["bin_num"].unique()
p = p[exploding_nums]

print("Starting galactic evolution")

p.perform_galactic_evolution()

p.initC["particle_id"] = particle_ids.loc[p.bin_nums]

print("Saving file now")

p.save(f"/mnt/home/twagg/ceph/pops/feedback-variations/singles", overwrite=True)
