#!/bin/bash

# don't source anything up
# the bash script does that

export tune="AR23_20i_00_000"

export target="1000180400"

export ghep_dir="/pnfs/sbnd/persistent/users/kplows/GENIE-GEN/for_PAC/v3_06_02_sbn2/from_histogram/${tune}/"
export outdir="/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/sbnd_fnal_pac_2026/${tune}/"

# Loop over all ghep.root files
#for infile in ${ghep_dir}/*.ghep.root; do

#    # Get filename without path
#    basefile=$(basename ${infile} .ghep.root)

#    echo "Processing ${basefile}..."

#    /exp/sbnd/app/users/gputnam/Ar23-knobs/update_ghep_reweight_anywhere.sh \
#        ${infile} \
#        ${outdir}/syst_${basefile}.root


#done

hadd -f ${outdir}/syst_${tune}.root ${outdir}/syst_*.root
