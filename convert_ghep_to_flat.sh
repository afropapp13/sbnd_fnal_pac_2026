#!/bin/bash

# don't forget to first source the setup script that includes nuisance

#export tune="AR23_20i_00_000"
#export tune="AR25_20i_00_000"
#export tune="AR25_20j_00_000"

export tune="AR25_20l_00_000"

export target="1000180400"

export fluxfile="/pnfs/sbnd/persistent/users/apapadop/Fluxes/sbnd_flux.root"
export fluxhisto="flux_sbnd_numu"

export ghep_dir="/pnfs/sbnd/persistent/users/kplows/GENIE-GEN/for_PAC/v3_06_02_sbn2/from_histogram/${tune}/"
export outdir="/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/sbnd_fnal_pac_2026/${tune}/"

# Loop over all ghep.root files
for infile in ${ghep_dir}/*.ghep.root; do

    # Get filename without path
    basefile=$(basename ${infile} .ghep.root)

    echo "Processing ${basefile}..."

    # Convert file from ghep to nuisance format
    PrepareGENIE \
        -i ${infile} \
        -t ${target}[1] \
        -o ${outdir}/${basefile}.gprep.root \
        -f ${fluxfile},${fluxhisto}

    # Convert to nuisance flat tree format
    nuisflat \
        -i GENIE:${outdir}/${basefile}.gprep.root \
        -o ${outdir}/${basefile}.flat.root

done

# Remove unnecessary files
rm -f input-flux.root

hadd -f ${outdir}/all_${tune}.flat.root ${outdir}/*.flat.root