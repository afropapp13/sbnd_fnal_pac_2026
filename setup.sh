#!/bin/bash

export TERM=screen

#source /grid/fermiapp/products/uboone/setup_uboone_mcc9.sh
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup uboonecode v08_00_00_52 -q e17:prof

htgettoken -a htvaultprod.fnal.gov -i uboone
