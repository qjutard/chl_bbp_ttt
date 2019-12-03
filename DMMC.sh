#!/bin/bash

usage() { 
	echo "Usage: $0 -W <WMO_number> | -L <profile_list> [-D <DEEP_EST>] [-cdfh]
Do '$0 -h' for help" 1>&2
	exit 1 
}
helprint() {
	echo "
#########################################################################################

DMMC does Delayed mode computing and writing following the work done by M. Cornec 
in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)

Usage: $0 -W <WMO_number> | -L <profile_list> [-D <DEEP_EST>] [-cdfh]

### Options

-W <WMO_number> : Do the delayed mode on all profiles of a float identified with its
                  7 digits WMO number
-L <profile_list> : Do the delayed mode on the profiles identified in a list (text file)
                    with the format 'XXXXXXX_YYYZ' where XXXXXXX is the WMO number,
                    YYY is the profile number, ans Z is one of '.' or 'D' depending on
                    whether the profile is ascending or descending
[-D <DEEP_EST>] : Use an already existing DEEP_EST table. This table can take some time
                  to be computed so if DMMC has already been used on this float it is
                  best practice to reuse the DEEP_EST table that has been created
[-c] : Just copy the profiles from the input directory to the output directory
[-d] : Accept descent profile
[-f] : Fill the delayed mode profiles with fill values and bad QC
[-h] : help

#########################################################################################
" 1>&2
	exit 0
}

WMO=NA
List=NA
DEEP=NA
copy=FALSE
fill=FALSE
descent=FALSE

while getopts W:L:D:cfdh option
do
case "${option}"
in
W) WMO=${OPTARG};;
L) List=${OPTARG};;
D) DEEP=${OPTARG};;
c) copy=TRUE;;
f) fill=TRUE;;
d) descent=TRUE;;
h) helprint;;
*) usage;;
esac
done


#echo ${List}
#echo ${copy}
#echo ${WMO}

Rscript start_DMMC.R $WMO $List $DEEP $copy $fill $descent
