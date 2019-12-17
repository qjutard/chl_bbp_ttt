#!/bin/bash

usage() { 
	echo "Usage: $0 -W <WMO_number> | -L <profile_list> | -P <profile_name> [-D <DEEP_EST>] [-p <position>] [-o <offset>] [-cdfhq]
Do '$0 -h' for help" 1>&2
	exit 1 
}
helprint() {
	echo "
#########################################################################################

DMMC does Delayed mode computing and writing following the work done by M. Cornec 
in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)

Usage: $0 -W <WMO_number> | -L <profile_list> | -P <profile_name> [-D <DEEP_EST>] [-p <position>] [-o <offset>] [-cdfhq]

### Options

-W <WMO_number> : Do the delayed mode on all profiles of a float identified with its
                  7 digits WMO number.
-L <profile_list> : Do the delayed mode on the profiles identified in a list (text file)
                    with the format 'XXXXXXX_YYYZ' where XXXXXXX is the WMO number,
                    YYY is the profile number, ans Z is one of '.' or 'D' depending on
                    whether the profile is ascending or descending.
-P <profile_name> : Do the delayed mode on the profile <profile_name> with the same
                    format as described in -L.
[-D <DEEP_EST>] : Use an already existing DEEP_EST table. This table can take some time
                  to be computed so if DMMC has already been used on this float it is
                  best practice to reuse the DEEP_EST table that has been created.
[-c] : Just copy the profiles from the input directory to the output directory.
[-d] : Accept descent profile.
[-f] : Fill the delayed mode profiles with fill values and bad QC.
[-h] : help
[-o <offset>] : Override the computed dark offset and pressure minimum for dark. <offset>
                should be formatted as 'OFF.off;MIN.min' with the single brackets, where
                where OFF.off is the desired offset (chl_dark=chl-offset) and MIN.min is
                a pressure from which all chl_dark values are set to 0. Use 'NA' in any
                or all of these parameters to ignore them (i.e. OFF=0 and MIN=+inf).
[-p <position>] : Override the profile position in the case of a bad QC flag ('3' or
                  '4'). <position> should be formatted as 'LAT.lat;LON.lon' with the
                  single brackets.
[-q] : Accept profiles on the greylist with QC='3'.

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
qc3=FALSE
Profile=NA
position=NA
offset=NA

while getopts W:L:D:cfdqP:p:o:h option
do
case "${option}"
in
W) WMO=${OPTARG};;
L) List=${OPTARG};;
D) DEEP=${OPTARG};;
c) copy=TRUE;;
f) fill=TRUE;;
d) descent=TRUE;;
q) qc3=TRUE;;
P) Profile=${OPTARG};;
p) position=${OPTARG};;
o) offset=${OPTARG};;
h) helprint;;
*) usage;;
esac
done

Rscript ~/Documents/cornec_chla_qc/chl_bbp_ttt/start_DMMC.R $WMO $List $DEEP $copy $fill $descent $qc3 $Profile $position $offset
