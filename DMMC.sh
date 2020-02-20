#!/bin/bash

usage() { 
	echo "Usage: $0 -W <WMO_number> | -L <profile_list> | -P <profile_name> [-D <DEEP_EST>] [-o <offset> | -O <offset_file>] [-p <position>] [-t <date>] [-B|-C] [-cdfhq]
Do '$0 -h' for help" 1>&2
	exit 1 
}
helprint() {
	echo "
#########################################################################################

DMMC does Delayed mode computing and writing following the work done by M. Cornec 
in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)

Usage: $0 -W <WMO_number> | -L <profile_list> | -P <profile_name> [-D <DEEP_EST>] [-o <offset> | -O <offset_file>] [-p <position>] [-t <date>] [-B|-C] [-cdfhq]

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
[-o <offset>] : Override the computed dark offset and pressure minimum for dark. <offset>
                should be formatted as 'OFF.off;MIN.min' with the single brackets, where
                OFF.off is the desired offset (chl_dark=chl-offset) and MIN.min is a 
                pressure from which all chl_dark values are set to 0. Use 'NA' in any or
                all of these parameters to ignore them (i.e. OFF=0 and MIN=+inf). The 
                option can also accept 'dmmc' as an argument to just force the use of the
                offset and min computed by DMMC.
[-O <offset_file>] : Override the computed dark offset with a list of offsets given in
                     <offset_file>. The file should contain in its first column profile 
                     names formatted as described in -L, and in the second column it
                     should contain corresponding offset values. Such files can be
                     obtained with the DARK software. The MIN.min argument described in
                     the -o option is not supported here.
[-p <position>] : Override the profile position or the case of a bad QC flag ('3' or
                  '4'). <position> should be formatted as 'LAT.lat;LON.lon' with the
                  single brackets. This does not change the position in the output file
[-t <date>] : Override the date or date QC. <date> should be formatted as 
              'yyyy-mm-dd;hh:mm:ss' in UTC. This does not change the date in the output
              file.
[-B] : Only do the delayed mode for BBP700, this only stops the writing of CHLA delayed
       mode, errors in the computation of CHLA delayed mode are still raised.
[-c] : Just copy the profiles from the input directory to the output directory.
[-C] : Only do the delayed mode for CHLA, this only stops the writing of BBP700 delayed
       mode, errors in the computation of BBP700 delayed mode are still raised.
[-d] : Accept descent profile.
[-f] : Fill the delayed mode profiles with fill values and bad QC, without running DMMC
       computations.
[-h] : help
[-q] : Accept profiles on the greylist with QC='3', this also limits the QC of adjusted
       parameters to '3' at best.

#########################################################################################
" 1>&2
	exit 0
}

WMO=NA
List=NA
Profile=NA
DEEP=NA
offset=NA
Offset_file=NA
position=NA
time_date=NA
BBP_only=FALSE
CHL_only=FALSE
copy=FALSE
descent=FALSE
fill=FALSE
qc3=FALSE

while getopts W:L:P:D:o:O:p:t:BcCdfqh option
do
case "${option}"
in
W) WMO=${OPTARG};;
L) List=${OPTARG};;
P) Profile=${OPTARG};;
D) DEEP=${OPTARG};;
o) offset=${OPTARG};;
O) Offset_file=${OPTARG};;
p) position=${OPTARG};;
t) time_date=${OPTARG};;
B) BBP_only=TRUE;;
c) copy=TRUE;;
C) CHL_only=TRUE;;
d) descent=TRUE;;
f) fill=TRUE;;
q) qc3=TRUE;;
h) helprint;;
*) usage;;
esac
done

Rscript ~/Documents/cornec_chla_qc/chl_bbp_ttt/start_DMMC.R $WMO $List $DEEP $copy $fill $descent $qc3 $Profile $position $offset $BBP_only $time_date $Offset_file $CHL_only
