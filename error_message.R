############################################################
# function giving the error message asociated to an index
#
# 1XX : errors in process_files.R
# 2XX : errors in write_DM_with_mcornec.R
# 3XX : errors in increment_N_CALIB.R
############################################################

error_message <- function(error_id) {
    switch(toString(error_id), 
           "0" = "Success",
           
           "101" = "101 : Descent profile, use -d to override", 
           "102" = "102 : Bad position (QC = 3 or 4), use -p to force position",
           "103" = "103 : No geolocalisation (lat or lon missing), use -p to force position",
           "104" = "104 : Date is missing",
           "105" = "105 : Bad date (QC = 3 or 4)",
           "106" = "106 : No CHLA in variables",
           "107" = "107 : Stuck pressure",
           "108" = "108 : Problem with RESO, possibly due to a partially stuck pressure",
           "109" = "109 : Profile is on the greylist with QC 4, do process_file(profile_actual, ...) for details",
           "110" = "110 : Profile is on the greylist with an unmanaged QC",
           "111" = "111 : Profile is on the greylist with QC 3, use -q to override",
		   "112" = "112 : The calculated dark offset is too far from factory",
           
           "201" = "201 : BD file given as input",
           "202" = "202 : Multiple profiles of CHLA/BBP700 in the profile",
           "203" = "203 : NPQ correction and dark correction cross pressure domain", 
           "204" = "204 : Dimensions of imported arrays are not aligned as expected",
           "205" = "205 : No CHLA or no BBP700 in variables",
           "206" = "206 : Constructed file names are not unique",
           
           "306" = "306 : Constructed file names are not unique",
           
           paste("No message associated to error_id =", error_id))
}
