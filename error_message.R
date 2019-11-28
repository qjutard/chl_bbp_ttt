############################################################
# function giving the error message asociated to an index
#
# 1XX : errors in process_files.R
# 2XX : errors in write_DM_with_mcornec.R
############################################################

error_message <- function(error_id) {
    switch(toString(error_id), 
           "0" = "Success",
           
           "101" = "101 : Descent profile", 
           "102" = "102 : Bad position (QC = 3 or 4)",
           "103" = "103 : No geolocalisation (lat or lon missing)",
           "104" = "104 : Date is missing",
           "105" = "105 : Bad date (QC = 3 or 4)",
           "106" = "106 : No CHLA in variables",
           "107" = "107 : Stuck pressure",
           "108" = "108 : Problem with RESO, possibly due to a partially stuck pressure",
           
           "201" = "201 : BD file given as input",
           "202" = "202 : Multiple profiles of CHLA/BBP700 in the profile",
           "203" = "203 : NPQ correction and dark correction cross pressure domain", 
           "204" = "204 : The CHLA in the computed NPQ zone is not unique", ###deprecated, to replace
           "205" = "205 : Unexpected formatting of SCIENTIFIC_CALIB_COEFFICIENT", ###deprecated, to replace
           "206" = "206 : Constructed file names are not unique",
           
           paste("No message associated to error_id =", error_id))
}
