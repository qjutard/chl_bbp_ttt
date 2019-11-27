############################################################
# function giving the error message asociated to an index
#
# 1XX : errors in process_files.R
# 2XX : errors in write_DM_with_mcornec.R
############################################################

error_message <- function(error_id) {
    switch(toString(error_id), 
           "0" = "Success",
           
           "101" = "Descent Profile", 
           "102" = "Bad position (QC = 3 or 4)",
           "103" = "No geolocalisation (lat or lon missing)",
           "104" = "Date is missing",
           "105" = "Bad date (QC = 3 or 4)",
           "106" = "No CHLA in variables",
           "107" = "Stuck pressure",
           
           "201" = "BD file given as input",
           "202" = "Multiple profiles of CHLA/BBP700 in the profile",
           
           paste("No message associated to error_id =", error_id))
}
