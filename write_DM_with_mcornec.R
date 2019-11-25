#########################
# This script allows to read BR files (later BD files) and to write a BD file with a delayed mode (DM)
# for CHLA and BBP700 according to the work done by M. Cornec
#
# The script is written as a function to ideally loop over the profiles of a float or parallelize eg with mcmapply
#########################

rm(list = ls())

library(ncdf4) #deal with netcdf format files
library(oce) #calculate density sigma
library(MASS)

source("process_files.R")

index_ifremer<-read.table("~/Documents/data/argo_merge-profile_index.txt", skip=9, sep = ",")
path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"

profile_list<-c("6901524_150.")
profile_actual = profile_list[1]

### DEEP_EST should be computed once per FLOAT 
DEEP_EST = Dark_MLD_table_coriolis(substr(profile_actual,1,7), path_to_netcdf, index_ifremer) 


write_DM_MC <- function(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=NULL){
    
    L = process_file(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=DEEP_EST)
    
    return(L)
    
}

M = write_DM_MC(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST = DEEP_EST)
