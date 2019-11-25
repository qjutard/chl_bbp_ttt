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
    
    files<-as.character(index_ifremer[,1]) #retrieve the path of each netcfd file
    ident<-strsplit(files,"/") #separate the different roots of the files paths
    ident<-matrix(unlist(ident), ncol=4, byrow=TRUE)
    prof_id<-ident[,4] #retrieve all profiles  name as a vector
    
    i <-which(substr(prof_id,3,14)==profile_actual) #identify profile position in the index
    
    ############################
    ### Get the chla and bbp corrections from the method
    ############################
    L = process_file(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=DEEP_EST)
    if (is.null(L)){
        print("process_file(...) did not end properly")
        return(NULL)
    }
    CHLA_ADJUSTED = L$CHLA_ADJUSTED
    BBP700_ADJUSTED = L$BBP700_ADJUSTED
    CHLA_ADJUSTED_QC = L$CHLA_ADJUSTED_QC
    BBP700_ADJUSTED_QC = L$BBP700_ADJUSTED_QC
    CHLA_ADJUSTED_ERROR = L$CHLA_ADJUSTED_ERROR
    BBP700_ADJUSTED_ERROR = L$BBP700_ADJUSTED_ERROR

    ############################
    ### Open input and output files
    ############################
    path_split = unlist( strsplit(files[i],"/") )
    path_to_profile = paste(path_split[1], path_split[2], path_split[3], sep="/")
    
    filenc_name_M = path_split[4]
    filenc_name_B = paste("B?",substring(filenc_name_M, 3),sep="")
    filenc_name_out = paste("BD",substring(filenc_name_M, 3,nchar(filenc_name_M)-3),"_test.nc",sep="")
    
    file_B = paste(path_to_netcdf, path_to_profile,"/", filenc_name_B, sep="") 
    file_B = system2("ls",file_B,stdout=TRUE) # identify R or D file 
    file_out = paste(path_to_netcdf, path_to_profile,"/", filenc_name_out, sep="") 
    
    ### check whether the existing file is a D file
    path_split_B = unlist( strsplit(file_B, "/") )
    B_type = substr(path_split_B[length(path_split_B)],1,2)
    if (B_type=="BD"){
        print("a BD-file has been given as an input, these are currently not accepted")
        return(NULL)
    }
    
    ### create the output file as a copy of the input
    system2("cp", paste(file_B, file_out))
    
    #filenc_in <- nc_open(file_B, readunlim=FALSE, write=FALSE)
    filenc_out <- nc_open(file_out, readunlim=FALSE, write=TRUE)
    
    ############################
    ### Write correction from method to BD-file
    ############################
    ncvar_put(filenc_out,"CHLA_ADJUSTED",CHLA_ADJUSTED)
    ncvar_put(filenc_out,"BBP700_ADJUSTED",BBP700_ADJUSTED)
    ncvar_put(filenc_out,"CHLA_ADJUSTED_QC",CHLA_ADJUSTED_QC)
    ncvar_put(filenc_out,"BBP700_ADJUSTED_QC",BBP700_ADJUSTED_QC)
    ncvar_put(filenc_out,"CHLA_ADJUSTED_ERROR",CHLA_ADJUSTED_ERROR)
    ncvar_put(filenc_out,"BBP700_ADJUSTED_ERROR",BBP700_ADJUSTED_ERROR)
    
    #nc_close(filenc_in)
    nc_close(filenc_out)
    
    return(0)
    
}

M = write_DM_MC(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST = DEEP_EST)
