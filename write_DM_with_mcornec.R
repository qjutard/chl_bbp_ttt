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
library(stringr)

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
    ### Find some parameter indices and define date
    ############################
    
    parameters=ncvar_get(filenc_out,"STATION_PARAMETERS")
    
    id_param_chl = idc=grep("CHLA                                                            ",parameters)
    id_param_bbp = idc=grep("BBP700                                                          ",parameters)
    id_prof = which(parameters=="CHLA                                                            ", arr.ind=TRUE)[2] # The CHL and BBP share the same profile
    n_prof = dim(parameters)[2]
    
    N_HISTORY=filenc_out$dim[['N_HISTORY']]$len
    i_history=N_HISTORY+1
    
    date = Sys.time()
    date = str_sub(date,1,19)
    date = unlist(strsplit(date, "-"))
    date = unlist(strsplit(date, " "))
    date = unlist(strsplit(date, ":"))
    DATE = paste(date, collapse="")
    
    ############################
    ### Write correction from method to BD-file
    ############################
    
    ncvar_put(filenc_out,"CHLA_ADJUSTED",CHLA_ADJUSTED)
    ncvar_put(filenc_out,"BBP700_ADJUSTED",BBP700_ADJUSTED)
    ncvar_put(filenc_out,"CHLA_ADJUSTED_QC",CHLA_ADJUSTED_QC)
    ncvar_put(filenc_out,"BBP700_ADJUSTED_QC",BBP700_ADJUSTED_QC)
    ncvar_put(filenc_out,"CHLA_ADJUSTED_ERROR",CHLA_ADJUSTED_ERROR)
    ncvar_put(filenc_out,"BBP700_ADJUSTED_ERROR",BBP700_ADJUSTED_ERROR)
    
    ############################
    ### Write scientific_calib
    ############################
    
    ### scientific calib information
    scientific_comment_chl = "sample scientific comment chla" # TODO fill comment
    scientific_comment_bbp = "sample scientific comment bbp700" # TODO fill comment
    scientific_coefficient_chl = "sample scientific coefficient chla" # TODO fill coefficient
    scientific_coefficient_bbp = "sample scientific coefficient bbp700" # TODO fill coefficient
    scientific_equation_chl = "sample scientific equation chla" # TODO fill equation
    scientific_equation_bbp = "sample scientific equation bbp700" # TODO fill equation
    scientific_date_chl = DATE
    scientific_date_bbp = DATE
    
    scientific_calib_chl = c(scientific_comment_chl, scientific_coefficient_chl, scientific_equation_chl, scientific_date_chl)
    scientific_calib_bbp = c(scientific_comment_bbp, scientific_coefficient_bbp, scientific_equation_bbp, scientific_date_bbp)
    SCIENTIFIC_CALIB_VARIABLE = paste("SCIENTIFIC_CALIB_",c("COMMENT", "COEFFICIENT", "EQUATION", "DATE"), sep="")
    
    for (i in seq(1,4)){
        SCIENTIFIC_CALIB_INFO_CHL = str_pad(scientific_calib_chl[i], 256, "right")
        SCIENTIFIC_CALIB_INFO_BBP = str_pad(scientific_calib_bbp[i], 256, "right")
        
        scientific_calib_info = ncvar_get(filenc_out, SCIENTIFIC_CALIB_VARIABLE[i])
        
        scientific_calib_info[id_param_chl] = SCIENTIFIC_CALIB_INFO_CHL
        scientific_calib_info[id_param_bbp] = SCIENTIFIC_CALIB_INFO_BBP
        
        ncvar_put(filenc_out , SCIENTIFIC_CALIB_VARIABLE[i], scientific_calib_info)
    
    }
    
    ############################
    ### Write history
    ############################
    
    ### Should this info only be written to the CHLA/BBP profile ?
    ### TODO answer this and fill history information
    
    HISTORY_INSTITUTION = rep("XXXX", n_prof)
    ncvar_put(filenc_out, "HISTORY_INSTITUTION", HISTORY_INSTITUTION, start=c(1,1,i_history), count=c(4,n_prof,1))
    
    HISTORY_STEP = rep("YYYY", n_prof)
    ncvar_put(filenc_out, "HISTORY_STEP", HISTORY_STEP, start=c(1,1,i_history), count=c(4,n_prof,1))
    
    HISTORY_SOFTWARE = rep("ZZZZ", n_prof)
    ncvar_put(filenc_out, "HISTORY_SOFTWARE", HISTORY_SOFTWARE, start=c(1,1,i_history), count=c(4,n_prof,1))
    
    HISTORY_SOFTWARE_RELEASE = rep("0000", n_prof)
    ncvar_put(filenc_out, "HISTORY_SOFTWARE_RELEASE", HISTORY_SOFTWARE_RELEASE, start=c(1,1,i_history), count=c(4,n_prof,1))
    
    HISTORY_DATE = rep(DATE, n_prof)
    ncvar_put(filenc_out, "HISTORY_DATE", HISTORY_DATE, start=c(1,1,i_history), count=c(14,n_prof,1))
    
    HISTORY_ACTION = rep("AAAA", n_prof)
    ncvar_put(filenc_out, "HISTORY_ACTION", HISTORY_ACTION, start=c(1,1,i_history), count=c(4,n_prof,1))
    

   
    #nc_close(filenc_in)
    nc_close(filenc_out)
    
    return(0)
    
}

M = write_DM_MC(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST = DEEP_EST)
