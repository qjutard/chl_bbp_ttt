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
library(parallel)
library(stringi)

source("process_files.R")
source("error_message.R")

write_DM_MC <- function(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=NULL, just_copy=FALSE){
    
    files<-as.character(index_ifremer[,1]) #retrieve the path of each netcfd file
    ident<-strsplit(files,"/") #separate the different roots of the files paths
    ident<-matrix(unlist(ident), ncol=4, byrow=TRUE)
    prof_id<-ident[,4] #retrieve all profiles  name as a vector
    
    i <-which(substr(prof_id,3,14)==profile_actual) #identify profile position in the index
    
    ############################
    ### Get the chla and bbp corrections from the method
    ############################
    
    L = try(process_file(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=DEEP_EST), silent=TRUE)
    if (!is.list(L)){
        print("process_file(...) did not end properly")
        return(L)
    }
    CHLA_ADJUSTED = L$CHLA_ADJUSTED
    BBP700_ADJUSTED = L$BBP700_ADJUSTED
    CHLA_ADJUSTED_QC = L$CHLA_ADJUSTED_QC
    BBP700_ADJUSTED_QC = L$BBP700_ADJUSTED_QC
    CHLA_ADJUSTED_ERROR = L$CHLA_ADJUSTED_ERROR
    BBP700_ADJUSTED_ERROR = L$BBP700_ADJUSTED_ERROR
    chl_dark_offset = L$chl_dark_offset
    chl_dark_min_pres = L$chl_dark_min_pres
    bbp_offset = L$bbp_offset
    is_npq = L$is_npq
    npq_depth = L$npq_depth
    npq_val = L$npq_val

    ############################
    ### Open input and output files
    ############################
    
    path_split = unlist( strsplit(files[i],"/") )
    path_to_profile = paste(path_split[1], path_split[2], path_split[3], sep="/")
    
    filenc_name_M = path_split[4]
    filenc_name_B = paste("B?",substring(filenc_name_M, 3),sep="")
    filenc_name_out = paste("BD",substring(filenc_name_M, 3,nchar(filenc_name_M)-3),".nc",sep="")
    
    file_B = paste(path_to_netcdf, path_to_profile,"/", filenc_name_B, sep="") 
    file_B = system2("ls",file_B,stdout=TRUE) # identify R or D file 
    file_out = paste(path_to_netcdf, path_to_profile,"/DM_cornec/", filenc_name_out, sep="") 
    
    if (length(file_B)!=1 | length(file_out)!=1) {
        print(error_message(206))
        return(206)
    }
    
    ### check whether the existing file is a D file
    path_split_B = unlist( strsplit(file_B, "/") )
    B_type = substr(path_split_B[length(path_split_B)],1,2)
    if (B_type=="BD"){
        print("a BD-file has been given as an input, these are currently not accepted")
        return(201)
    }
    
    ### create the output file as a copy of the input
    system2("cp", c(file_B, file_out))
    if (just_copy) {
        return(0)
    }
    
    #filenc_in <- nc_open(file_B, readunlim=FALSE, write=FALSE)
    filenc_out <- nc_open(file_out, readunlim=FALSE, write=TRUE)
    
    ############################
    ### Find some parameter indices and define date
    ############################
    
    parameters=ncvar_get(filenc_out,"STATION_PARAMETERS")
    
    id_param_chl = grep("CHLA                                                            ",parameters)
    id_param_bbp = grep("BBP700                                                          ",parameters)
    id_prof = which(parameters=="CHLA                                                            ", arr.ind=TRUE)[2] # The CHL and BBP share the same profile
    id_param_chl_arr = which(parameters=="CHLA                                                            ", arr.ind=TRUE)
    id_param_bbp_arr = which(parameters=="BBP700                                                          ", arr.ind=TRUE)
    n_prof = dim(parameters)[2]
    
    if ( length(id_param_chl)!=1 | length(id_param_bbp)!=1 ){
        print("several profiles of chl or bbp detected")
        nc_close(filenc_out)
        system2("rm", file_out)
        return(202)
    } 
    
    N_HISTORY=filenc_out$dim[['N_HISTORY']]$len
    i_history=N_HISTORY+1
    
    date_update = Sys.time()
    
    DATE = stri_datetime_format(date_update, format="uuuuMMddHHmmss", tz="UTC")
    
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
    
    ### scientific_coefficient
    if (!is.na(chl_dark_offset)) {
        scientific_coefficient_chl = paste("CHLA_OFFSET =", chl_dark_offset)
    } else {
        scientific_coefficient_chl = ""
    }
    
    if (!is.na(bbp_offset)) {
        scientific_coefficient_bbp = paste("BBP700_OFFSET =", bbp_offset)
    } else {
        scientific_coefficient_bbp = ""
    }
    
    ### scientific equation
    if (!is.na(chl_dark_offset)) {
        scientific_equation_chl = "CHLA_ADJUSTED = (CHLA-CHLA_OFFSET)/2"
        if (!is.na(chl_dark_min_pres)) {
            scientific_equation_chl = paste("CHLA_ADJUSTED = 0 for PRES in [", chl_dark_min_pres, ",+inf], ", scientific_equation_chl, sep="")
        }
    } else {
        scientific_equation_chl = "CHLA_ADJUSTED = CHLA/2"
    }
    if (is_npq) {
        scientific_equation_chl = paste("CHLA_ADJUSTED = ", npq_val, " for PRES in [0, ", npq_depth,"], ", scientific_equation_chl, sep="")
    }
    if (is_npq | !is.na(chl_dark_min_pres)){
        scientific_equation_chl = paste(scientific_equation_chl, " otherwise", sep="")
    }
    if (is_npq & !is.na(chl_dark_min_pres)){
        if (npq_depth>chl_dark_min_pres){
            print(error_message(203))
            nc_close(filenc_out)
            system2("rm", file_out)
            return(203)
        }
    }
    
    if (!is.na(bbp_offset)){
        scientific_equation_bbp = "BBP700_ADJUSTED = BBP700-BBP700_OFFSET"
    } else {
        scientific_equation_bbp = "BBP700_ADJUSTED = BBP700"
    }
    
    ### scientific comment
    scientific_comment_chl = "CHLA delayed mode adjustment following the work done by M. Cornec in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)"
    scientific_comment_bbp = "BBP700 delayed mode adjustment following the work done by M. Cornec in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)"
    
    ### calib date
    scientific_date_chl = DATE
    scientific_date_bbp = DATE
    
    
    ### write on file
    scientific_calib_chl = c(scientific_comment_chl, scientific_coefficient_chl, scientific_equation_chl, scientific_date_chl)
    scientific_calib_bbp = c(scientific_comment_bbp, scientific_coefficient_bbp, scientific_equation_bbp, scientific_date_bbp)
    SCIENTIFIC_CALIB_VARIABLE = paste("SCIENTIFIC_CALIB_", c("COMMENT", "COEFFICIENT", "EQUATION", "DATE"), sep="")
    
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

    HISTORY_INSTITUTION = "VF  "
    ncvar_put(filenc_out, "HISTORY_INSTITUTION", HISTORY_INSTITUTION, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_STEP = "ARSQ"
    ncvar_put(filenc_out, "HISTORY_STEP", HISTORY_STEP, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_SOFTWARE = rep("DMMC", n_prof)
    ncvar_put(filenc_out, "HISTORY_SOFTWARE", HISTORY_SOFTWARE, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_SOFTWARE_RELEASE = rep("0000", n_prof)
    ncvar_put(filenc_out, "HISTORY_SOFTWARE_RELEASE", HISTORY_SOFTWARE_RELEASE, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_DATE = rep(DATE, n_prof)
    ncvar_put(filenc_out, "HISTORY_DATE", HISTORY_DATE, start=c(1,id_prof,i_history), count=c(14,1,1))
    
    HISTORY_ACTION = rep("CV  ", n_prof)
    ncvar_put(filenc_out, "HISTORY_ACTION", HISTORY_ACTION, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    ############################
    ### Write profile_QC
    ############################
    
    chl_QC = unlist(strsplit(CHLA_ADJUSTED_QC[id_prof],""))
    bbp_QC = unlist(strsplit(BBP700_ADJUSTED_QC[id_prof],""))

    n_QC_chl = sum( chl_QC!=" " )
    n_QC_bbp = sum( bbp_QC!=" " )
    
    n_good_chl = 100 * sum( chl_QC=="1" | chl_QC=="2" | chl_QC=="5" | chl_QC=="8" ) / n_QC_chl
    n_good_bbp = 100 * sum( bbp_QC=="1" | bbp_QC=="2" | bbp_QC=="5" | bbp_QC=="8" ) / n_QC_bbp
    
    # Write CHLA profile QC
    if ( n_good_chl == 0) ncvar_put(filenc_out,"PROFILE_CHLA_QC","F",start=id_prof, count=1)
    if ( n_good_chl > 0 && n_good_chl < 25 ) ncvar_put(filenc_out,"PROFILE_CHLA_QC","E",start=id_prof, count=1)
    if ( n_good_chl >= 25 && n_good_chl < 50 ) ncvar_put(filenc_out,"PROFILE_CHLA_QC","D",start=id_prof, count=1)
    if ( n_good_chl >= 50 && n_good_chl < 75 ) ncvar_put(filenc_out,"PROFILE_CHLA_QC","C",start=id_prof, count=1)
    if ( n_good_chl >= 75 && n_good_chl < 100 ) ncvar_put(filenc_out,"PROFILE_CHLA_QC","B",start=id_prof, count=1)
    if ( n_good_chl == 100 ) ncvar_put(filenc_out,"PROFILE_CHLA_QC","A",start=id_prof, count=1)
    
    # Write BBP700 profile QC
    if ( n_good_bbp == 0) ncvar_put(filenc_out,"PROFILE_BBP700_QC","F",start=id_prof, count=1)
    if ( n_good_bbp > 0 && n_good_bbp < 25 ) ncvar_put(filenc_out,"PROFILE_BBP700_QC","E",start=id_prof, count=1)
    if ( n_good_bbp >= 25 && n_good_bbp < 50 ) ncvar_put(filenc_out,"PROFILE_BBP700_QC","D",start=id_prof, count=1)
    if ( n_good_bbp >= 50 && n_good_bbp < 75 ) ncvar_put(filenc_out,"PROFILE_BBP700_QC","C",start=id_prof, count=1)
    if ( n_good_bbp >= 75 && n_good_bbp < 100 ) ncvar_put(filenc_out,"PROFILE_BBP700_QC","B",start=id_prof, count=1)
    if ( n_good_bbp == 100 ) ncvar_put(filenc_out,"PROFILE_BBP700_QC","A",start=id_prof, count=1)
        
    ############################
    ### Write other variables
    ############################
    
    # DATA_MODE
    ncvar_put(filenc_out, "DATA_MODE", "D", start=c(id_prof), count=c(1))
    
    # DATA_STATE_INDICATOR
    DATA_STATE_INDICATOR = rep("2C  ", n_prof)
    ncvar_put(filenc_out, "DATA_STATE_INDICATOR", DATA_STATE_INDICATOR, start=c(1,1), count=c(4,n_prof))
    
    # PARAMETER_DATA_MODE
    PARAMETER_DATA_MODE = ncvar_get(filenc_out,"PARAMETER_DATA_MODE")
    str_sub(PARAMETER_DATA_MODE[id_param_chl_arr[2]], id_param_chl_arr[1],id_param_chl_arr[1]) = "D"
    str_sub(PARAMETER_DATA_MODE[id_param_bbp_arr[2]], id_param_bbp_arr[1],id_param_bbp_arr[1]) = "D"
    ncvar_put(filenc_out, "PARAMETER_DATA_MODE", PARAMETER_DATA_MODE)
    
    # DATE_UPDATE
    ncvar_put(filenc_out, "DATE_UPDATE", DATE)
    
    ############################
    ### Change attributes
    ############################
    
    comment_dmqc_operator1 = "PRIMARY | https://orcid.org/16-digit-number | operator name, institution" ### TODO fill
    comment_dmqc_operator2 = "CHLA | https://orcid.org/16-digit-number | operator name, institution" ### TODO fill
    comment_dmqc_operator3 = "BBP700 | https://orcid.org/16-digit-number | operator name, institution" ### TODO fill
    ncatt_put(filenc_out, varid=0, "comment_dmqc_operator1", comment_dmqc_operator1)
    ncatt_put(filenc_out, varid=0, "comment_dmqc_operator2", comment_dmqc_operator2)
    ncatt_put(filenc_out, varid=0, "comment_dmqc_operator3", comment_dmqc_operator3)
    
    #TODO change history attribute, is it necessary ?
    history = ncatt_get(filenc_out, varid=0, "history")$value
    
    DATE_history_day = stri_datetime_format(date_update, format="uuuu-MM-dd", tz="UTC")
    DATE_history_hour = stri_datetime_format(date_update, format="HH:mm:ss", tz="UTC")
    
    history_last_update = paste(DATE_history_day, "T", DATE_history_hour, "Z", " last update (personal code)", sep="")
    history_new = paste(unlist(strsplit(history,";"))[1], "; ", history_last_update, sep="")
    
    ncatt_put(filenc_out, varid=0, "history", history_new)
    
    
   
    #nc_close(filenc_in)
    nc_close(filenc_out)
    
    return(0)
    
}

### Paths and profile definition
index_ifremer<-read.table("~/Documents/data/argo_merge-profile_index.txt", skip=9, sep = ",")
path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"
#profile_WMO = "6901524"
profile_WMO = "6901524"


### build list of profiles from float WMO
files<-as.character(index_ifremer[,1]) #retrieve the path of each netcfd file
ident<-strsplit(files,"/") #separate the different roots of the files paths
ident<-matrix(unlist(ident), ncol=4, byrow=TRUE)
prof_id<-ident[,4] #retrieve all profiles  name as a vector
prof_id_WMO = substr(prof_id, 3, 9)
profile_list_all = substr(prof_id[which(prof_id_WMO==profile_WMO)], 3, 14)
profile_actual = profile_list_all[1]

# Test values
#profile_list<-c("6901524_150.","6901524_151.","6901524_152.","6901524_153.","6901524_154.","6901524_155.","6901524_001D")
#profile_actual = profile_list[1]
#profile_actual = "6901524_233."
#profile_actual = "6901524_087."
#profile_actual = "6901527_213."
#profile_actual = "6901527_277."

### DEEP_EST should be computed once per FLOAT 
DEEP_EST = Dark_MLD_table_coriolis(substr(profile_actual,1,7), path_to_netcdf, index_ifremer) 

# Test the func
#M = write_DM_MC(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST = DEEP_EST)

numCores = detectCores()
M = mcmapply(write_DM_MC, profile_list_all, MoreArgs=list(index_ifremer, path_to_netcdf, DEEP_EST = DEEP_EST), mc.cores=numCores)

errors = as.numeric(M)

n_profiles = length(profile_list_all)
n_success = sum(errors==0, na.rm=TRUE) 
n_errors = sum(is.na(errors))
n_fails = n_profiles - n_success - n_errors

print(paste(n_profiles,"profiles were found and treated, of which",n_success,"succesfully finished,",n_fails,"were stopped, and",n_errors,"returned an unmanaged error"))

is_error = which(errors!=0 | is.na(errors))

messages = M[is_error]
is_managed_error = which(!is.na(errors[is_error]))
messages[is_managed_error] = lapply(messages[is_managed_error], error_message)
View(messages)