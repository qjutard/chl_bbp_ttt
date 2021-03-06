######################################################################
# Script to be called from bash by DMMC.sh
# Starts the Delayed Mode with the adequate options
######################################################################

uf = commandArgs(trailingOnly = TRUE)

print(uf)

### get all arguments
profile_WMO = uf[1]
List = uf[2]
DEEP_EST_table = uf[3]
just_copy = as.logical(uf[4])
fill_value = as.logical(uf[5])
accept_descent = as.logical(uf[6])
accept_QC3 = as.logical(uf[7])
Profile = uf[8]
position_override_call = uf[9]
offset_override_call = uf[10]
only_BBP = as.logical(uf[11])
date_override_call = uf[12]
offset_override_file = uf[13]
only_CHL = as.logical(uf[14])
test_env = as.logical(uf[15])
multi_core = uf[16]

### Check conflicting options
exists_WMO = as.numeric(profile_WMO!="NA")
exists_List = as.numeric(List!="NA")
exists_Profile = as.numeric(Profile!="NA")

if ( (exists_WMO + exists_List + exists_Profile)!=1 ) {
	print("Please give exactly one of 'WMO' (-W), 'List' (-L), or 'Profile' (-P) as argument")
	stop()
}

if (offset_override_call!="NA" & offset_override_file!="NA") {
    print("Cannot override offsets with -o and -O at the same time")
    stop()
}

if (only_BBP & only_CHL) {
    print("-B and -C are conflicting, use neither if you want to write delayed mode of both variables")
    stop()
}


### import pathways
source("~/Documents/cornec_chla_qc/chl_bbp_ttt/pathways.R")
if (test_env) {
	path_to_netcdf = path_to_netcdf_test_env
}

### Source the DM functions and subfunctions, and libraries
library(parallel)
library(ncdf4)
library(stringr)
library(stringi)
library(oce)
library(MASS)

source(paste(path_to_source, "DM_writing/write_DM_with_mcornec.R", sep=""))
source(paste(path_to_source, "DM_writing/process_files.R", sep=""))
source(paste(path_to_source, "DM_writing/write_DM.R", sep=""))
source(paste(path_to_source, "DM_writing/error_message.R", sep=""))
source(paste(path_to_source, "DM_writing/increment_N_CALIB.R", sep=""))

dir_function = paste(path_to_source,"Functions/", sep="")
source(paste(dir_function,"NPQ_cor_X12_XB18.R",sep="")) # Needed to correct the profile from the NPQ (in presence of PAR measured in situ)
source(paste(dir_function,"NPQ_cor_P18.R",sep="")) # Needed to correct the profile from the NPQ (in absence of PAR measured in situ)
source(paste(dir_function,"RunningFilter.R",sep="")) # Needed to smooth some profiles (running median/mean on a given window)
source(paste(dir_function,"MLD_calc.R",sep="")) #Needed to calculate the MLD
source(paste(dir_function,"Outliars_med.R",sep="")) #Needed to test outliars in the Dark offset time series
source(paste(dir_function,"Dark_MLD_table_coriolis.R",sep="")) #Needed to correct the dark offset in cases of deep vertical mixing
source(paste(dir_function,"Zone.R",sep="")) # Needed to attribute a regional criterion on the profile for the cases where a regional correction is needed
source(paste(dir_function,"Darkoz.R",sep="")) # Needed to calculate a dark offset in the regions of Oxygen Minimum Zone
source(paste(dir_function,"DarkXing.R",sep="")) # Needed to calculate a dark offset in the regions where the chla show an increase at depth 
source(paste(dir_function,"Dark_Fchla_Corr.R",sep="")) # Needed to correct the dark offset of the chla (depends on the float location, and if there are cases of deep vertical mixing)


### import tables
index_ifremer = read.table(path_to_index_ifremer, sep=",", header = T)
index_greylist = read.csv(path_to_index_greylist, sep = ",")

### build profile_list from WMO, list, or profile name
if (profile_WMO!="NA") {
    files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
    ident = strsplit(files,"/") #separate the different roots of the files paths
    ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
    prof_id = ident[,4] #retrieve all profiles  name as a vector
    prof_id_WMO = substr(prof_id, 3, 9)
    profile_list_all = substr(prof_id[which(prof_id_WMO==profile_WMO)], 3, 14)
} else if (List!="NA") {
    profile_list_all = read.table(List, colClasses="character")$V1
} else {
    profile_list_all = c(Profile)
}



### Treat optional arguments

if (multi_core == "NA") {
    num_cores = detectCores()
} else {
    num_cores = as.numeric(multi_cores)
}

if (position_override_call == "NA") {
	position_override = NULL
} else {
	position_override = as.numeric(unlist(strsplit(position_override_call, ";")))
}

if (offset_override_call == "NA") {
	offset_override = list(NULL)
} else if (offset_override_call=="dmmc") {
    offset_override = list(offset_override_call)
} else {
	offset_override = list(as.numeric(unlist(strsplit(offset_override_call, ";"))))
}
if (offset_override_file != "NA") {
    offset_override = NULL 
    offset_table = read.table(offset_override_file)
    for (i in 1:length(profile_list_all)) {
        id_prof = which(offset_table$V1 == profile_list_all[i])
        if (length(id_prof)!=1) {
            print("One or more profiles has a non unique or non existing corresponding offset in the given offset file (-O)")
            stop()
        }
        offset_override[[i]] = c(offset_table$V2[i], NA) # currently the min argument can only be NA
    }
}

if (DEEP_EST_table == "NA") { # calculate deep est if it is not given
    if (!just_copy & !fill_value & offset_override_call=="NA" & offset_override_file=="NA") { # DEEP_EST is not necessary if we just want to copy or fill the files, or if offset values are given
        profile_actual = profile_list_all[1]
        DEEP_EST = Dark_MLD_table_coriolis(substr(profile_actual,1,7), path_to_netcdf, index_ifremer, n_cores=num_cores)
        write.table(DEEP_EST, "DEEP_EST.t", row.names=F)
    } else {
        DEEP_EST = NULL
    }
} else {
    DEEP_EST = try(read.table(DEEP_EST_table, header=T), silent=TRUE)
    if (inherits(DEEP_EST, "try-error")) {
        print(paste("Warning :", DEEP_EST_table, "could not be opened as a table"))
        DEEP_EST = NULL
    }
}

if (date_override_call == "NA") {
    date_override = NULL
} else {
    date_override = date_override_call
}



### Compute and write delayed modes
M = mcmapply(write_DM_MC, profile_actual=profile_list_all, offset_override=offset_override, mc.cores=num_cores,
             MoreArgs=list(index_ifremer=index_ifremer, path_to_netcdf=path_to_netcdf, DEEP_EST=DEEP_EST, index_greylist=index_greylist, accept_descent=accept_descent,
                           just_copy=just_copy, fill_value=fill_value, accept_QC3=accept_QC3, position_override=position_override, only_BBP=only_BBP, 
                           date_override=date_override, only_CHL=only_CHL))

### assess error messages
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

write.table(unlist(messages), "list_errors.t", col.names=FALSE)

if (test_env) {
	write.table(c(n_profiles, n_success, n_fails, n_errors), "test_n_errors.t", col.names=FALSE, row.names=FALSE)
}
