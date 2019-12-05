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

### Check that exactly one of WMO, List, or Profile is given
exists_WMO = as.numeric(profile_WMO!="NA")
exists_List = as.numeric(List!="NA")
exists_Profile = as.numeric(Profile!="NA")

if ( (exists_WMO + exists_List + exists_Profile)!=1 ) {
	print("Please give exactly one of 'WMO' (-W), 'List' (-L), or 'Profile' (-P) as argument")
	stop()
}

### import pathways
source("~/Documents/cornec_chla_qc/chl_bbp_ttt/pathways.R")
### Source the DM function and subfunctions
source("~/Documents/cornec_chla_qc/chl_bbp_ttt/write_DM_with_mcornec.R")

### import tables
index_ifremer<-read.table(path_to_index_ifremer, skip=9, sep = ",")
index_greylist<-read.csv(path_to_index_greylist, sep = ",")

### build profile_list from WMO
if (profile_WMO!="NA") {
    files<-as.character(index_ifremer[,1]) #retrieve the path of each netcfd file
    ident<-strsplit(files,"/") #separate the different roots of the files paths
    ident<-matrix(unlist(ident), ncol=4, byrow=TRUE)
    prof_id<-ident[,4] #retrieve all profiles  name as a vector
    prof_id_WMO = substr(prof_id, 3, 9)
    profile_list_all = substr(prof_id[which(prof_id_WMO==profile_WMO)], 3, 14)
} else if (List!="NA") {
    profile_list_all = read.table(List, colClasses="character")$V1
} else {
    profile_list_all = c(Profile)
}

if (DEEP_EST_table=="NA") { # calculate deep est if it is not given
	if (!just_copy & !fill_value) { # DEEP_EST is not necessary if we just want to copy or fill the files
    	profile_actual = profile_list_all[1]
    	DEEP_EST = Dark_MLD_table_coriolis(substr(profile_actual,1,7), path_to_netcdf, index_ifremer)
    	write.table(DEEP_EST, "DEEP_EST.t")
	} else {
		DEEP_EST = NULL
	}
} else {
    DEEP_EST = try(read.table(DEEP_EST_table), silent = TRUE)
    if (inherits(DEEP_EST, "try-error")) {
        DEEP_EST = NULL
    }
}

### Compute and write delayed modes
numCores = detectCores()
M = mcmapply(write_DM_MC, profile_list_all, MoreArgs=list(index_ifremer, path_to_netcdf, DEEP_EST = DEEP_EST, index_greylist=index_greylist, 
                                                          accept_descent=accept_descent, just_copy=just_copy, fill_value=fill_value, 
                                                          accept_QC3=accept_QC3), mc.cores=numCores)

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

write.table(unlist(messages), "list_errors.t", col.names = FALSE)
