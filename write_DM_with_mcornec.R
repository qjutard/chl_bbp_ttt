#########################
# This script allows to read BR files (later BD files) and to write a BD file with a delayed mode (DM)
# for CHLA and BBP700 according to the work done by M. Cornec
#
# The script is written as a function to ideally loop over the profiles of a float or parallelize eg with mcmapply
#########################

require(ncdf4) #deal with netcdf format files
require(oce) #calculate density sigma
require(MASS)
require(stringr)
require(parallel)
require(stringi)

source(paste(path_to_source, "process_files.R", sep=""))
source(paste(path_to_source, "write_DM.R", sep=""))
source(paste(path_to_source, "error_message.R", sep=""))
source(paste(path_to_source, "increment_N_CALIB.R", sep=""))

write_DM_MC <- function(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=NULL, index_greylist=NULL, 
                        accept_descent=FALSE, just_copy=FALSE, fill_value=FALSE, position_override=NULL, 
                        offset_override=NULL, accept_QC3=FALSE, only_BBP=FALSE, date_override=NULL,
                        only_CHL=FALSE){
    
	print(profile_actual)
	
    files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
    ident = strsplit(files,"/") #separate the different roots of the files paths
    ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
    prof_id = ident[,4] #retrieve all profiles  name as a vector
    
    iii = which(substr(prof_id,3,14)==profile_actual) #identify profile position in the index
    
    ############################
    ### Get the chla and bbp corrections from the meihod
    ############################
    
    CHLA_ADJUSTED = NULL
    BBP700_ADJUSTED = NULL
    CHLA_ADJUSTED_QC = NULL
    BBP700_ADJUSTED_QC = NULL
    CHLA_ADJUSTED_ERROR = NULL
    BBP700_ADJUSTED_ERROR = NULL
    
    chl_dark_offset = NA
    bbp_offset = NA
    chl_dark_min_pres = NA
    
    is_npq = FALSE
    npq_depth = NA
    npq_val = NA
	is_XB18 = FALSE
    
    if (!just_copy & !fill_value){
        
        L = try(process_file(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=DEEP_EST, index_greylist=index_greylist, 
                             accept_descent=accept_descent, accept_QC3=accept_QC3, position_override=position_override, 
                             offset_override=offset_override, date_override=date_override), silent=TRUE)
        
        if (inherits(L, "try-error")) {
            print("process_file(...) did not end properly")
            return(L)
        }
        if (is.numeric(L)) { # L is the error index of a managed error
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
		is_XB18 = L$is_XB18
    }

    ############################
    ### Create output file
    ############################
    
    path_split = unlist( strsplit(files[iii],"/") )
    path_to_profile = paste(path_split[1], path_split[2], path_split[3], sep="/")
    
    file_B = paste(path_to_netcdf, files[iii], sep="") # input file
    
    filenc_name_out = paste("BD",substring(path_split[4], 3),sep="")
    file_out = paste(path_to_netcdf, path_to_profile,"/DMMC/DMMC_profiles/", filenc_name_out, sep="") # output file
    dir_out = paste(path_to_netcdf, path_to_profile,"/DMMC/DMMC_profiles", sep="") # output directory
    
    # create directories if they do not exist
    system2("mkdir", c("-p", dir_out))
    
    if (length(file_B)!=1 | length(file_out)!=1) {
        print(error_message(206))
        return(206)
    }
    
    if (just_copy) {
        file_out_just_copy = unlist(strsplit(file_B,"/"))
        file_out_just_copy[length(file_out_just_copy)+1] = file_out_just_copy[length(file_out_just_copy)]
        file_out_just_copy[length(file_out_just_copy)-1] = "DMMC/DMMC_profiles"
        file_out_just_copy = paste(file_out_just_copy, collapse = "/")
        system2("cp", c(file_B, file_out_just_copy))
        return(0)
    }
    
    ### create the output file as a copy of the input
    system2("cp", c(file_B, file_out))
    
    
    ############################
    ### Create DATE and sctientific informations
    ############################
    
    date_update = Sys.time()
    DATE = stri_datetime_format(date_update, format="uuuuMMddHHmmss", tz="UTC")
    
    ### scientific_coefficient
    if (!is.na(chl_dark_offset)) {
        scientific_coefficient_chl = paste("CHLA_OFFSET =", round(chl_dark_offset,4))
    } else {
        scientific_coefficient_chl = paste("CHLA_OFFSET =", 0)
    }
    if (is_XB18) {
        scientific_coefficient_chl = paste(scientific_coefficient_chl, ", r = 0.092, iPARmid = 261, e = 2.2")
    }
    
    if (!is.na(bbp_offset)) {
        scientific_coefficient_bbp = paste("BBP700_OFFSET =", round(bbp_offset,4))
    } else {
        scientific_coefficient_bbp = ""
    }
    
    ### scientific equation
    if (is_XB18) {
        scientific_equation_chl = "CHLA_ADJUSTED = ((CHLA-CHLA_OFFSET)/2)/(r+(1-r)/(1+(DOWNWELLING_PAR/iPARmid)^e)), where DOWNWELLING_PAR is interpolated to align with CHLA"
    } else {
        scientific_equation_chl = "CHLA_ADJUSTED = (CHLA-CHLA_OFFSET)/2"
    }
    if (is_npq | !is.na(chl_dark_min_pres)){
        scientific_equation_chl = paste("otherwise ", scientific_equation_chl, sep="")
    }
    if (!is.na(chl_dark_min_pres)) {
        scientific_equation_chl = paste("CHLA_ADJUSTED = 0 for PRES in [", round(chl_dark_min_pres,4), ",+inf], ", scientific_equation_chl, sep="")
    }
    if (is_npq) {
        scientific_equation_chl = paste("CHLA_ADJUSTED = ", round(npq_val,4), " for PRES in [0, ", round(npq_depth,4),"], ", scientific_equation_chl, sep="")
    }
    if (is_npq & !is.na(chl_dark_min_pres)){
        if (npq_depth>chl_dark_min_pres){
            print(error_message(203))
            system2("rm", file_out)
            return(203)
        }
    }
    
    if (!is.na(bbp_offset)){
        scientific_equation_bbp = "BBP700_ADJUSTED = BBP700-BBP700_OFFSET"
    } else {
        scientific_equation_bbp = "BBP700_ADJUSTED = BBP700, no adjustment needed"
    }
    
    
    if (fill_value) {
        scientific_equation_chl = "CHLA_ADJUSTED = _FillValue"
        scientific_equation_bbp = "BBP700_ADJUSTED = _FillValue"
        scientific_coefficient_chl = ""
        scientific_coefficient_bbp = ""
    }
    
    ### scientific comment
    scientific_comment_chl = "CHLA delayed mode adjustment following the work done by M. Cornec in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)"
    scientific_comment_bbp = "BBP700 delayed mode adjustment following the work done by M. Cornec in Bellacicco et al. 2019 (http://dx.doi.org/10.1029/2019GL084078)"
    comment_dmqc_operator_PRIMARY = "PRIMARY | https://orcid.org/0000-0001-9992-5334 | Raphaelle Sauzede, CNRS" 
    comment_dmqc_operator_BBP700 = "BBP700 | https://orcid.org/0000-0002-1230-164X | Catherine Schmechtig, CNRS" 
    comment_dmqc_operator_CHLA = "CHLA | https://orcid.org/0000-0002-1230-164X | Catherine Schmechtig, CNRS" 
    #comment_dmqc_operator_PARAM = paste(param_name, " | https://orcid.org/0000-0002-1230-164X | Catherine Schmechtig, CNRS", sep="")
    
    if (!only_CHL) {
        exit = write_DM(file_out=file_out, param_name="BBP700", DATE=DATE, scientific_comment=scientific_comment_bbp, scientific_coefficient=scientific_coefficient_bbp, 
                 scientific_equation=scientific_equation_bbp, comment_dmqc_operator_PRIMARY=comment_dmqc_operator_PRIMARY, 
                 comment_dmqc_operator_PARAM=comment_dmqc_operator_BBP700, param_adjusted=BBP700_ADJUSTED, param_adjusted_qc=BBP700_ADJUSTED_QC, 
                 param_adjusted_error=BBP700_ADJUSTED_ERROR, fill_value=fill_value)
        if (exit!=0) {
            return(exit)
        }
    }
    
    if (!only_BBP) {
        exit = write_DM(file_out=file_out, param_name="CHLA", DATE=DATE, scientific_comment=scientific_comment_chl, scientific_coefficient=scientific_coefficient_chl, 
                 scientific_equation=scientific_equation_chl, comment_dmqc_operator_PRIMARY=comment_dmqc_operator_PRIMARY, 
                 comment_dmqc_operator_PARAM=comment_dmqc_operator_CHLA, param_adjusted=CHLA_ADJUSTED, param_adjusted_qc=CHLA_ADJUSTED_QC, 
                 param_adjusted_error=CHLA_ADJUSTED_ERROR, fill_value=fill_value)
        if (exit!=0) {
            return(exit)
        }
    }

    return(0)
    
}

