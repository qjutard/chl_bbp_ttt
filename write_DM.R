################################################################
# Function to write DM information of a parameter onto a B file
################################################################

write_DM <- function(file_out, param_name, DATE, scientific_comment, scientific_coefficient, scientific_equation, 
                     param_adjusted=NULL, param_adjusted_qc=NULL, param_adjusted_error=NULL, fill_value=FALSE) {
    
    filenc_out <- nc_open(file_out, readunlim=FALSE, write=TRUE)
    
    ############################
    ### Find some parameter indices
    ############################
    
    parameters = ncvar_get(filenc_out,"STATION_PARAMETERS")
    
    id_param = grep(str_pad(param_name, 64, side="right"), parameters)
    id_prof = which(parameters==str_pad(param_name, 64, side="right"), arr.ind=TRUE)[2]
    id_param_arr = which(parameters==str_pad(param_name, 64, side="right"), arr.ind=TRUE)
    n_prof = dim(parameters)[2]
    
    if ( length(id_param)<1 ){ 
        print(error_message(205))
        nc_close(filenc_out)
        system2("rm", file_out)
        return(205)
    } 
    
    if ( length(id_param)>1 ){ 
        print(error_message(202))
        nc_close(filenc_out)
        system2("rm", file_out)
        return(202)
    } 
    
    N_HISTORY = filenc_out$dim[['N_HISTORY']]$len
    i_history = N_HISTORY+1
    
    
    ### find adequate N_CALIB
    calib_date = ncvar_get(filenc_out, "SCIENTIFIC_CALIB_DATE")
    
    N_PARAM = filenc_out$dim[['N_PARAM']]$len
    N_CALIB = filenc_out$dim[['N_CALIB']]$len
    N_PROF = filenc_out$dim[['N_PROF']]$len
    
    # check that dimensions are aligned correctly
    dim_par = dim(parameters)
    dim_cal = dim(calib_date)
    if ( dim_par[1]!=N_PARAM | dim_par[2]!=N_PROF | dim_cal[1]!=N_PARAM | dim_cal[2]!=N_CALIB | dim_cal[3]!=N_PROF ){
        print(error_message(204))
        return(204)
    }
    
    # get vectors of dates corresponding to the parameter
    calib_date = calib_date[id_param_arr[1],,id_param_arr[2]]
    
    # find the n_calib of the latest calibration
    last_cal = max(which(calib_date!="              "),0)
    
    # increment n_calib for the parameter
    new_cal = last_cal + 1
    
    # position to write calibration to
    id_calib = c(id_param_arr[1], new_cal, id_param_arr[2])
    
    ############################
    ### Increment N_CALIB if necessary
    ############################
    
    if (new_cal > N_CALIB) {
        
        nc_close(filenc_out)
        
        inc_ret = increment_N_CALIB(file_out=file_out, file_out_copy=paste(file_out, "_copy", sep=""))
        if (inc_ret!=0) { 
            system2("rm", file_out)
            return(inc_ret) 
        }
        
        filenc_out = nc_open(file_out, readunlim=FALSE, write=TRUE)
        
        # Write the new level of PARAMETER
        ncvar_put(filenc_out, "PARAMETER", str_pad(param_name, 64, side="right"), start=c(1, id_calib), count=c(64,1,1,1))
        
    }
    
    
    ############################
    ### Write correction from method to BD-file
    ############################
    
    PARAM_ADJUSTED_NAME = paste(param_name, "_ADJUSTED", sep="")
    PARAM_ADJUSTED_QC_NAME = paste(param_name, "_ADJUSTED_QC", sep="")
    PARAM_ADJUSTED_ERROR_NAME = paste(param_name, "_ADJUSTED_ERROR", sep="")
    
    if (!fill_value) { 
        ncvar_put(filenc_out, PARAM_ADJUSTED_NAME, param_adjusted)
        ncvar_put(filenc_out, PARAM_ADJUSTED_QC_NAME, param_adjusted_qc)
        ncvar_put(filenc_out, PARAM_ADJUSTED_ERROR_NAME, param_adjusted_error)
    } else {
        
        fill_params = c(PARAM_ADJUSTED_NAME, PARAM_ADJUSTED_ERROR_NAME)
        for (fill_name in fill_params) {
            fill_var = ncvar_get(filenc_out, fill_name)
            fill_var[] = NA
            ncvar_put(filenc_out, fill_name, fill_var)
        }

        fill_test = ncvar_get(filenc_out, param_name, start=c(1,id_prof)) #get the original profile
        fill_space = which(is.na(fill_test))
        fill_test[] = "4"
        fill_test[fill_space] = " "
        fill_test = paste(fill_test, collapse="")
        ncvar_put(filenc_out, PARAM_ADJUSTED_QC_NAME, fill_test, start=c(1,id_prof))

    }
    
    ############################
    ### Write scientific_calib
    ############################
    
    ### calib date
    scientific_date = DATE
    
    ### write on file
    scientific_calib = c(scientific_comment, scientific_coefficient, scientific_equation, scientific_date)
    SCIENTIFIC_CALIB_VARIABLE = paste("SCIENTIFIC_CALIB_", c("COMMENT", "COEFFICIENT", "EQUATION", "DATE"), sep="")
    
    for (i in seq(1,4)){
        
        SCIENTIFIC_CALIB_INFO = str_pad(scientific_calib[i], 256, "right")
        
        scientific_calib_info = ncvar_get(filenc_out, SCIENTIFIC_CALIB_VARIABLE[i])
        
        scientific_calib_info[id_calib[1], id_calib[2], id_calib[3]] = SCIENTIFIC_CALIB_INFO
        
        ncvar_put(filenc_out , SCIENTIFIC_CALIB_VARIABLE[i], scientific_calib_info)
        
    }
    
    ############################
    ### Write history
    ############################
    
    HISTORY_INSTITUTION = "VF  "
    ncvar_put(filenc_out, "HISTORY_INSTITUTION", HISTORY_INSTITUTION, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_STEP = "ARSQ"
    ncvar_put(filenc_out, "HISTORY_STEP", HISTORY_STEP, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_SOFTWARE = "DMMC"
    ncvar_put(filenc_out, "HISTORY_SOFTWARE", HISTORY_SOFTWARE, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_SOFTWARE_RELEASE = "1.00"
    ncvar_put(filenc_out, "HISTORY_SOFTWARE_RELEASE", HISTORY_SOFTWARE_RELEASE, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_DATE = DATE
    ncvar_put(filenc_out, "HISTORY_DATE", HISTORY_DATE, start=c(1,id_prof,i_history), count=c(14,1,1))
    
    HISTORY_ACTION = "CV  "
    ncvar_put(filenc_out, "HISTORY_ACTION", HISTORY_ACTION, start=c(1,id_prof,i_history), count=c(4,1,1))
    
    HISTORY_PARAMETER = str_pad(param_name, 64, side="right")
    ncvar_put(filenc_out, "HISTORY_PARAMETER", HISTORY_PARAMETER, start=c(1,id_prof,i_history), count=c(64,1,1))
    
    ############################
    ### Write profile_QC
    ############################
    
    if (!fill_value) {
        axis_QC = unlist(strsplit(param_adjusted_qc[id_prof],""))
    } else {
        axis_QC = "4"
    }
    
    n_QC = sum( axis_QC!=" " )
    n_good = 100 * sum( axis_QC=="1" | axis_QC=="2" | axis_QC=="5" | axis_QC=="8" ) / n_QC
    
    PROFILE_PARAM_QC_NAME = paste("PROFILE_", param_name, "_QC", sep="")

    # Write BBP700 profile QC
    if ( n_good == 0) ncvar_put(filenc_out, PROFILE_PARAM_QC_NAME, "F", start=id_prof, count=1)
    if ( n_good > 0 && n_good < 25 ) ncvar_put(filenc_out, PROFILE_PARAM_QC_NAME, "E", start=id_prof, count=1)
    if ( n_good >= 25 && n_good < 50 ) ncvar_put(filenc_out, PROFILE_PARAM_QC_NAME, "D", start=id_prof, count=1)
    if ( n_good >= 50 && n_good < 75 ) ncvar_put(filenc_out, PROFILE_PARAM_QC_NAME, "C", start=id_prof, count=1)
    if ( n_good >= 75 && n_good < 100 ) ncvar_put(filenc_out, PROFILE_PARAM_QC_NAME, "B", start=id_prof, count=1)
    if ( n_good == 100 ) ncvar_put(filenc_out, PROFILE_PARAM_QC_NAME, "A", start=id_prof, count=1)
    
    ############################
    ### Write other variables
    ############################
    
    # DATA_MODE
    ncvar_put(filenc_out, "DATA_MODE", "D", start=c(id_prof), count=c(1))
    
    # DATA_STATE_INDICATOR
    DATA_STATE_INDICATOR = "2C  "
    ncvar_put(filenc_out, "DATA_STATE_INDICATOR", DATA_STATE_INDICATOR, start=c(1,id_prof), count=c(4,1))
    
    # PARAMETER_DATA_MODE
    PARAMETER_DATA_MODE = ncvar_get(filenc_out,"PARAMETER_DATA_MODE")
    str_sub(PARAMETER_DATA_MODE[id_param_arr[2]], id_param_arr[1],id_param_arr[1]) = "D"
    ncvar_put(filenc_out, "PARAMETER_DATA_MODE", PARAMETER_DATA_MODE)
    
    # DATE_UPDATE
    ncvar_put(filenc_out, "DATE_UPDATE", DATE)
    
    ############################
    ### Change attributes
    ############################
    
    comment_dmqc_operator_PRIMARY = "PRIMARY | https://orcid.org/0000-0001-9992-5334 | Raphaelle Sauzede, CNRS" 
    #comment_dmqc_operator_BBP700 = "BBP700 | https://orcid.org/0000-0002-1230-164X | Catherine Schmechtig, CNRS" 
    #comment_dmqc_operator_CHLA = "CHLA | https://orcid.org/0000-0002-1230-164X | Catherine Schmechtig, CNRS" 
    comment_dmqc_operator_PARAM = paste(param_name, " | https://orcid.org/0000-0002-1230-164X | Catherine Schmechtig, CNRS", sep="")
    
    all_att = ncatt_get(filenc_out, varid=0)
    
    if (is.null(all_att[["comment_dmqc_operator1"]])) {
        ncatt_put(filenc_out, varid=0, "comment_dmqc_operator1", comment_dmqc_operator_PRIMARY)
    }
    
    for (i in 2:length(all_att)) { #find the first empty comment_dmqc_operatorX (X will always be lower than length(all_att))
        att_name = paste("comment_dmqc_operator", i, sep="")
        if ( is.null(all_att[[att_name]]) ) {
            ncatt_put(filenc_out, varid=0, att_name, comment_dmqc_operator_PARAM)
            break
        }
    }
    
    nc_close(filenc_out)
    
    return(0)
}