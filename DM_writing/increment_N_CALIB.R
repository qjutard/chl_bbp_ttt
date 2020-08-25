##############################################################################################
# This function takes a netcdf file name (with its path) and adds 1 to its N_CALIB dimension
# by writing a temporary file of which the name is also an input
##############################################################################################

require(ncdf4)

increment_N_CALIB <- function(file_out, file_out_copy) {
    
    if (length(file_out)!=1 | length(file_out_copy)!=1) {
        return(306)
    }
    
    filenc_out = nc_open(file_out, readunlim=FALSE, write=TRUE)
    
    ### This requires creating a new file and copying everything to it with a new N_calib
    
    len_old = filenc_out$dim$N_CALIB$len
    vals_old = filenc_out$dim$N_CALIB$vals
    id_ncalib = filenc_out$dim$N_CALIB$id
    len_new = len_old + 1
    vals_new = 1:len_new
    
    old_var = filenc_out$var
    new_var = filenc_out$var
    
    ### create new variables with updated N_CALIB
    for (j in seq(1, filenc_out$nvars)) { # iterate on variables
        for (k in seq(1, old_var[[j]]$ndims)) { # iterate on dimensions in variable
            if (old_var[[j]]$dim[[k]]$id == id_ncalib) { # id the dimension id corresponds to N_CALIB
                new_var[[j]]$dim[[k]]$len = len_new
                new_var[[j]]$dim[[k]]$vals = vals_new
            }
        }
    }
    
    filenc_copy = nc_create(file_out_copy, new_var) # create copy ncfile
    
    ### write data to copy
    for (j in seq(1, filenc_out$nvars)) { # iterate on variables
        
        var_name = old_var[[j]]$name
        loop_var_old = ncvar_get(filenc_out, var_name)
        
        if (any(old_var[[j]]$dimids == id_ncalib)) { # if N_CALIB is a dimension
            
            loop_var_new = ncvar_get(filenc_copy, var_name)
            arr_indices = which(array(TRUE,dim(loop_var_old)), arr.ind = TRUE) # where do loop_var_old data exist
            loop_var_new[arr_indices] = loop_var_old
            
            ncvar_put(filenc_copy, var_name, loop_var_new)
            
        } else {
            count = old_var[[j]]$varsize
            start = rep(1,length(count))
            ncvar_put(filenc_copy, var_name, loop_var_old, start = start, count = count)
        }
        
        # copy variable attributes
        var_atts_old = attributes(ncatt_get(filenc_out, varid=old_var[[j]]))$names
        var_atts_new = attributes(ncatt_get(filenc_copy, varid=new_var[[j]]))$names
        for (loop_att in var_atts_old) {
            if (!any(var_atts_new==loop_att)) { #if the attribute does not already exist
                ncatt_put(filenc_copy, varid=new_var[[j]], loop_att, ncatt_get(filenc_out, varid=new_var[[j]], attname=loop_att)$value)
            }
        }
        
    }
    
    ### copy global attributes
    all_atts = attributes(ncatt_get(filenc_out, varid=0))$names
    for (loop_att in all_atts) {
        ncatt_put(filenc_copy, varid=0, loop_att, ncatt_get(filenc_out, varid=0, attname=loop_att)$value)
    }
    
    # close files
    nc_close(filenc_out)
    nc_close(filenc_copy)
    
    # overwrite original file
    system2("mv", c(file_out_copy, file_out))
    
    return(0)
}