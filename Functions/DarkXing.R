# DARK OFFSET XING CORRECTION
#
#
# M.Cornec 20/11/19
#
# Adapted from Xing et al., 2017
#
# Correct chla profile from the dark offset in the subtropical gyres and Black Sea, where the chla increases at depth
# Identifies the depth of the minimum chla value (used as offset value to correct the whole profile) and set all the chl values below this depth to 0
#
# Inputs:
# - dep_chl : a vector of depth values associated with the chl ones
# - chl : a vector of chl values 
# 
# Output: a vector of chla values corrected from the dark offset (same length as the input chla vector)


DarkXing<-function(dep_chl,chl) {
  
  chl_dark<-chl # set the default output vector
  offset = NA
  min_dep = NA
  
  # calculate the mean vertical resolution of measurements in the 500 first meters
  step_res<-NA
  step_res<-mean(abs(diff(dep_chl[which(dep_chl<=500)])),na.rm = T)
  
  if (step_res!=0) { # test if the resolution is non null
    dep_reg_res<-NA
    dep_reg_res <-seq(0,500,step_res) # set a interpolated depth vector from surface to 500m with steps set as the mean vertical resolution
    chl_interp_res<-NA
    chl_interp_res = approx(dep_chl,chl,dep_reg_res)$y # interpolate the chl values on the interpolated depth vector
    
    # Smooth the interpolated chla profile with a median running filter with a window depending on the mean vertical resolution of the profile
    chl_filt_res<-NA
    if (step_res < 3) {
      chl_filt_res = RunningFilter(5,chl_interp_res, Method="Median")
    }
    if (step_res > 3) {
      chl_filt_res = RunningFilter(1,chl_interp_res, Method="Median")
    }
    
    # Identify the main deep chl maximum (DCM) depth
    Fchla_max<-NA
    Fchla_max<-dep_reg_res[which.max(chl_filt_res)]
    
    # Retrieve the min chla value below the DCM and correct the chla profile from this offset  
    min_chla<-NA
    min_dep<-NA
    chl_sub<-NA
    dep_chl_sub<-NA
    chl_sub<-chl_filt_res[which(dep_reg_res>Fchla_max)] # retrieve the part of the chla profile below the maximum of chla depth
    dep_chl_sub<-dep_reg_res[which(dep_reg_res>Fchla_max)] # retrieve the part of the depth profile below the maximum of chla depth
    min_dep<-dep_chl_sub[min(which.min(chl_sub))] # indentify the depth where the chla is at its minimum value
    min_chla<-min(chl_sub,na.rm=T) # retrieve the min chl value (dark offset)
    chl_dark<-chl_dark-min_chla # correct the chla profile from the offset
    offset = min_chla
    
    chl_dark[which(dep_chl>=min_dep)]<-0 # chl profile below the min chl depth is set to 0
  }
  return(list("chl_dark"=chl_dark, "offset"=offset, "min_dep"=min_dep))
}