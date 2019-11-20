#### DARK OFFSET CALCULATION IN THE MINIMUM OXYGEN ZONE
#
#
# M.Cornec 20/11/19
#
# Adapted from Wojtasiewicz et al., 2018
#
# Correct chla profile from the dark offset in the Oxygen Minimum Zone, without removing a secondary
# peak of chla if existing
#
# Inputs:
# - dep_chl : a vector of depth values associated with the chl ones
# - chl : a vector of chl values 
# - index: a value of the criterion for testing the presence a secondary peak, based on the coefficient of variation of the chla profile (recommanded 0.5)
#
# Output: a vector of chla values corrected from the dark offset (same length as the input chla vector)

Darkoz<-function(dep_chl,chl,index) {
  
  chl_dark<-chl # set a default chl_dark vector
  
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
    chl_sub1<-NA
    dep_chl_sub1<-NA
    chl_sub1<-chl_filt_res[which(dep_reg_res>Fchla_max)] # retrieve the part of the chla profile below the maximum of chla depth
    dep_chl_sub1<-dep_reg_res[which(dep_reg_res>Fchla_max)] # retrieve the part of the depth profile below the maximum of chla depth
    min_dep<-dep_chl_sub1[min(which.min(chl_sub1))] # indentify the depth where the chla is at its minimum value
    min_chla<-min(chl_sub1,na.rm=T) # retrieve the min chl value (dark offset)
    chl_dark<-chl_dark-min_chla # correct the chla profile from the offset
    
    # Test if second peak of Fchla below the min chla depth (based on a the coefficient of variation of the chla profile below the min chl depth)
    sub_set_chla<-NA
    sub_set_chla<-chl_filt_res[which(dep_reg_res<=500 & dep_reg_res>=min_dep)] # select the portion of the chl profile between the min chl depth and 500m
    sub_set_dep<-NA
    sub_set_dep<-dep_reg_res[which(dep_reg_res<=500 & dep_reg_res>=min_dep)] # select the portion of the depth profile between the min chl depth and 500m
    coef_var<-NA
    coef_var<- sd(sub_set_chla,na.rm=T)/mean(sub_set_chla,na.mr=T)# calculate the coefficient of variation of the sub_profile
    
    if(is.na(coef_var)==F) {
      if(coef_var > index) { # if the coef var is higher than the index criterion, a secondary peak is present
        sub_max<-NA
        sub_max<-sub_set_dep[which.max(sub_set_chla)] # identify the secondary peak depth
        chl_sub2<-NA
        chl_sub2<-chl_filt_res[which(dep_reg_res > sub_max  & dep_reg_res <= 500)] # select the portion of the chl profile between the secondaty chla peak and 500m
        dep_chl_sub2<-NA
        dep_chl_sub2<-dep_reg_res[which(dep_reg_res > sub_max & dep_reg_res <= 500)] # select the portion of the depth profile between the secondaty chla peak and 500m
        min_dep<-dep_chl_sub2[min(which.min(chl_sub2))] # identify a new minimum chla depth below the secondary peak
      }
    }
    chl_dark[which(dep_chl>=min_dep)]<-0 # chl profile below the min chl depth is set to 0
  }
  return(chl_dark)
}