######## NPQ CORRECTION XING 2018 (WITHOUT MEASURMENT)
#
# Derived from Xing et al., 2018
# M.Cornec 20/11/19
#
#
# Correct the NPQ using the PAR profile and the MLD:
# Extrapolation of max Fchl value between min(MLD,ziPAR15) and surface
#
# Inputs:
# - chl: a vector with chl values
# - dep_chl: a vector with depth values corresponding to the chl vector
# - MLD: the MLD depth
#
# Output: a vector of chl corrected from the NPQ (same length as the input chl vector)

####### NPQ_COR_P18 
NPQ_cor_P18 <- function (chl,dep_chl,MLD) {
  
  chl_npq<-chl # set the default output vector
  
  MLD<-0.9*MLD # set the MLD value to 90% of its depth
  
  chl_max_npq_depth<-NA
  which_is_npq_changed = NA
  
  if (length(chl[!is.na(chl)])>5) { # at least 5 no-NA values needed in Fchl profile
    
    #create interpolated chl profile at a resolution of 1m
    chl_int<-NA
    chl_int<-approx(dep_chl,chl,seq(0,500,1))$y
    
    if (is.na(MLD)==F) { 
      
      # Estimate a KD from Kim et al., 2005 (relation between chl and attenuation), from the chl profile
      KD_PAR<-NA     
      KD_PAR<-(0.0232 + (0.074*(chl_int^0.674)))
      
      ### Put 0 values for NA
      KD_PAR[which(is.na(KD_PAR))]<-0
      
      ### Compute the integral of PAR
      int_PAR <- cumsum(KD_PAR)
      
      ### Finally, get the euphotic depth
      Zeu<-NA
      Zeu<-seq(0,500,1)[which(int_PAR==min(int_PAR[which(int_PAR > 4.6)], na.rm=T))]
      
      #vertical resolution of the Fchl profile (above 250m)
      RESO<-NA
      RESO<-mean(abs(diff(dep_chl[which(dep_chl<=250)])),na.rm=T) 
      
      #create interpolated depth and chl profile (fonction of mean vertical resolution)
      dep_reg_res <-seq(0,500,RESO)
      chl_interp_res<-NA
      chl_interp_res = approx(dep_chl,chl,dep_reg_res)$y
      
      # Filter to remove spike (according to vertical resolution)
      chl_filt_res<-NA
      if (RESO < 3) {
        chl_filt_res = RunningFilter(5,chl_interp_res, Method="Median")
      }
      if (RESO > 3) {
        chl_filt_res = RunningFilter(1,chl_interp_res, Method="Median")
      }
      
      ### shallower depth btw MLD and Zeu
      min_depth_npq<-min(MLD,Zeu,na.rm=T)
      
      # identify the value of max Fchla bewteen the surface and the min depth 
      chl_max_npq<-NA
      chl_max_npq<-max(chl_filt_res[which(dep_reg_res<=min_depth_npq)], na.rm=T)
      
      #identify the depth related to the max Fchla
      chl_max_npq_depth<-dep_reg_res[which(chl_filt_res==chl_max_npq & dep_reg_res<=min_depth_npq)]
      chl_max_npq_depth<-max(chl_max_npq_depth,na.rm=T)
      
      # correct the profile with the max Fchla
      if (chl_max_npq!=-Inf) {
        chl_npq[which(dep_chl <= chl_max_npq_depth )]<-chl_max_npq
        which_is_npq_changed = which(dep_chl <= chl_max_npq_depth )
      }
    }
  }
  return(list("chl_npq"=chl_npq, "which_is_npq_changed"=which_is_npq_changed))
}