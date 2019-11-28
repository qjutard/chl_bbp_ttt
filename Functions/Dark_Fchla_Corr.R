##### DARK OFFSET CORRECTION
#
# M.Cornec 20/11/19
#
# Correct the profile of Fchla from the dark offset 
# Default correction: Schmechtig et al., 2014 (substracting the median botom chl values)
# Oxygen Minimum Zone: Wojtasiewic et al., 2018 (local chla minimum as offset value, keep secondary chla peak if existing, set deeper value to 0)
# Subtropical Gyre/Black Sea: Apply Xing et al., 2017 (local chla minimum as offset value, set chl values below this depth to 0)
#
# Inputs: 
# - WMO_pro: reference of the profile, WMO+profile number (11 character, format: "XXXXXXX_XXX")
# - chl: vector of chl values
# - dep_chl: vector of depth values associated to the chl values
# - MLD : MLD value
# - zone : profile regional location (determined with Zone function)
# - DEEP_EST : dataframe with the profiles references needing deep vertical mixing correction and corresponding offsets (determined with Dark_MLD_table_coriolis function)
#
# Output: a vector with the chl values corrected from the dark offset (same length as the input chl vector)
#
Dark_Fchla_Corr <- function (WMO_pro,chl,dep_chl,MLD,zone,DEEP_EST) {
  
  
  chl_dark<-chl #set default vector to be returned
  offset = NA
  min_dep = NA
  
  # Substract a deep offset on the profile, calculated as the median of the values of chl between the max depth (below 900m) and 50m above
  if (max(dep_chl,na.rm=T)>900) { # test if the profile max is deeper than 900m
    dark_deep<-NA
    dark_deep<-median(chl[which(dep_chl<=max(dep_chl,na.rm=T) & 
                                  dep_chl>= (max(dep_chl,na.rm=T)-50))],na.rm=T) # calculate the offset (median of deep values)
    chl_dark<-chl-dark_deep #substract the offset
    offset = dark_deep
  }
  
 if(is.na(zone)==F) { # test if the profile has a regional location
   
   #Test if the profile is located in Oxygen Minimum Zones
   if( zone=="IOMZ" |  zone=="ASEW" | zone=="PSEW") { 
     chl_dark<-chl*NA
     chl_Darkoz = Darkoz(dep_chl,chl,0.5) # Apply Wojtasiewic et al., 2018 correction (OMZ with keeping the secondary peaks)
     chl_dark<-chl_Darkoz$chl_dark
     offset = chl_Darkoz$offset
     min_dep = chl_Darkoz$min_dep
   }
   
   #Test if the profile is located in Subtropical Gyre Zone/ Black Sea
   if( zone=="NPSTG" |  zone=="SPSTG" | zone=="NASTG" | zone=="SASTG" | zone=="SISTG" | zone=="BKS") {
     chl_dark<-chl*NA
     chl_DarkXing = DarkXing(dep_chl,chl) # Apply Xing et al., 2017 correction
     chl_dark<-chl_DarkXing$chl_dark
     offset = chl_DarkXing$offset
     min_dep = chl_DarkXing$min_dep
   }
 }
  
  #Test if the profile is concerned by deep vertical mixing correction
  if(is.null(DEEP_EST)==F) {
    if(WMO_pro %in% DEEP_EST$sub_ref){ # test if the profile ref is in the deep vertical mixing dataframe
      DEEP_EST$dark_est<-as.numeric(paste(DEEP_EST$dark_est)) #convert the dark value from the dataframe into a numeric object
      chl_dark<-chl*NA
      chl_dark<-chl-DEEP_EST$dark_est[which(DEEP_EST$sub_ref==WMO_pro)] # correct the chl profile with the dark mixing value
      offset = DEEP_EST$dark_est[which(DEEP_EST$sub_ref==WMO_pro)]
    }
  }
  
  return(list("chl_dark"=chl_dark, "offset"=offset, "min_dep"=min_dep))
  
  
}