########## CODE TTT BGC-ARGO
#
# M. Cornec 20/11/19
#
# Treatment of the Fchla and bbp profiles
# 
# Structure of the code:
#   
## Set the libraries
##
## Set the functions
##
## Read the merge index file
##
## Set the profiles to be corrected
##
## Loop the treatment:
### A) set the id of the profiles
### B) dark test to prepare further special deep vertical mixing offset on the float time serie
### C) open the file
### D) position informations (lon, lat, time)
### E) physics informations (pressure, depth, temperature, salinity, potential density, Mixed Layer Depth, profile regional location)
### F) PAR profile retrieval (if available)
### G) chl profile retrieval and treatment:
##### 1) range test
##### 2) dark offset
##### 3) & 4) NPQ correction and Factor 2
### H) bbp profile retrieval and treatment:
##### 1) range test
##### 2) special manual drift correction
### I) plot the Fchla and bbp profiles before and after correction
### J) close the netcdf profile

#rm(list = ls())

####################
############### LIBRARIES
###################

require(ncdf4) #deal with netcdf format files
require(oce) #calculate density sigma
require(MASS)

####################
############### FUNCTIONS
###################

#dir_function <- "/PATH_FUNCTIONS/"
dir_function = "~/Documents/cornec_chla_qc/chl_bbp_ttt/Functions/"

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

####################
############### READ THE MERGE FILE INDEX
###################

# Read the merge file index
#index_ifremer<-read.table("~/Documents/data/argo_merge-profile_index.txt", skip=9, sep = ",")

####################
###############  INITIALIZE THE VECTOR OF PROFILES TO TREAT
###################


## !!!!!!!! Initiate the profiles list to be treated
# Must be a vector of profiles ID at the format "WMO_profilenumber.", example : "6901495_025." (12 elements)
# The "." is used to tell that the profile is an Ascent one (should be replaced by "D" if it is a Descent one) 

#profile_list<-c("6901524_150.", # case of deep vertical mixing in the north atlantic subpolar gyre 
#                "6901472_024.", # case of subtropical gyre with increase of Fchla at depth
#                "6901527_040.", # case of NPQ correction with the PAR profile
#                "5904686_040."# case of NPQ correction without the PAR profile
#                )
#profile_list<-c("6901524_150.")
#path_to_netcdf = "~/Documents/data/chla_night_profiles/"
#path_to_netcdf = "/DATA/ftp.ifremer.fr/ifremer/argo/dac/"


########################################################################################################
####################
#################### LOOP THE TREATMENT
#################### 
########################################################################################################


#dark_old<-"XXXXXXX" # initiate the dark marker (use to calculate a dark correction on a float time serie) 

#for (profile_actual in profile_list) {

process_file <- function(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=NULL, index_greylist=NULL, accept_descent=FALSE, plot_chla=FALSE){ 
 
  print(profile_actual)
    
  #################
  ############# A) SET THE ID
  #################
    
  files<-as.character(index_ifremer[,1]) #retrieve the path of each netcfd file
  ident<-strsplit(files,"/") #separate the different roots of the files paths
  ident<-matrix(unlist(ident), ncol=4, byrow=TRUE)
  dac<-ident[,1] #retrieve the DAC of all profiles as a vector
  wod<-ident[,2] #retrieve the WMO of all profiles as a vector
  prof_id<-ident[,4] #retrieve all profiles  name as a vector
  variables<-as.character(index_ifremer[,8]) #retrieve the list of variables available in each file
  variables<-strsplit(variables," ") #separate the different available variables of each profile
  lat<-index_ifremer[,3] #retrieve the latitude of all profiles as a vector
  lon<-index_ifremer[,4] #retrieve the longitude of all profiles as a vector
  prof_date = index_ifremer$V2 #retrieve the date of all profiles as a vector
  
  www<-substr(profile_actual,1,11)
  i <-which(substr(prof_id,3,14)==profile_actual) #identify profile position in the index
  dac_prof<-dac[i] #identify the dac
  
  #################
  ############# A.1) GREYLIST
  #################
  
  chl_greylist_qc = NA
  bbp_greylist_qc = NA
  
  if (!is.null(index_greylist)) {
      
      indices_greylist = which( index_greylist$PLATFORM_CODE==wod[i] & (index_greylist$PARAMETER_NAME=="CHLA" | index_greylist$PARAMETER_NAME=="BBP700") )
      
      prof_date_trunc = as.numeric( substr(as.character(prof_date[i]), 1, 8) )
      
      for (j in indices_greylist) {
          
          is_greylist = FALSE
          
          if (is.na(index_greylist$END_DATE[j])) { # exact date
              if (prof_date_trunc==index_greylist$START_DATE[j]) {
                  is_greylist = TRUE
              } 
          } else { # date interval
              if (index_greylist$START_DATE[j]<=prof_date_trunc & prof_date_trunc<=index_greylist$END_DATE[j]){
                  is_greylist = TRUE
              }
          }
          
          if (is_greylist){
              if (index_greylist$QUALITY_CODE[j] == 4) {
                  print(paste("profile on the greylist with QC 4 at index ", j, " with comment : ", index_greylist$COMMENT[j], sep=""))
                  return(109)
              } else if (index_greylist$QUALITY_CODE == 3){
                  if (index_greylist$PARAMETER_NAME[j]=="CHLA") { chl_greylist_qc = "3" }
                  if (index_greylist$PARAMETER_NAME[j]=="BBP700") { bbp_greylist_qc = "3" }
              } else {
                  print(error_message(110))
                  return(110)
              }
          }
          
          
          
      }
  }
  
  
  
  # Skip if the profile is a descent one (optional)
  if (substr(profile_actual,12,12)=="D" & !accept_descent) {
    print("Descent Profile")
    return(101)
  } 
  
  
  
  #################
  ############# B) DARK TEST FOR DEEP VERTICAL MIXING FURTHER CORRECTION
  #################
  
  # Calculation of the Dark time series to identify if deep vertical mixing offset will be needed on profiles of this float time serie
  # This is done once per float
  #dark_new<-NA
  #dark_new<-substr(profile_actual,1,7) # attribute the actual wmo to the dark marker
  #if (dark_new!=dark_old) { #test if the profile is from a new WMO or not
  #if (is.null(DEEP_EST)){ ### if DEEP_EST is not given as argument, compute it (if several profiles from the same float will be used, it should be calculated once and given as argument)
    #print("dark TS calc")
    #DEEP_EST<-NULL
    #DEEP_EST<-Dark_MLD_table_coriolis(substr(profile_actual,1,7), # calculation of the dark time serie (Dark_MLD_table_coriolis function)
    #                                  path_to_netcdf,index_ifremer)
  #}
  #dark_old<-dark_new #attribute the actual wmo to the dark marker for the next profile
  
  #################
  ############# C) OPEN THE FILE
  #################
  
  path_split = NA
  path_split = unlist( strsplit(files[i],"/") )
  path_to_profile = NA
  path_to_profile = paste(path_split[1], path_split[2], path_split[3], sep="/")
  
  filenc_name_M = NA
  filenc_name_C = NA
  filenc_name_B = NA
  filenc_name_M = path_split[4]
  filenc_name_C = paste("?",substring(filenc_name_M, 3),sep="")
  filenc_name_B = paste("B", filenc_name_C, sep="")
  
  file_M = NA
  file_C = NA
  file_B = NA
  file_M = files[i]
  file_C = paste(path_to_netcdf, path_to_profile,"/", filenc_name_C, sep="") 
  file_C = system2("ls",file_C,stdout=TRUE) # identify R or D file 
  file_B = paste(path_to_netcdf, path_to_profile,"/", filenc_name_B, sep="") 
  file_B = system2("ls",file_B,stdout=TRUE) # identify R or D file 
  
  profile_C<-NULL
  profile_B<-NULL
  profile_C <- nc_open(file_C, readunlim=FALSE, write=FALSE)
  profile_B <- nc_open(file_B, readunlim=FALSE, write=FALSE)
  
  #################
  ############# D) POSITION INFORMATIONS : LON / LAT / DATE
  #################
  
  position_qc<-NA
  position_qc<-substr(ncvar_get(profile_B,"POSITION_QC"),1,1) # read position QC
  
  # skip the profile if the position QC is bad
  if (position_qc == 3 | position_qc==4) {
    print("bad position")
    #nc_close(profile) #close the netcdf
    nc_close(profile_C)
    nc_close(profile_B)
    return(102)
  }
  
  lat<-NA
  lat<- ncvar_get(profile_B,"LATITUDE")[1]
  lon<-NA
  lon<- ncvar_get(profile_B,"LONGITUDE")[1]
  
  # skip the profile if one (or both) coordinate(s) is(are) missing
  if (is.na(lat) | is.na(lon)) {
    print("no geoloc")
    #nc_close(profile) #close the netcdf
    nc_close(profile_C)
    nc_close(profile_B)
    return(103)
  }
  
  jd<-NA
  jd <- ncvar_get(profile_B,"JULD")[1] #read julian day
  origin<-NA # set the origin date
  origin<-as.POSIXct("1950-01-01 00:00:00", order="ymdhms") #convert juld->time
  time<-NA
  time<-origin + jd*3600*24 #calculate the time (format POSIXct yyyy-mm-dd hh:mm:ss)
  jd_qc<-NA
  jd_qc<-substr(ncvar_get(profile_B,"JULD_QC"),1,1) # read julian day qc
  
  # skip the profile if the date is missing
  if (is.na(time)) {
    print("bad date")
    #nc_close(profile) #close the netcdf
    nc_close(profile_C)
    nc_close(profile_B)
    return(104)
  }
  
  # skip the profile if julian date qc is bad
  if (jd_qc == 3 | jd_qc==4) {
    print("bad date")
    #nc_close(profile) #close the netcdf
    nc_close(profile_C)
    nc_close(profile_B)
    return(105)
  }
  
  ###################
  ############# E) PHYSICS INFORMATIONS : DEPTH / TEMP / SAL / POTENTIAL DENSITY / MLD
  ##################
  
  pres<-NA 
  pres <- as.vector(ncvar_get(profile_B,"PRES")) #read the pressure variable as one unique vector
  pres_qc<-NULL #set the qc pressure 
  
  parameters<-ncvar_get(profile_B,"STATION_PARAMETERS") #read the parameters variable (indicating the variables corresponding to each column of the profile file)
  
  # Change the QC pressure to 1 for the bio-optic parameters (remove some Argo processing issue putting wrong QC pressure)
  for (ik in 1:dim(ncvar_get(profile_C,"PRES_QC"))) {
    if (length(grep("CHLA",parameters[,ik]))==1 | #identify the columns where bio-optic parameters are measured
        length(grep("BBP700",parameters[,ik]))==1) {
      optic_depth_qc<-NA
      optic_depth_qc<-paste(rep(1,nchar(ncvar_get(profile_C,"PRES_QC")[1])),collapse="") #create vector with QC 1 
      pres_qc<-paste(pres_qc,optic_depth_qc,sep="") # bind the qc vectors per column into one
      next
    }
    pres_qc<-paste(pres_qc,ncvar_get(profile_C,"PRES_QC")[ik],sep="") # bind the qc vectors per column into one
  }
  
  # put NA to pressure with a bad QC
  for (jj in 1:nchar(pres_qc)) {
    if (substr(pres_qc,jj,jj)==3 | substr(pres_qc,jj,jj)==4) {
      pres[jj]<-NA
    }
  }
  
  pres_na<-na.omit(pres) #remove NAs
  pres_order<-pres_na[order(pres_na)] # order pres
  depth<-swDepth(pres_order,lat) #convert pressure into depth (package oce) according to the latitude
  DEPTH<-unique(depth) # remove duplicata
  
  
  temp_get<-NA
  temp_get <- as.vector(ncvar_get(profile_C,"TEMP")) #read the temperature variable as one unique vector
  
  # read the qc temp
  temp_qc<-NULL 
  for (ik in 1:dim(ncvar_get(profile_C,"TEMP_QC"))) {
    temp_qc<-paste(temp_qc,ncvar_get(profile_C,"TEMP_QC")[ik],sep="")
  }
  
  # attribute NA to temp values with bad qc
  temp_all<-NA
  temp_all<-temp_get
  for (jj in 1:nchar(temp_qc)) {
    if (substr(temp_qc,jj,jj)==3 | substr(temp_qc,jj,jj)==4) {
      temp_all[jj]<-NA
    }
  }
  
  pres_temp<-NA
  pres_temp<-pres[which(!is.na(pres)==T & !is.na(temp_all)==T)] # attribute pressure vector corresponding to the temperature values (with no NAs)
  temp<-NA
  temp<-temp_all[which(!is.na(pres)==T & !is.na(temp_all)==T)] # remove NAs from the temp vector
  temp<-temp[order(pres_temp)] # order the temp vector according to the increasing pressure
  pres_temp<-pres_temp[order(pres_temp)] # order the pressure vector corresponding to the temperature
  
  sal_get<-NA
  sal_get <- as.vector(ncvar_get(profile_C,"PSAL")) #read the salinity variable as one unique vector
  
  # read the qc sal 
  sal_qc<-NULL
  for (ik in 1:dim(ncvar_get(profile_C,"PSAL_QC"))) {
    sal_qc<-paste(sal_qc,ncvar_get(profile_C,"PSAL_QC")[ik],sep="")
  }
  
  # attribute NA to sal values with bad qc
  sal_all<-NA
  sal_all<-sal_get
  for (jj in 1:nchar(sal_qc)) {
    if (substr(sal_qc,jj,jj)==3 | substr(sal_qc,jj,jj)==4) {
      sal_all[jj]<-NA
    }
  }
  
  pres_sal<-NA
  pres_sal<-pres[which(!is.na(pres)==T & !is.na(sal_all)==T)] # attribute pressure vector corresponding to the salinity values (with no NAs)
  sal<-NA
  sal<-sal_all[which(!is.na(pres)==T & !is.na(sal_all)==T)]# remove NAs from the salinity vector
  sal<-sal[order(pres_sal)] # order the sal vector according to the increasing pressure
  pres_sal<-pres_sal[order(pres_sal)] # order the pressure vector corresponding to the salinity
  
  
  # Calculate sigma (potential density)
  for (wii in 1:length(sal_all)) {
    if (is.na(sal_all[wii])==T | is.na(temp_all[wii])==T) { #attribute NA to salinity initial vector where there is corresponding NA in the initial temp vector (and inverse)
      sal_all[wii]<-NA
      temp_all[wii]<-NA
    }
  }
  sigma_all<-NA
  sigma_all <- swSigmaTheta(sal_all,temp_all,pres) # calculation of sigma (package oce)
  pres_sigma<-NA
  pres_sigma<-pres[which(!is.na(pres)==T & !is.na(sigma_all)==T)] # attribute pressure vector corresponding to the sigmainity values (with no NAs)
  sigma<-NA
  sigma<-sigma_all[which(!is.na(pres)==T & !is.na(sigma_all)==T)]# remove NAs from the sigmainity vector
  sigma<-sigma[order(pres_sigma)] # order the sigma vector according to the increasing pressure
  pres_sigma<-pres_sigma[order(pres_sigma)] # order the pressure vector corresponding to the sigmainity
  
  # calculate depth from pressure for the different depth vector (package oce) according to the latitude
  dep_temp<-NA
  dep_temp = swDepth(pres_temp,lat)
  dep_sal<-NA
  dep_sal = swDepth(pres_sal,lat)
  dep_sigma<-NA
  dep_sigma = swDepth(pres_sigma,lat)
  
  # Zone calculation from Zone function
  zone<-NA
  zone<-Zone(lat,lon,temp,dep_temp)
  
  # Calculate of MLD (crit diff dty 0.03) from MLD_calc function
  MLD<-NA
  MLD<-MLD_calc(sigma,dep_sigma)
  
  
  ##############################
  ################### F) PAR RETRIEVAL (IF AVAILABLE)
  ##############################
  # Needed for the NPQ correction 
  
  if ("DOWNWELLING_PAR" %in% names(profile_B$var)==T) { # test if the variable is present in the netcdf file
    light_get<-NA
    light_get <- as.vector(ncvar_get(profile_B,"DOWNWELLING_PAR"))  #read the PAR variable as one unique vector
    
    # read the qc varaible as a vector
    light_qc<-NULL
    for (ik in 1:dim(ncvar_get(profile_B,"DOWNWELLING_PAR_QC"))) {
      light_qc<-paste(light_qc,ncvar_get(profile_B,"DOWNWELLING_PAR_QC")[ik],sep="")
    }
    
    # attribute NA to par values with bad qc
    light_all<-NA
    light_all<-light_get
    for (jj in 1:nchar(light_qc)) {
      if (substr(light_qc,jj,jj)==3 | substr(light_qc,jj,jj)==4) {
        light_all[jj]<-NA
      }
    }
    
    pres_light<-NA
    pres_light<-pres[which(!is.na(pres)==T & !is.na(light_all)==T)] # attribute pressure vector corresponding to the light values (with no NAs)
    light<-NA
    light<-light_all[which(!is.na(pres)==T & !is.na(light_all)==T)]# remove NAs from the light vector
    light<-light[order(pres_light)] # order the light vector according to the increasing pressure
    pres_light<-pres_light[order(pres_light)] # order the pressure vector corresponding to the light
    
    # calculate depth from pressure (package oce) according to the latitude
    dep_light<-NA
    dep_light = swDepth(pres_light,lat)
  }
  
  ###############################
  ######################### G) CHL RETRIEVAL AND TREATMENT 
  ###############################
  
  if ("CHLA" %in% names(profile_B$var)==F) { # test if the chla variable is present in the netcdf file (skip if not)
    print("no chl available")
    #nc_close(profile) #close the netcdf
    nc_close(profile_C)
    nc_close(profile_B)
    return(106)
  }
  
  chl_get<-NA
  chl_get <- as.vector(ncvar_get(profile_B,"CHLA"))  #read the chla variable as one unique vector
  chl_all<-NA
  chl_all<-chl_get
  
  ############ 1) RANGE TEST : attribute NA to values out of range ###########################
  chl_all[which((chl_all > 50) | (chl_all < - 0.1))]<-NA
  
  chl_not_isna = NA
  chl_not_isna = which(!is.na(pres)==T & !is.na(chl_all)==T) # keep the information on NA values to later add them
  
  pres_chl<-NA
  pres_chl<-pres[chl_not_isna] # attribute pressure vector corresponding to the chl values (with no NAs)
  chl<-NA
  chl<-chl_all[chl_not_isna]# remove NAs from the chl vector
  pres_chl_unsorted = NA
  pres_chl_unsorted = pres_chl # save the vector before it was sorted to allow for unsorting
  chl<-chl[order(pres_chl)] # order the chl vector according to the increasing pressure
  pres_chl<-pres_chl[order(pres_chl)] # order the pressure vector corresponding to the chl
  
  # calculate depth from pressure (package oce) according to the latitude
  dep_chl<-NA
  dep_chl = swDepth(pres_chl,lat)
  
  # Test if the chl values are associated to only one depth (error of the measurement: "stuck pressure); if so, skip
  if(length(unique(dep_chl))==1) {
    print("stuck pressure")
    #nc_close(profile) #close the netcdf file
    nc_close(profile_C)
    nc_close(profile_B)
    return(107)
  }
  
  
  ################ 2) DARK OFFSET #######################
  # correct the vertical profile from a deep offset (Dark_Fchla_Corr function)
  chl_dark<-NA
  list_dark = Dark_Fchla_Corr(substr(profile_actual,1,11),chl,dep_chl,MLD,zone,DEEP_EST)
  chl_dark<-list_dark$chl_dark
  
  chl_dark_offset = NA 
  chl_dark_min_pres = NA
  chl_dark_offset = list_dark$offset # QJ : for scientific_calib_*
  chl_dark_min_pres = list_dark$min_dep 
  if (!is.na(chl_dark_min_pres)) {
      chl_dark_min_pres = swPressure(chl_dark_min_pres, lat) # invert swDepth
  }
  
  ############ 3) & 4) NPQ and FACTOR 2 ########################
  # Correct the chla profile from the:
  # - Non Photochemical Quenching: Xing et al., 2018
  # - Factor 2: Roesler et al., 2017
  
  #Pre-calculate RESO which will cause an error in the chl_npq computation if RESO=0
  RESO = NA
  RESO = mean(abs(diff(dep_chl[which(dep_chl<=250)])),na.rm=T) 
  if (RESO==0 | is.na(RESO)){
      print("RESO = 0")
      nc_close(profile_C)
      nc_close(profile_B)
      return(108)
  }
  
  if ("DOWNWELLING_PAR" %in% names(profile_B$var)==T) { # test if there is a PAR in situ measured
    chl_npq_list<-NPQ_cor_X12_XB18(chl_dark/2,dep_chl,dep_light,light,MLD)
  } else {
    chl_npq_list<-NPQ_cor_P18(chl_dark/2,dep_chl,MLD)
  }
  print(chl_npq_list)
  
  print(chl_dark/2)
  
  chl_npq=chl_npq_list$chl_npq
  
  flag_NPQ_changed = NA
  flag_NPQ_changed = rep(FALSE, length(chl_npq))
  chl_max_npq_depth = chl_npq_list$chl_max_npq_depth
  if (!is.na(chl_max_npq_depth)) {
    flag_NPQ_changed[which(dep_chl <= chl_max_npq_depth )] = TRUE # save the values that were changed by the NPQ correction to later write QC
  }
  ############################
  ############ H) BBP700 RETRIEVAL AND TREATMENT
  ############################
  
  bbp_get<-NA
  bbp_get <- as.vector(ncvar_get(profile_B,"BBP700")) #read the chla variable as one unique vector
  
  bbp_all<-NA
  bbp_all<-bbp_get
  
  ################### 1) RANGE TEST ###############################
  bbp_all[which((bbp_all > 0.1) | (bbp_all <  (-0.000005)))]<-NA
  
  bbp_not_isna = NA
  bbp_not_isna = which(!is.na(pres)==T & !is.na(bbp_all)==T) # keep the information on NA values to later add them
  
  pres_bbp<-NA
  pres_bbp<-pres[bbp_not_isna] # attribute pressure vector corresponding to the bbp values (with no NAs)
  bbp<-NA
  bbp<-bbp_all[bbp_not_isna]# remove NAs from the bbp vector
  pres_bbp_unsorted = NA
  pres_bbp_unsorted = pres_bbp # save the vector before it was sorted to allow for unsorting
  bbp<-bbp[order(pres_bbp)] # order the bbp vector according to the increasing pressure
  pres_bbp<-pres_bbp[order(pres_bbp)] # order the pressure vector corresponding to the bbp
  
  # calculate depth from pressure (package oce) according to the latitude
  dep_bbp<-NA
  dep_bbp = swDepth(pres_bbp,lat)
  
  
  #### 2) SPECIAL DRIFT CORRECTION (MANUALLY CALCULATED) #############################
  
  diff_bottom = NA
  
  if (substr(profile_actual,1,7)==5904218 & substr(profile_actual,1,11) < 660 & substr(profile_actual,1,11) > 456 ) {
    med_bottom<-NULL
    med_bottom<-median(bbp[which(dep_bbp < 300 & dep_bbp > 250)], na.rm=T)
    diff_bottom<-0-med_bottom
    bbp<-bbp+diff_bottom
  }
  
  if (substr(profile_actual,1,7)==2902092 & substr(profile_actual,1,11) < 125 & substr(profile_actual,1,11) > 97 ) {
    med_bottom<-NULL
    med_bottom<-median(bbp[which(dep_bbp < 1600 & dep_bbp > 1400)], na.rm=T)
    diff_bottom<-0.000035-med_bottom
    bbp<-bbp+diff_bottom
  }
  
  if (substr(profile_actual,1,7)==6901174 & format(time,"%Y-%m-%d") >= as.POSIXct("06/03/18",format="%d/%m/%Y")) {
    med_bottom<-NULL
    med_bottom<-median(bbp[which(dep_bbp < 950 & dep_bbp > 850)], na.rm=T)
    diff_bottom<-0.0011-med_bottom
    bbp<-bbp+diff_bottom
  }
  
  if (substr(profile_actual,1,7)==2902118 ) {
    med_bottom<-NULL
    med_bottom<-median(bbp[which(dep_bbp < 950 & dep_bbp > 850)], na.rm=T)
    diff_bottom<-0.00035-med_bottom
    bbp<-bbp+diff_bottom
  }
  
  if (substr(profile_actual,1,7)==6901485 & format(time,"%Y-%m-%d") >= as.POSIXct("08/05/15",format="%d/%m/%Y") &
      format(time,"%Y-%m-%d")  <= as.POSIXct("28/07/15",format="%d/%m/%Y")) {
    med_bottom<-NULL
    med_bottom<-median(bbp[which(dep_bbp < 1600 & dep_bbp > 1400)], na.rm=T)
    diff_bottom<-0.0012-med_bottom
    bbp<-bbp+diff_bottom
  }
  
  if (substr(profile_actual,1,7)==6901485 & format(time,"%Y-%m-%d") >= as.POSIXct("01/09/15",format="%d/%m/%Y") &
      format(time,"%Y-%m-%d")  <= as.POSIXct("04/03/16",format="%d/%m/%Y")) {
    med_bottom<-NULL
    med_bottom<-median(bbp[which(dep_bbp < 1600 & dep_bbp > 1400)], na.rm=T)
    diff_bottom<-0.0012-med_bottom
    bbp<-bbp+diff_bottom
  }
  
  ############################
  ############ I) PLOT THE PROFILES BEFORE/AFTER CORRECTION
  ############################
  if (plot_chla) {
      par(mar=c(5,6.5,5,1.5))
      plot(chl_npq,-dep_chl,
           xlab = NA, ylab=NA,col="red",pch=16,xaxt="n",yaxt="n",cex=0.8,type="p",
           ylim=c(-400,0),
           # xlim=(c(min(chl,na.rm = T),
           #         max(chl,na.rm = T))),las=1)
           xlim=(range(chl)),las=1)
  
      axis(3,col.axis = "black",cex.axis=2)
      mtext(expression(paste("Fchl",italic(a),sep="")~(mg~chl~m^{-3})), side=3, line=2.7, cex = 2, col="black")
      axis(2,cex.axis=2, las=1)
      mtext("Depth (m)", side=2, line=4.2, cex = 2)
      par(new=T)
      plot(chl,-dep_chl,cex=0.8,type="p",pch=16,
           xlab = NA, ylab=NA,xaxt="n",yaxt="n",col="black",
           ylim=c(-400,0),
           xlim=(range(chl)),las=1)
      mtext(www, side=1, line=2, cex = 1.7)
  }
  
  ############################
  ############ J) REFORMAT TO NETCDF ORIGINAL SIZES
  ############################
  
  chl_fin = chl_npq # chl after all corrections were applied
  bbp_fin = bbp     # bbp after all corrections were applied
  
  ### unsort vectors
  chl_unsorted = chl_fin
  bbp_unsorted = bbp_fin
  chl_unsorted[order(pres_chl_unsorted)] = chl_fin
  bbp_unsorted[order(pres_bbp_unsorted)] = bbp_fin
  
  ### add NA values
  chl_with_na = rep(NA,length(chl_get))
  bbp_with_na = rep(NA,length(bbp_get))
  chl_with_na[chl_not_isna] = chl_unsorted
  bbp_with_na[bbp_not_isna] = bbp_unsorted
  
  ### separate the profiles
  chl_get_array = ncvar_get(profile_B, "CHLA")
  bbp_get_array = ncvar_get(profile_B, "BBP700")
  chl_array = array(chl_with_na, dim(chl_get_array))
  bbp_array = array(bbp_with_na, dim(bbp_get_array))
  
  ############################
  ############ K) create flags
  ############################
  
  chl_qc = array(" ", dim(chl_get_array))
  bbp_qc = array(" ", dim(bbp_get_array))
  
  ### chl flags
  chl_qc[which(!is.na(chl_array))] = "2" # start with a 2 if a value exists
  chl_qc[which( is.na(chl_array) & !is.na(chl_get_array) )] = "4" # write a "4" if the method removed a value
  
  # place the NPQ flags on the original array
  flag_unsorted = flag_NPQ_changed
  flag_unsorted[order(pres_chl_unsorted)] = flag_NPQ_changed
  flag_with_na = rep(NA,length(chl_get))
  flag_with_na[chl_not_isna] = flag_unsorted
  flag_array = array(flag_with_na, dim(chl_get_array))
  
  chl_qc[which( flag_array )] = "5" # write a "5" where the NPQ correction changed a value
  
  if (!is.na(chl_greylist_qc)){
      chl_qc[which(chl_qc=="2")] = chl_greylist_qc
  }
  
  ### bbp flags
  bbp_qc[which(!is.na(bbp_array))] = "1" # start with a 1 if a value exists
  bbp_qc[which( is.na(bbp_array) & !is.na(bbp_get_array) )] = "4" # write a "4" if the method removed a value
  
  if (!is.na(bbp_greylist_qc)){
      bbp_qc[which(bbp_qc=="1")] = bbp_greylist_qc
  }
  
  ### convert to string-words
  chl_adjusted_qc = rep(NA, dim(chl_qc)[2])
  bbp_adjusted_qc = rep(NA, dim(bbp_qc)[2]) #dim(bbp_qc)==dim(chl_qc)
  for (j in seq(1,dim(chl_qc)[2])){
      chl_adjusted_qc[j] = paste(chl_qc[,j],collapse="")
      bbp_adjusted_qc[j] = paste(bbp_qc[,j],collapse="")
  }
  
  ############################
  ############ K) create error arrays
  ############################
  
  ### open the metadata
  path_to_meta = paste(path_to_netcdf, path_split[1], "/", path_split[2], "/", substring(profile_actual,1,7), "_meta.nc",sep="") # build meta file adress
  
  meta_nc = nc_open(path_to_meta) # open meta file
  
  calib = ncvar_get(meta_nc,"PREDEPLOYMENT_CALIB_COEFFICIENT") # get the cailbration coefficients
  
  ### get the scales
  
  # chla scale
  id_chla = grep("CHLA", calib) # find chla index
  chla_calib = calib[id_chla] # get chla calibration
  chla_calib = unlist(strsplit(chla_calib,",")) # separate coefficients
  chla_calib_scale = chla_calib[grep("SCALE_CHLA",chla_calib)] # get the scale information
  chla_calib_scale = unlist(strsplit(chla_calib_scale,"=")) # separate the name from the number
  
  chla_scale = as.numeric(chla_calib_scale[2]) # get the scale coefficient as a number
  
  # bbp scale
  id_bbp = grep("BACKSCATTERING700", calib) # find bbp index
  bbp_calib = calib[id_bbp] # get chla calibration
  bbp_calib = unlist(strsplit(bbp_calib,",")) # separate coefficients
  bbp_calib_scale = bbp_calib[grep("SCALE_BACKSCATTERING700",bbp_calib)] # get the scale information
  bbp_calib_scale = unlist(strsplit(bbp_calib_scale,"=")) # separate the name from the number
  
  bbp_scale = as.numeric(bbp_calib_scale[2]) # get the scale coefficient as a number
  
  nc_close(meta_nc) # close meta
  
  ### create error arrays
  chl_error = pmax( 3 * array(chla_scale, dim(chl_array)), 1.5 * abs(chl_array))
  bbp_error = pmax( 3 * array(bbp_scale, dim(bbp_array)), 0.2 * abs(bbp_array))
  
  ############################
  ############ L) GET NPQ RELATED INFO
  ############################
  
  is_npq = NA
  npq_depth = NA
  npq_val = NA
  is_npq = any(flag_NPQ_changed)
  if (is_npq) {
      npq_depth = max(pres_chl[which(flag_NPQ_changed)])
      npq_val = chl_npq[which(flag_NPQ_changed)][1]
  }
  
  ############################
  ############ M) CLOSE THE NETCDF PROFILE
  ############################
  #nc_close(profile)
  nc_close(profile_C)
  nc_close(profile_B)
  
  return(list("CHLA_ADJUSTED"=chl_array, "BBP700_ADJUSTED"=bbp_array, 
              "CHLA_ADJUSTED_QC"=chl_adjusted_qc, "BBP700_ADJUSTED_QC"=bbp_adjusted_qc,
              "CHLA_ADJUSTED_ERROR"=chl_error, "BBP700_ADJUSTED_ERROR"=bbp_error,
              "chl_dark_offset"=chl_dark_offset, "chl_dark_min_pres"=chl_dark_min_pres,"bbp_offset"=-diff_bottom,
              "is_npq"=is_npq, "npq_depth"=npq_depth, "npq_val"=npq_val))
  
}

#profile_actual = profile_list[1]
#DEEP_EST = Dark_MLD_table_coriolis(substr(profile_actual,1,7), path_to_netcdf, index_ifremer)

#L = process_file(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=DEEP_EST)

