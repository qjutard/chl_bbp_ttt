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

process_file <- function(profile_actual, index_ifremer, path_to_netcdf, DEEP_EST=NULL, index_greylist=NULL, 
                         accept_descent=FALSE, accept_QC3=FALSE, position_override=NULL, offset_override=NULL, 
                         date_override=NULL, plot_chla=FALSE){ 
 
  #print(profile_actual)
    
  #################
  ############# A) SET THE ID
  #################
    
  files = as.character(index_ifremer$file) #retrieve the path of each netcfd file
  ident = strsplit(files,"/") #separate the different roots of the files paths
  ident = matrix(unlist(ident), ncol=4, byrow=TRUE)
  dac = ident[,1] #retrieve the DAC of all profiles as a vector
  wod = ident[,2] #retrieve the WMO of all profiles as a vector
  prof_id = ident[,4] #retrieve all profiles  name as a vector
  variables = as.character(index_ifremer$parameters) #retrieve the list of variables available in each file
  variables = strsplit(variables," ") #separate the different available variables of each profile
  lat = index_ifremer$latitude #retrieve the latitude of all profiles as a vector
  lon = index_ifremer$longitude #retrieve the longitude of all profiles as a vector
  prof_date = index_ifremer$date #retrieve the date of all profiles as a vector
  
  www = substr(profile_actual,1,11)
  iii = which(substr(prof_id,3,14)==profile_actual) #identify profile position in the index
  dac_prof = dac[iii] #identify the dac
 
  # Skip if the profile is a descent one (optional)
  if (substr(profile_actual,12,12)=="D" & !accept_descent) {
    print(error_message(101))
    return(101)
  } 
  
  #################
  ############# C) OPEN THE FILE
  #################
  
  file_B = NA
  file_B = paste(path_to_netcdf, files[iii], sep="") 
  
  path_split = NA
  path_split = unlist( strsplit(files[iii], "/") )
  path_to_profile = NA
  path_to_profile = paste(path_split[1], path_split[2], path_split[3], sep="/")
  
  filenc_name_C = NA
  file_C = NA
  filenc_name_C = paste("?",substring(path_split[4], 3),sep="")
  file_C = paste(path_to_netcdf, path_to_profile, "/", filenc_name_C, sep="") 
  file_C = system2("ls", file_C, stdout=TRUE) # identify R or D file 
  
  if (length(file_C)==2) { # if both R and D files exist
      file_C = file_C[1] # use the D file which is first in alphabetical order
  }
  
  profile_C = NULL
  profile_B = NULL
  profile_C = nc_open(file_C, readunlim=FALSE, write=FALSE)
  profile_B = nc_open(file_B, readunlim=FALSE, write=FALSE)
  
  #################
  ############# D) POSITION INFORMATIONS : LON / LAT / DATE
  #################
  
  lat = NA
  lat = ncvar_get(profile_B,"LATITUDE")[1]
  lon = NA
  lon = ncvar_get(profile_B,"LONGITUDE")[1]
  
  position_qc = NA
  position_qc = substr(ncvar_get(profile_B,"POSITION_QC"),1,1) # read position QC
  
  # skip the profile if the position QC is bad
  if ( (position_qc == 3 | position_qc==4) & is.null(position_override) ) {
    print(error_message(102))
    nc_close(profile_C)
    nc_close(profile_B)
    return(102)
  }
   
  # skip the profile if one (or both) coordinate(s) is(are) missing
  if ( (is.na(lat) | is.na(lon)) & is.null(position_override) ) {
    print(error_message(103))
    nc_close(profile_C)
    nc_close(profile_B)
    return(103)
  }
  
  if (!is.null(position_override)) {
    lat = position_override[1]
    lon = position_override[2]
  }
  
  jd = NA
  jd = ncvar_get(profile_B,"JULD")[1] #read julian day
  origin = NA # set the origin date
  origin = as.POSIXct("1950-01-01 00:00:00", order="ymdhms", tz="UTC") #convert juld->time
  time = NA
  time = origin + jd*3600*24 #calculate the time (format POSIXct yyyy-mm-dd hh:mm:ss)
  jd_qc = NA
  jd_qc = substr(ncvar_get(profile_B,"JULD_QC"),1,1) # read julian day qc
  
  # skip the profile if the date is missing
  if (is.na(time) & is.null(date_override)) {
      print(error_message(104))
      nc_close(profile_C)
      nc_close(profile_B)
      return(104)
  }
  
  # skip the profile if julian date qc is bad
  if ((jd_qc == 3 | jd_qc==4) & is.null(date_override)) {
      print(error_message(105))
      nc_close(profile_C)
      nc_close(profile_B)
      return(105)
  }
  
  if (!is.null(date_override)) {
      time = as.POSIXct(date_override, format="%Y-%m-%d;%H:%M:%S", tz="UTC") #format POSIXct yyyy-mm-dd hh:mm:ss
  }
  
  #################
  ############# D.1) GREYLIST
  #################
  
  chl_greylist_qc = NA
  bbp_greylist_qc = NA
  
  if (!is.null(index_greylist)) {
      
      indices_greylist = which( index_greylist$PLATFORM_CODE==wod[iii] & (index_greylist$PARAMETER_NAME=="CHLA" | index_greylist$PARAMETER_NAME=="BBP700") )
      
      prof_date_trunc = as.numeric( substr(as.character(prof_date[iii]), 1, 8) )
      
      for (j in indices_greylist) {
          
  		  ## is the profile on the greylist ?
          if (is.na(index_greylist$END_DATE[j])) { # all past that date
              is_greylist = (prof_date_trunc>=index_greylist$START_DATE[j])
          } else { # date interval
              is_greylist = (index_greylist$START_DATE[j]<=prof_date_trunc & prof_date_trunc<=index_greylist$END_DATE[j])
          }
          
		  ## if so, what is the QC and what to do ?
          if (is_greylist){
              
              if (index_greylist$QUALITY_CODE[j] == 4) {
                  
                  print(paste("profile on the greylist with QC 4 at index ", j, " with comment : ", index_greylist$COMMENT[j], sep=""))
                  return(109)
                  
              } else if (index_greylist$QUALITY_CODE[j] == 3) {
                  
                  if (!accept_QC3) {
                      print(paste("profile on the greylist with QC 3 at index ", j, " with comment : ", index_greylist$COMMENT[j], sep=""))
                      return(111)
                  }
                  
                  if (index_greylist$PARAMETER_NAME[j]=="CHLA") { chl_greylist_qc = "3" }
                  if (index_greylist$PARAMETER_NAME[j]=="BBP700") { bbp_greylist_qc = "3" }
                  
              } else {
                  print(error_message(110))
                  return(110)
              }
          }     
      }
  }
  
  ###################
  ############# E) PHYSICS INFORMATIONS : DEPTH / TEMP / SAL / POTENTIAL DENSITY / MLD
  ##################
  
  pres = NA 
  pres = as.vector(ncvar_get(profile_B,"PRES")) #read the pressure variable as one unique vector
  pres_qc = NULL #set the qc pressure 
  
  parameters = ncvar_get(profile_B,"STATION_PARAMETERS") #read the parameters variable (indicating the variables corresponding to each column of the profile file)
  
  # Change the QC pressure to 1 for the bio-optic parameters (remove some Argo processing issue putting wrong QC pressure)
  for (ik in 1:dim(ncvar_get(profile_C,"PRES_QC"))) {
    if (length(grep("CHLA",parameters[,ik]))==1 | #identify the columns where bio-optic parameters are measured
        length(grep("BBP700",parameters[,ik]))==1) {
      optic_depth_qc = NA
      optic_depth_qc = paste(rep(1,nchar(ncvar_get(profile_C,"PRES_QC")[1])),collapse="") #create vector with QC 1 
      pres_qc = paste(pres_qc,optic_depth_qc,sep="") # bind the qc vectors per column into one
      next
    }
    pres_qc = paste(pres_qc,ncvar_get(profile_C,"PRES_QC")[ik],sep="") # bind the qc vectors per column into one
  }
  
  # put NA to pressure with a bad QC
  for (jj in 1:nchar(pres_qc)) {
    if (substr(pres_qc,jj,jj)==3 | substr(pres_qc,jj,jj)==4) {
      pres[jj] = NA
    }
  }
  
  pres_na = na.omit(pres) #remove NAs
  pres_order = pres_na[order(pres_na)] # order pres
  depth = swDepth(pres_order,lat) #convert pressure into depth (package oce) according to the latitude
  DEPTH = unique(depth) # remove duplicata
  
  
  temp_get = NA
  temp_get  =  as.vector(ncvar_get(profile_C,"TEMP")) #read the temperature variable as one unique vector
  
  # read the qc temp
  temp_qc = NULL 
  for (ik in 1:dim(ncvar_get(profile_C,"TEMP_QC"))) {
    temp_qc = paste(temp_qc,ncvar_get(profile_C,"TEMP_QC")[ik],sep="")
  }
  
  # attribute NA to temp values with bad qc
  temp_all = NA
  temp_all = temp_get
  for (jj in 1:nchar(temp_qc)) {
    if (substr(temp_qc,jj,jj)==3 | substr(temp_qc,jj,jj)==4) {
      temp_all[jj] = NA
    }
  }
  
  pres_temp = NA
  pres_temp = pres[which(!is.na(pres)==T & !is.na(temp_all)==T)] # attribute pressure vector corresponding to the temperature values (with no NAs)
  temp = NA
  temp = temp_all[which(!is.na(pres)==T & !is.na(temp_all)==T)] # remove NAs from the temp vector
  temp = temp[order(pres_temp)] # order the temp vector according to the increasing pressure
  pres_temp = pres_temp[order(pres_temp)] # order the pressure vector corresponding to the temperature
  
  sal_get = NA
  sal_get  =  as.vector(ncvar_get(profile_C,"PSAL")) #read the salinity variable as one unique vector
  
  # read the qc sal 
  sal_qc = NULL
  for (ik in 1:dim(ncvar_get(profile_C,"PSAL_QC"))) {
    sal_qc = paste(sal_qc,ncvar_get(profile_C,"PSAL_QC")[ik],sep="")
  }
  
  # attribute NA to sal values with bad qc
  sal_all = NA
  sal_all = sal_get
  for (jj in 1:nchar(sal_qc)) {
    if (substr(sal_qc,jj,jj)==3 | substr(sal_qc,jj,jj)==4) {
      sal_all[jj] = NA
    }
  }
  
  pres_sal = NA
  pres_sal = pres[which(!is.na(pres)==T & !is.na(sal_all)==T)] # attribute pressure vector corresponding to the salinity values (with no NAs)
  sal = NA
  sal = sal_all[which(!is.na(pres)==T & !is.na(sal_all)==T)]# remove NAs from the salinity vector
  sal = sal[order(pres_sal)] # order the sal vector according to the increasing pressure
  pres_sal = pres_sal[order(pres_sal)] # order the pressure vector corresponding to the salinity
  
  
  # Calculate sigma (potential density)
  for (wii in 1:length(sal_all)) {
    if (is.na(sal_all[wii])==T | is.na(temp_all[wii])==T) { #attribute NA to salinity initial vector where there is corresponding NA in the initial temp vector (and inverse)
      sal_all[wii] = NA
      temp_all[wii] = NA
    }
  }
  sigma_all = NA
  sigma_all = swSigmaTheta(sal_all,temp_all,pres) # calculation of sigma (package oce)
  pres_sigma = NA
  pres_sigma = pres[which(!is.na(pres)==T & !is.na(sigma_all)==T)] # attribute pressure vector corresponding to the sigmainity values (with no NAs)
  sigma = NA
  sigma = sigma_all[which(!is.na(pres)==T & !is.na(sigma_all)==T)]# remove NAs from the sigmainity vector
  sigma = sigma[order(pres_sigma)] # order the sigma vector according to the increasing pressure
  pres_sigma = pres_sigma[order(pres_sigma)] # order the pressure vector corresponding to the sigmainity
  
  # calculate depth from pressure for the different depth vector (package oce) according to the latitude
  dep_temp = NA
  dep_temp = swDepth(pres_temp,lat)
  dep_sal = NA
  dep_sal = swDepth(pres_sal,lat)
  dep_sigma = NA
  dep_sigma = swDepth(pres_sigma,lat)
  
  # Zone calculation from Zone function
  zone = NA
  zone = Zone(lat,lon,temp,dep_temp)
  
  # Calculate of MLD (crit diff dty 0.03) from MLD_calc function
  MLD = NA
  MLD = MLD_calc(sigma,dep_sigma)
  
  
  ##############################
  ################### F) PAR RETRIEVAL (IF AVAILABLE)
  ##############################
  # Needed for the NPQ correction 
  
  if ("DOWNWELLING_PAR" %in% names(profile_B$var)==T) { # test if the variable is present in the netcdf file
    light_get = NA
    light_get = as.vector(ncvar_get(profile_B,"DOWNWELLING_PAR"))  #read the PAR variable as one unique vector
    
    # read the qc varaible as a vector
    light_qc = NULL
    for (ik in 1:dim(ncvar_get(profile_B,"DOWNWELLING_PAR_QC"))) {
      light_qc = paste(light_qc,ncvar_get(profile_B,"DOWNWELLING_PAR_QC")[ik],sep="")
    }
    
    # attribute NA to par values with bad qc
    light_all = NA
    light_all = light_get
    for (jj in 1:nchar(light_qc)) {
      if (substr(light_qc,jj,jj)==3 | substr(light_qc,jj,jj)==4) {
        light_all[jj] = NA
      }
    }
    
    pres_light = NA
    pres_light = pres[which(!is.na(pres)==T & !is.na(light_all)==T)] # attribute pressure vector corresponding to the light values (with no NAs)
    light = NA
    light = light_all[which(!is.na(pres)==T & !is.na(light_all)==T)]# remove NAs from the light vector
    light = light[order(pres_light)] # order the light vector according to the increasing pressure
    pres_light = pres_light[order(pres_light)] # order the pressure vector corresponding to the light
    
    # calculate depth from pressure (package oce) according to the latitude
    dep_light = NA
    dep_light = swDepth(pres_light,lat)
  }
  
  ###############################
  ######################### G) CHL RETRIEVAL AND TREATMENT 
  ###############################
  
  if ("CHLA" %in% names(profile_B$var)==F) { # test if the chla variable is present in the netcdf file (skip if not)
    print(error_message(106))
    nc_close(profile_C)
    nc_close(profile_B)
    return(106)
  }
  
  chl_get = NA
  chl_get = as.vector(ncvar_get(profile_B,"CHLA"))  #read the chla variable as one unique vector
  chl_all = NA
  chl_all = chl_get
  
  ############################
  ############ G.1) Gather meta information
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
  
  ### get the chla dark
  chla_calib_dark = chla_calib[grep("DARK_CHLA",chla_calib)] # get the dark information
  chla_calib_dark = unlist(strsplit(chla_calib_dark,"=")) # separate the name from the number
  chla_dark = as.numeric(chla_calib_dark[2]) # get the scale coefficient as a number
  
  nc_close(meta_nc) # close meta
  
  ############ 1) RANGE TEST : attribute NA to values out of range ###########################
  chl_all[which((chl_all > 50) | (chl_all < - 0.1))] = NA
  
  chl_not_isna = NA
  chl_not_isna = which(!is.na(pres)==T & !is.na(chl_all)==T) # keep the information on NA values to later add them
  
  pres_chl = NA
  pres_chl = pres[chl_not_isna] # attribute pressure vector corresponding to the chl values (with no NAs)
  chl = NA
  chl = chl_all[chl_not_isna]# remove NAs from the chl vector
  pres_chl_unsorted = NA
  pres_chl_unsorted = pres_chl # save the vector before it was sorted to allow for unsorting
  chl = chl[order(pres_chl)] # order the chl vector according to the increasing pressure
  pres_chl = pres_chl[order(pres_chl)] # order the pressure vector corresponding to the chl
  
  # calculate depth from pressure (package oce) according to the latitude
  dep_chl = NA
  dep_chl = swDepth(pres_chl,lat)
  
  # Test if the chl values are associated to only one depth (error of the measurement: "stuck pressure); if so, skip
  if(length(unique(dep_chl))==1) {
    print(error_message(107))
    nc_close(profile_C)
    nc_close(profile_B)
    return(107)
  }
  
  ################ 2) DARK OFFSET #######################
  # correct the vertical profile from a deep offset (Dark_Fchla_Corr function)
  chl_dark = NA
  chl_dark_offset = NA 
  chl_dark_min_pres = NA
  
  override_is_dmmc = (offset_override=="dmmc")
  if (length(override_is_dmmc)!=1) {override_is_dmmc=FALSE} #if offset override is NULL or is a vector, it is not =="dmmc"
  
  if (is.null(offset_override) | override_is_dmmc) { # if no override instruction is given or if instruction is to accept the offset from dmmc
      list_dark = Dark_Fchla_Corr(profile_actual, chl, dep_chl, MLD, zone, DEEP_EST)
      chl_dark = list_dark$chl_dark
      
      chl_dark_offset = list_dark$offset
      chl_dark_min_pres = list_dark$min_dep 
      if (!is.na(chl_dark_min_pres)) {
          chl_dark_min_pres = swPressure(chl_dark_min_pres, lat) # invert swDepth
      }
    
      if (!is.na(chl_dark_offset)) {
    	  factory_offset = chla_dark*chla_scale
    	  if ( abs(chl_dark_offset)>0.2*factory_offset & !override_is_dmmc) { # only check for error 112 if not instructed to ignore it
    	  	  print(error_message(112))
    		  return(112)
    	  }
      }
  } else { # if on override instruction is given
	  chl_dark_offset = offset_override[1]
	  chl_dark_min_pres = offset_override[2]
	  chl_dark = chl
	  if (!is.na(chl_dark_offset)) { chl_dark = chl_dark - chl_dark_offset }
	  if (!is.na(chl_dark_min_pres)) { chl_dark[which(dep_chl>=swDepth(chl_dark_min_pres,lat))] = 0 }
  }
  
  ############ 3) & 4) NPQ and FACTOR 2 ########################
  # Correct the chla profile from the:
  # - Non Photochemical Quenching: Xing et al., 2018
  # - Factor 2: Roesler et al., 2017
  
  #Pre-calculate RESO which will cause an error in the chl_npq computation if RESO=0
  RESO = NA
  RESO = mean(abs(diff(dep_chl[which(dep_chl<=250)])),na.rm=T) 
  if (RESO==0 | is.na(RESO)){
      print(error_message(108))
      nc_close(profile_C)
      nc_close(profile_B)
      return(108)
  }
  
  if ("DOWNWELLING_PAR" %in% names(profile_B$var)==T) { # test if there is a PAR in situ measured
    chl_npq_list = NPQ_cor_X12_XB18(chl_dark/2,dep_chl,dep_light,light,MLD)
  } else {
    chl_npq_list = NPQ_cor_P18(chl_dark/2,dep_chl,MLD)
  }
  
  chl_npq=chl_npq_list$chl_npq
  
  flag_NPQ_changed = NA
  flag_NPQ_changed = rep(FALSE, length(chl_npq))
  which_is_npq_changed = chl_npq_list$which_is_npq_changed
  flag_NPQ_changed[which_is_npq_changed] = TRUE # save the values that were changed by the NPQ correction to later write QC
 
  is_XB18 = chl_npq_list$is_XB18 
  
  ############################
  ############ H) BBP700 RETRIEVAL AND TREATMENT
  ############################
  
  bbp_get = NA
  bbp_get = as.vector(ncvar_get(profile_B,"BBP700")) #read the chla variable as one unique vector
  
  bbp_all = NA
  bbp_all = bbp_get
  
  ################### 1) RANGE TEST ###############################
  bbp_all[which((bbp_all > 0.1) | (bbp_all <  (-0.000025)))] = NA
  
  bbp_not_isna = NA
  bbp_not_isna = which(!is.na(pres)==T & !is.na(bbp_all)==T) # keep the information on NA values to later add them
  
  pres_bbp = NA
  pres_bbp = pres[bbp_not_isna] # attribute pressure vector corresponding to the bbp values (with no NAs)
  bbp = NA
  bbp = bbp_all[bbp_not_isna]# remove NAs from the bbp vector
  pres_bbp_unsorted = NA
  pres_bbp_unsorted = pres_bbp # save the vector before it was sorted to allow for unsorting
  bbp = bbp[order(pres_bbp)] # order the bbp vector according to the increasing pressure
  pres_bbp = pres_bbp[order(pres_bbp)] # order the pressure vector corresponding to the bbp
  
  # calculate depth from pressure (package oce) according to the latitude
  dep_bbp = NA
  dep_bbp = swDepth(pres_bbp,lat)
  
  
  #### 2) SPECIAL DRIFT CORRECTION (MANUALLY CALCULATED) #############################
  
  diff_bottom = NA
  
  if (substr(profile_actual,1,7)==5904218 & substr(profile_actual,1,11) < 660 & substr(profile_actual,1,11) > 456 ) {
    med_bottom = NULL
    med_bottom = median(bbp[which(dep_bbp < 300 & dep_bbp > 250)], na.rm=T)
    diff_bottom = 0-med_bottom
  }
  
  if (substr(profile_actual,1,7)==2902092 & substr(profile_actual,1,11) < 125 & substr(profile_actual,1,11) > 97 ) {
    med_bottom = NULL
    med_bottom = median(bbp[which(dep_bbp < 1600 & dep_bbp > 1400)], na.rm=T)
    diff_bottom = 0.000035-med_bottom
  }
  
  if (substr(profile_actual,1,7)==6901174 & format(time,"%Y-%m-%d", tz="UTC") >= as.POSIXct("06/03/2018",format="%d/%m/%Y", tz="UTC")) {
    med_bottom = NULL
    med_bottom = median(bbp[which(dep_bbp < 950 & dep_bbp > 850)], na.rm=T)
    diff_bottom = 0.0011-med_bottom
  }
  
  if (substr(profile_actual,1,7)==2902118 ) {
    med_bottom = NULL
    med_bottom = median(bbp[which(dep_bbp < 950 & dep_bbp > 850)], na.rm=T)
    diff_bottom = 0.00035-med_bottom
  }
  
  if (substr(profile_actual,1,7)==6901485 & format(time,"%Y-%m-%d", tz="UTC") >= as.POSIXct("08/05/2015",format="%d/%m/%Y", tz="UTC") &
      format(time,"%Y-%m-%d", tz="UTC")  <= as.POSIXct("28/07/2015",format="%d/%m/%Y", tz="UTC")) {
    med_bottom = NULL
    med_bottom = median(bbp[which(dep_bbp < 950 & dep_bbp > 850)], na.rm=T)
    diff_bottom = 0.0012-med_bottom
  }
  
  if (substr(profile_actual,1,7)==6901485 & format(time,"%Y-%m-%d", tz="UTC") >= as.POSIXct("01/09/2015",format="%d/%m/%Y", tz="UTC") &
      format(time,"%Y-%m-%d", tz="UTC")  <= as.POSIXct("04/03/2016",format="%d/%m/%Y", tz="UTC")) {
    med_bottom = NULL
    med_bottom = median(bbp[which(dep_bbp < 950 & dep_bbp > 850)], na.rm=T)
    diff_bottom = 0.0012-med_bottom
  }
  
  if (!is.na(diff_bottom)) {
    bbp = bbp + diff_bottom
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
  
  CHLA_QC = array(unlist(strsplit(ncvar_get(profile_B, "CHLA_QC"),"")), dim(chl_get_array))
  BBP700_QC = array(unlist(strsplit(ncvar_get(profile_B, "BBP700_QC"),"")), dim(bbp_get_array))
  
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
  
  chl_qc[which( CHLA_QC == "9")] = "9" # write 9 in place (missing values flagged)
  
  if (!is.na(chl_greylist_qc)){
      chl_qc[which(chl_qc=="2")] = chl_greylist_qc
  }
  
  ### bbp flags
  bbp_qc[which(!is.na(bbp_array))] = "1" # start with a 1 if a value exists
  bbp_qc[which( is.na(bbp_array) & !is.na(bbp_get_array) )] = "4" # write a "4" if the method removed a value
  bbp_qc[which( BBP700_QC == "9")] = "9" # write 9 in place (missing values flagged)
  
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
  ############ K.1) create error arrays 
  ############################
  
  chl_error = pmax( 3 * array(chla_scale, dim(chl_array)), 1.0 * abs(chl_array))
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
              "is_npq"=is_npq, "npq_depth"=npq_depth, "npq_val"=npq_val, "is_XB18"=is_XB18))
  
}


