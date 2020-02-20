# DARK OFFSET CALCULATION IN CASES OF DEEP VERTICAL MIXING
#
# M.Cornec 20/11/19
#
# Calculate a specific offset for the profiles where a deep vertical mixing leads to non-null deep chla values
# Those profiles are identified as showing an MLD deeper than 600m and an anomaly of the deep values of chla regardinf to the time serie of the float
# A dark offset is then estimated from a robust linear regression bewteen the time series of deep chla values as function of the number of profiles/time of the time serie for the concerned profiles.
#
# Inputs:
# - WMO : the WMO number of the float (7 digits)
# - path_to_netcdf: the path to the netcdf files 
# - index_ifremer: the merge index file
#
# Output: a dataframe with the reference of the profiles needing deep vertical mixing correction (WMO number of the float + profile number,format: "XXXXXXX_XXX"), 
#         and the assiocated dark offsets. If no profiles need correction, the dataframe is returned empty (NULL). 




Dark_MLD_table_coriolis <- function (WMO,path_to_netcdf,index_ifremer) {
  
  
  DEEP_EST = NULL # initiate the output
  
  # Read the index_ifremer file
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
  
  float_list = unique(wod)
  
  deep_table = NULL # set the time serie dataframe
  
  
  for(i in files[which(wod==WMO)])  { #loop on the profiles corresponding to the WMO
    
    # skip the profile if no "CHLA" is measured 
    if ("CHLA" %in% variables[[which(files==i)]] == F) {
      next
    }  
    
    iii<-NA
    iii<-prof_id[which(files==i)] # read the profile reference
    
    # skip if the profile is a descent one
    if(substr(iii,14,14)=="D") {
      next
    } 
    
    file_B = NA
    file_B = paste(path_to_netcdf, i, sep="") 
    
    path_split = NA
    path_split = unlist( strsplit(i, "/") )
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
    ############# POSITION : LON / LAT / DATE
    #################
    
    position_qc<-NA
    position_qc<-substr(ncvar_get(profile_B,"POSITION_QC"),1,1) # read position QC
    
    # skip the profile if the position QC is bad
    if (position_qc == 3 | position_qc==4) {
        #nc_close(profile)
        nc_close(profile_C)
        nc_close(profile_B)
        next
    }
    
    lat<-NA
    lat<- ncvar_get(profile_B,"LATITUDE")[1]
    lon<-NA
    lon<- ncvar_get(profile_B,"LONGITUDE")[1]
    
    # skip the profile if one (or both) coordinate(s) is(are) missing
    if (is.na(lat) | is.na(lon)) {
        #nc_close(profile)
        nc_close(profile_C)
        nc_close(profile_B)
        next
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
        #nc_close(profile)
        nc_close(profile_C)
        nc_close(profile_B)
        next
    }
    
    # skip the profile if julian date qc is bad
    if (jd_qc == 3 | jd_qc==4) {
        #nc_close(profile)
        nc_close(profile_C)
        nc_close(profile_B)
        next
    }
    
    #########################
    ############## DEF EMPTY DTFR
    ########################
    
    profile_dark<-NA
    profile_dark<-data.frame(profile_dark) # set the empty dataframe
    profile_dark$WMO<-NA
    profile_dark$JULD<-NA
    profile_dark$TIME<-NA
    profile_dark$MLD<-NA
    profile_dark$dark_deep<-NA
    profile_dark$profile<-NA
    
    profile_dark$WMO<-substr(iii,3,13) # add the profile id (wmo number+profile number)
    profile_dark$profile<-as.numeric(paste(substr(iii,11,13))) # add the profile number only
    profile_dark$TIME<-time # add the time (posixct format)
    profile_dark$JULD<-jd #add the time (julian day format)
    
    ###################
    ############# PHYSIQUE : TEMP / SAL / DTY / DEPTH / MLD
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
    
    # Calculate of MLD (crition of density of 0.03 (de Boyer Montegut 2004))
    MLD<-NA
    MLD<-MLD_calc(sigma,dep_sigma)
    
    profile_dark$MLD<-MLD # add the MLD in the profile
    
    ###############################
    ######################### CHL 
    ###############################
    
    if ("CHLA" %in% names(profile_B$var)==F) { # test if the chla variable is present in the netcdf file (skip if not)
        #nc_close(profile)
        nc_close(profile_C)
        nc_close(profile_B)
        next
    }
    
    chl_get<-NA
    chl_get <- as.vector(ncvar_get(profile_B,"CHLA"))  #read the PAR variable as one unique vector
    chl_all<-NA
    chl_all<-chl_get
    
    
    # Range test : attribute NA to values out of range
    chl_all[which((chl_all > 50) | (chl_all < - 0.1))]<-NA
    
    pres_chl<-NA
    pres_chl<-pres[which(!is.na(pres)==T & !is.na(chl_all)==T)] # attribute pressure vector corresponding to the chl values (with no NAs)
    chl<-NA
    chl<-chl_all[which(!is.na(pres)==T & !is.na(chl_all)==T)]# remove NAs from the chl vector
    chl<-chl[order(pres_chl)] # order the chl vector according to the increasing pressure
    pres_chl<-pres_chl[order(pres_chl)] # order the pressure vector corresponding to the chl
    
    # calculate depth from pressure (package oce) according to the latitude
    dep_chl<-NA
    dep_chl = swDepth(pres_chl,lat)
    
    # Test if the chl values are associated to only one depth (error of the measurement: "stuck pressure); if so, skip
    if(length(unique(dep_chl))==1) {
        #nc_close(profile)
        nc_close(profile_C)
        nc_close(profile_B)
        next
    }
    
    
    ###################################
    ######################### DARK OFFSET CALCULATION
    ###################################
    
    if(is.na(MLD)==F) { # test if there is an MLD value
      if (max(dep_chl,na.rm=T)>900) { # test if the chla profile present values deeper than 900m
        # Calculate the median of the chla values between the max depth and 50m above the max depth
        dark_deep<-NA
        dark_deep<-median(chl[which(dep_chl<=max(dep_chl,na.rm=T) & 
                                      dep_chl>= (max(dep_chl,na.rm=T)-50))],na.rm=T)
        
        profile_dark$dark_deep<-dark_deep #add the dark_deep median value in the profile dataframe 
        deep_table<-rbind(deep_table,profile_dark) # bind the profile dataframe to the time serie dataframe
      }
    }
    #nc_close(profile)
    nc_close(profile_C)
    nc_close(profile_B)
  }
  
  ###################################
  ######################### DEEP VERTICAL MIXING CASES
  ###################################
  
  if(all(is.na(deep_table$dark_deep))==F) { # test if the time serie dataframe present non NA values
    if(length(deep_table[,1]) > 5) { # test if there is at least 5 deep_dark profiles values
      
      # add a column in the dataframe indicating if there is profile with an MLD deeper than 600m
      deep_table$MLD_deep<-NA
      deep_table$MLD_deep[which(deep_table$MLD > 600)]<-1 
      
      # test if there is an outliar in the median deep chla values time serie (Outliars_med function)
      deep_table$flag_out<-NA
      deep_table$flag_out<-Outliars_med(deep_table$dark_deep,deep_table$profile)
      
      # retrieve the coefficients from a robust linear regression of the time serie of deep values as function of the profile number
      coeff_corr<-NA
      coeff_corr<-rlm(deep_table$dark_deep~deep_table$profile)$coef[2]
      offset_corr<-NA
      offset_corr<-rlm(deep_table$dark_deep~deep_table$profile)$coef[1]
      
      # Test if a profile shows both a deep outliar regarding to the time serie, and a deep MLD
      deep_table$corr_need<-0
      deep_table$corr_need[which(deep_table$flag_out==1 & deep_table$MLD_deep ==1 )]<-1
      
      
      if(length(deep_table[which(deep_table$corr_need==1),1]) > 0) { # test one or more profiles need the deep mixing correction
        for (kk in deep_table$profile[which(deep_table$corr_need==1)]) {
          
          # calculate the deep mixing offset fot the concerned profile (from the robust linear regression coefficients)
          dark_est<-NA
          dark_est<-(kk*coeff_corr)+offset_corr
          
          sub_ref<-NA
          sub_ref<-deep_table$WMO[which(deep_table$profile==kk)] #select the WMO reference of the profile (WMO + profile number, 11 elements, format: "XXXXXXX_XXX")
          mer_est<-NA
          mer_est<-as.data.frame(cbind(sub_ref,dark_est)) # create a dataframe with the WMO reference and the dark deep mixing coefficient
          DEEP_EST<-rbind(DEEP_EST,mer_est) # bind all the profiles needing the dark ddep vertical mixing correction
        }
      }
      
    }
  }
  return(DEEP_EST)
}