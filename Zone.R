############################### ZONE DELIMITATION 
# 
# 
# M.Cornec 20/11/19
#
# Attribute a regional qualification to the profile depending on the GPS position
# and/or the temperature profile
#
# BAFF: Baffin Bay
# NASPG: North Atlantic SubPolar Gyre
# NPSPG: North Pacific SubPolar Gyre
# ARCT: Arctic waters 
# NS: Northern Seas (Groenland, Barents and Norway seas)
# NASTG: North Atlantic SubTropical Gyre
# NAC: North Atlantic Current
# KURIO: Kuroshio Current
# NPSTG: North Pacific SubTropical Gyre
# ASEW: Atlantic SubEquatorial Waters
# PSEW: Pacific SubEquatorial Waters
# SPSTG: South Pacific SubTropical Gyre
# SASTG: South Atlantic SubTropical Gyre
# SISTG: South Indian SubTropical Gyre
# IEQ: Indian Equatorial waters
# IOMZ: Indian Oxygen Minimum Zones (Arabian Sea and Bengal Bay)
# MOONS: Mounsoon region
# UPW: Upwelling Zone
# WMS & EMS: Western & Eastern Mediterranean Sea
# BKS: Black Sea
# RDS: Red Sea
# CHINS: Chineese Sea
# AUS: Australian waters
# CAL: Californian current
# ARCH: Archipelagos waters
# STZ: Sub Tropical Zone
# SAZ: Sub Antarctic Zone 
# PFZ: Polar Frontal Zone
# ASZ_SIZ: Seasonal Ice Zone
#
#  The Southern ocean delimitation is based on temperature profile from Gray et al., 2018.
#  Temperature is also used to delineate the Kuroshio Current, the North Atlantic Current,
#  and the Arctic waters.
#
# Inputs:
# - lon: longitude
# - lat: latitude
# - temp: vector of temperature values the vertical profile
# - dep_temp: vector of depth values associated to the temperature vector
#
# Output: character variable with the acronym of the region 
#

Zone<-function(lat,lon,temp,dep_temp) {
  zone<-NA
  
  if(lat > 45 & lon < -100  ) {
    zone<-"NPSPG"
  }
  
  if(lat > 60 & lat < 80 & lon > -12.5  & lon < 20  ) {
    zone<-"NS"
  }
  
  if(lat > 17.5 & lat < 32.5 & lon > -75  & lon < -20  ) {
    zone<-"NASTG"
  }
  
  if(lat > 17.5 & lat < 32.5 & lon > -175  & lon < -140  ) {
    zone<-"NPSTG"
  }
  
  if(lat > -30 & lat < -14 & lon > -180  & lon < -90  ) {
    zone<-"SPSTG"
  }
  
  if(lat > -14 & lat < 6 & lon > -125  & lon < -75  ) {
    zone<-"UPW"
  }
  
  if(lat > 0 & lat < 17.5 & lon > -35  & lon < -10  ) {
    zone<-"ASEW"
  }
  
  if(lat > 0 & lat < 17.5 & lon > -140  & lon < -100  ) {
    zone<-"PSEW"
  }
  
  if(lat > -30 & lat < -15 & lon > -45  & lon < 0  ) {
    zone<-"SASTG"
  }
  
  if(lat > -30 & lat < -15 & lon > 50  & lon < 110  ) {
    zone<-"SISTG"
  }
  
  if(lat > -5 & lat < 7.5 & lon > 40  & lon < 100  ) {
    zone<-"IEQ"
  }
  
  if(lat > -15 & lat < -5 & lon > 50  & lon < 110  ) {
    zone<-"MOONS"
  }
  
  if(lat >= 7.5 & lat < 24 & lon > 43  & lon < 100  ) {
    zone<-"IOMZ"
  }
  
  if(lat > 35 & lat < 43 & lon > -5  & lon < 6  ) {
    zone<-"WMS"
  }
  
  if(lat > 35 & lat < 45 & lon > 5  & lon < 16 ) {
    zone<-"WMS"
  }
  
  if(lat > 30 & lat < 45 & lon > 16  & lon < 22.5  ) {
    zone<-"EMS"
  }
  
  if(lat > 30 & lat < 40 & lon > 22.5  & lon < 35  ) {
    zone<-"EMS"
  }
  
  if(lat > 41 & lat < 46 & lon > 12  & lon < 15.5  ) {
    zone<-"EMS"
  }
  
  if(lat > 30 & lat < 37.5 & lon > 10  & lon < 15.5  ) {
    zone<-"EMS"
  }
  
  if(lat > 40.98 & lat < 47.35 & lon > 27.67  & lon < 42.5  ) {
    zone<-"BKS"
  }
  
  if(lat > 10.78 & lat < 30.86 & lon > 32.78  & lon < 42.4  ) {
    zone<-"RDS"
  }
  
  if(lat > -30 & lat < -7 & lon > 110  & lon < 160  ) {
    zone<-"AUS"
  }
  
  if(lat > 0 & lat < 24 & lon > 100  & lon < 120  ) {
    zone<-"CHINS"
  }
  
  if(lat > -30 & lat < -5 & lon > 142  & lon < 180  ) {
    zone<-"ARCH"
  }
  
  if(lat > 30 & lat < 50 & lon > -135  & lon < -115  ) {
    zone<-"CAL"
  }
  
  if(lat > 0 & lat < 17.5 & lon > -140  & lon < -100  ) {
    zone<-"PSEW"
  }
  
  if (lat < 65 & lat > 30 & lon < -10 & lon > -80) { # North Atlantic delimitation
    
    if(length(dep_temp) > 5 & 
       max(dep_temp,na.rm=T) >= 100 & 
       min(dep_temp,na.rm=T) <= 100) {
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      if (t_100 > 17.5) {
        zone<-"NASTG"
      }
      
      if (t_100 < 17.5 & t_100 > 8 ) {
        zone<-"NAC"
      }
      
      if (t_100 < 8 ) {
        zone<-"NASPG"
      }
    }
  }
  
  if (lat > 65 ) { # Arctic delimitation
    
    if(length(dep_temp) > 5 & 
       max(dep_temp,na.rm=T) >= 100 & 
       min(dep_temp,na.rm=T) <= 100) {
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      if (t_100 < 3) {
        zone<-"ARCT"
      }
    }
  }
  
  if (lat < 65 & lat > 0 & lon > 120) { # North Pacific delimitation
    
    if(length(dep_temp) > 5 & 
       max(dep_temp,na.rm=T) >= 100 & 
       min(dep_temp,na.rm=T) <= 100) {
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      if (t_100 > 15) {
        zone<-"NPSTG"
      }
      
      if (t_100 < 15 & t_100 > 7.5 ) {
        zone<-"KURIO"
      }
      
      if (t_100 < 7.5 ) {
        zone<-"NPSPG"
      }
    }
  }
  
  
  if(lat > 65 & lat < 70 & lon > -72  & lon < -50  ) {
    zone<-"BAFF"
  }
  
  if(lat > 70 & lat < 80 & lon > -80  & lon < -50  ) {
    zone<-"BAFF"
  }
  
  if(is.na(zone) &  lat > -30 & lat < -25 & lon > 0  & lon < 50  ) {
    zone<-"STZ"
  }
  
  
  if(lat < - 30) { # Southern Ocean temperature delimitation
    
    if (max(dep_temp,na.rm=T)>=400 & min(dep_temp,na.rm=T)<=100) {
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      t_400<-NA
      t_400<-approx(dep_temp,temp,400)$y
      
      t_0_200<-NA
      t_0_200<-min(temp[which(dep_temp<=200)],na.rm=T)
      
      if (t_100 >= 11) {
        zone<-"STZ"
      }
      
      if (t_100 < 11 & t_400 >= 5) {
        zone<-"SAZ"
      }
      
      if (t_400 < 5 & t_0_200 >= 2) {
        zone<-"PFZ"
      }
      
      if (t_0_200 < 2) {
        zone<-"ASZ_SIZ"
      }
    }
  }
  
  if(lat > -30 & lat < -14 & lon > -180  & lon < -90  ) {
    zone<-"SPSTG"
  }
  
  return(zone)
}