## OUTLIARS DETECTION
# 
# M.Cornec 20/11/19
#
# Detect and attribute a flag to outliars in a serie of values
# the criterion is twice the quantile 85 (tested on the residus from the substraction of the baseline on the serie of values)
# Baseline is calculted as a robust linear regression on the serie of values as function of the profile number
#
# Inputs:
# - Varia: vector of the serie of values
# - Leng: vector of index corresponding to the serie of values (profile number or time)
#
# Output: a binary vector (0,1) with same length of the input, 1=outliar 
#
#

Outliars_med<-function(Varia,Leng) {
  a<-NA
  b<-NA
  
  a<-rlm(Varia~Leng)$coef[2] # retrieve the slope of the robust linear regression of the serie of values as function of the index
  b<-rlm(Varia~Leng)$coef[1] # retrieve the offset of the robust linear regression of the serie of values as function of the index

  
  baseline<-NA
  baseline<-(Leng*a)+b # calculate the baseline from the robust linear regression coefficients
  
  res_var<-Leng*0
  res_var<-Varia-baseline # calculate the residus from the substraction of the baseline to the serie of values
  Q85<-NA
  Q85<-rep(2*quantile(res_var,0.85,na.rm = T),length(Leng)) # create a vector with the criterion value
  flag_out<-Varia*0 # create default vector of 0 (no outliar)
  flag_out[which(res_var>Q85)]<-1 # attribute a value of one the residus detected as outliars
  
  return(flag_out)
  
}