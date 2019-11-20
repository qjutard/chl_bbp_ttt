############## MLD CALCULATION
#
#
# M.Cornec 20/11/19
#
# Calculate the MLD based on a difference in sigma (potential density) of 0.03 with reference to a "surface" value
# (de Boyer Montegut, 2003)
#
# Inputs:
# - sigma: vector of sigma values
# - dep_sigma: vector of depths assiociated to the sigma values
#
# Output: the value of the MLD depth
#


MLD_calc <- function(sigma,dep_sigma) { 
  MLD<-NA # set the MLD variable
  if(length(sigma[!is.na(sigma)==TRUE])>=2) { # test if there is at least 2 NA values of sigma to calculate the MLD
    sigmaSurface<-NA
    sigmaSurface <- approx(dep_sigma,sigma,10)$y # estimation of the reference density (approximation) at 10m (~"surface")
    if (is.na(sigmaSurface)==FALSE) { # test if there is a "surface" value
      MLD <- max(dep_sigma[sigma <= (sigmaSurface + 0.03)],na.rm = T) #identify the maximum depth where the sigma is <= to the surface value+0.03
    }
  }
  return(MLD[1])
}