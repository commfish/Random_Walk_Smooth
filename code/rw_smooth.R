# Random Walk Smooth for Survey Estimates
# ---------------------------------------
# Author: William Gaeuman
# Contact: william.gaeuman@alaska.gov
# Last updated: 2017-11-02
# 
# Assumptions:
# 1) log observed survey estimate = log true survey value + N(0, sqrt(log(1 + cv^2));
# 2) true log survey value in year t = true (hypothetical) log survey value in year t - 1 + N(0, sigma), 
# where sigma is the "process error".
#
# Uses function optim() with assumed process error to minimize negative log likelihood wrt true survey values.
# Returns a list whose components include results of optim() and all function inputs.
# 
# Function inputs: 
# s.yr		= survey years
# b		= observed survey value, i.e. abunance/biomass
# cv.b		= survey cvs
# sigma		= assumed process error (random walk standard deviation on log scale)
#
# Defaults are PIGKC total biomass estimates from "biennial" NMFS shelf slope survey.
# ----

rw.smooth <- function(s.yr = c(2002, 2004, 2008, 2010, 2012, 2016),
	b = c(816, 989, 1216, 1776, 1444, 1464),
	cv.b = c(.19, .32, .26, .29, .35, .26), 
	sigma = 0.15){

n <- s.yr[length(s.yr)] - s.yr[1] + 1 		# number of years, e.g. 15
s.index <- s.yr - s.yr[1] + 1			# survey year index, e.g. 1, 3, 7, 9, 11, 15	

# Negative log likelihood function wrt predicted survey abundance/biomass b.hat
negll <- function(b.hat, sigma, b, cv.b, n, s.index){
	var.log.b <- log(1 + cv.b^2)
	sigma2 <- sigma^2
	
	# survey contribution
	value <- sum(log(2 * pi * var.log.b) + (log(b) - log(b.hat[s.index]))^2 / var.log.b)
	
	# random walk contribution
	for(t in 2:n) value <- value + log(2 * pi * sigma2) + (log(b.hat[t]) - log(b.hat[t - 1]))^2 / sigma2
	
	return(0.5 * value)
	}

# Starting parameter values
b.hat.init <- rep(mean(b), n)
b.hat.init[s.index] <- b

# Minimization
out <- optim(b.hat.init, negll, gr = NULL, sigma = sigma, b = b, cv.b = cv.b, n = n, s.index = s.index, 
	method = c("L-BFGS-B"), hessian = TRUE)

return(list(optim = out, yr = yr,  s.yr = s.yr, b = b, cv.b = cv.b, sigma = sigma))
}




