library("MASS")

# relation between ERQ and psychopathology is around .3
# https://onlinelibrary.wiley.com/doi/10.1002/jclp.23206

# define correlation between simulated variables
correlXY = 0.3
correlThird = 0.3
N = 1000



# load function to correct variance restriction and unreliability
correctedCorrelations <- function(r, relX = 1, relY = 1, restrFactorX = 1, restrFactorY = 1) {
  # correct r for measurement error
  r_disatten <- r / sqrt(relX * relY)
  
  # correct r for restricted variance
  uX <- 1 / restrFactorX
  uY <- 1 / restrFactorY
  
  uT <- sqrt((relX * uX^2) / (1 + relX * uX^2 - uX^2))
  uP <- sqrt((relY * uY^2) / (1 + relY * uY^2 - uY^2))
  
  r_corrected <- r_disatten * uT * uP + sqrt((1 - uT^2) * (1 - uP^2))
  
  return(r_corrected)
}

# simulate data with defined correlation
simData = as.data.frame(mvrnorm(N, c(0,0, 0), Sigma = matrix(c(1,correlXY,correlThird,correlXY, 1, correlThird, correlThird, correlThird, 1), ncol=3)))
corrBefore <- cor(simData[,1], simData[,2])
sd(simData[,1])
sd(simData[,2])

# select subsample
simData_sub = simData[simData$V3 < quantile(simData$V3, 0.25), ]

corrAfter <- cor(simData_sub[,1], simData_sub[,2])

corrBefore
corrAfter

sd(simData[,1])
sd(simData[,2])

sd(simData_sub[,1])
sd(simData_sub[,2])

correctedCorrelations(corrAfter, restrFactorX=sd(simData_sub[,1])/sd(simData[,1]))
