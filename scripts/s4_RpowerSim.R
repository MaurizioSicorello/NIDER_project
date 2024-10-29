

################################
# Load packages and functions

library("here")
library("readxl")
library("poolr")
library("HiClimR")
library("RSpectra")
library("corrplot")
library("MASS")
library("clusterGeneration")
library("corpcor")
library("Hmisc")
library("ggplot2")
source(here("functions", "metaPoweR.R"))
df <- read_excel(here::here("data", "studyInformation", "studyInformation_publication.xlsx"))

# Load necessary package
library(Matrix)



################################
# test power function
metaPower(effectSize=0.05, sampleSizes=df$`sample size quest`)


################################
# number of effective tests

effectiveN <- data.frame(voxels=c(1, 10, 100, 1000, 10000), effective = rep(1,5))

effMethod = "liji"

# for "small" correlation matrices
for(i in 2:4){
  
  datTemp <- read.csv(here("results", "power", paste0("data", effectiveN[i,1], "voxels.csv")))
  effectiveN[i,2] <- meff(fastCor(datTemp),method = effMethod)
  
}


# for the largest correlation matrix
# takes 30 minutes to run. Resulting value was 23 for the galway method. 52 for Li
dat10000 <- read.csv(here("results", "power", "data10000voxels.csv"))
corr10000 <- fastCor(dat10000)

t1 <- Sys.time()
effectiveN[5, 2] <- meff(corr10000, method = effMethod)
t2 <- Sys.time()
t2-t1

# Galway
if(effMethod=="liji"){
  effectiveN[5,2] <- 52
}else if(effMethod=="galway"){
  effectiveN[5,2] <- 23
}




################################
# simulate data

dat100 <- read.csv(here("results", "power", "data100voxels.csv"))
dat100corr <- cor(dat100)

# matrix is not positive definite without changes
corrplot(dat100corr)
is.positive.definite(dat100corr)
dat100corr <- make.positive.definite(dat100corr)
dat100corr <- dat100corr[1:20, 1:20]
corrplot(dat100corr)

effectSize = 0.1
interCorr = 0.25
N = 250
iterations = 10000

numVars = nrow(dat100corr)


# add variable with preset correlations
new_matrix <- matrix(0, nrow = nrow(dat100corr) + 1, ncol = ncol(dat100corr) + 1)
new_matrix[1:nrow(dat100corr), 1:ncol(dat100corr)] <- dat100corr
new_matrix[nrow(new_matrix), 1:ncol(dat100corr)] <- rep(effectSize, numVars) # rep(c(effectSize, 0), each=numVars/2)
new_matrix[1:nrow(dat100corr), ncol(new_matrix)] <- rep(effectSize, numVars) # rep(c(effectSize, 0), each=numVars/2)
new_matrix[nrow(new_matrix), ncol(new_matrix)] <- 1
corrplot(new_matrix)

# output dataframe
powerSimOut = as.data.frame(
  matrix(nrow=iterations, ncol=6,
         dimnames = list(NULL, c("iterations", "sig", "powerNy", "powerLi", "powerGao", "powerGal"))))


set.seed(100)
# main loop
for(i in 1:iterations){
  
  # simulate data
  simDat <- mvrnorm(n = N, 
                    mu = rep(0, numVars+1),
                    Sigma = new_matrix)
  
  # test if significant after fdr
  pTemp <- rcorr(simDat)$P[ncol(simDat), 1:(ncol(simDat)-1)]
  powerSimOut[i,"sig"] <- ifelse(sum(p.adjust(pTemp, method = "fdr") <= 0.05) > 0, 1, 0) 
  
  # power based on effective number of tests from different methods
  nEffTemp <- meff(cor(simDat[, 1:(ncol(simDat)-1)]),method = "nyholt")
  powerSimOut[i,"powerNy"] <- metaPower(effectSize=effectSize, sampleSizes=N, tau = 0, alpha = 0.05/numVars, numTruePositives = nEffTemp)$fixed
  
  nEffTemp <- meff(cor(simDat[, 1:(ncol(simDat)-1)]),method = "liji")
  powerSimOut[i,"powerLi"] <- metaPower(effectSize=effectSize, sampleSizes=N, tau = 0, alpha = 0.05/numVars, numTruePositives = nEffTemp)$fixed
  
  nEffTemp <- meff(cor(simDat[, 1:(ncol(simDat)-1)]),method = "gao")
  powerSimOut[i,"powerGao"] <- metaPower(effectSize=effectSize, sampleSizes=N, tau = 0, alpha = 0.05/numVars, numTruePositives = nEffTemp)$fixed
  
  nEffTemp <- meff(cor(simDat[, 1:(ncol(simDat)-1)]),method = "galwey")
  powerSimOut[i,"powerGal"] <- metaPower(effectSize=effectSize, sampleSizes=N, tau = 0, alpha = 0.05/numVars, numTruePositives = nEffTemp)$fixed
}
  
sum(powerSimOut[,2])/iterations
colMeans(powerSimOut)


metaPower(effectSize=0.1, sampleSizes=250, tau = 0, alpha = 0.05/numVars, numTruePositives = meff(cor(dat100), method="nyholt"))

################################
# calculate power for different scenarios

# parameters
effVec <- c(seq(0,0.4,0.005))
tauVec <- c(0, 0.15, 0.3)
voxTruePos <- effectiveN$effective
voxFDR <- c(10000, 20000, 180000)

# output table
outTable <- expand.grid(effVec, tauVec, voxTruePos, voxFDR)
nConst <- nrow(outTable)
outTable <- data.frame(outTable, numeric(nConst))
colnames(outTable) <- c("effectSize", "tau", "numTruePos", "voxelFDR", "power")

# loop through parameter constellations (mapply probably also possible)
for(i in 1:nConst){
  
  outTable[i,"power"] <- metaPower(effectSize=outTable$effectSize[i], 
                                   sampleSizes=df$`sample size quest`, 
                                   tau = outTable$tau[i], 
                                   alpha = 0.05/outTable$voxelFDR[i], 
                                   numTruePositives = outTable$numTruePos[i])$random
  
}

# plot results
ggplot(data=outTable, aes(x=effectSize, y=power, linetype=as.factor(numTruePos))) +
  geom_line() +
  geom_hline(yintercept=0.8, color="blue") +
  facet_grid(tau ~ voxelFDR, scales = "free", labeller = labeller(tauVec = tau_labels, voxFDR = voxFDR_labels)) +
  labs(x="Correlation", y="Power", linetype="# True Positive Voxels") +
  scale_linetype_discrete(labels = c("1", "10", "100", "1000", "10000")) +
  
  theme_bw() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "bottom"
  ) 

ggsave(here("results", "figures", "powerPlot.png"))



metaPower(effectSize=0.125, sampleSizes=df$`sample size quest`, tau = 0, alpha = 0.05/20000, numTruePositives = 1)
