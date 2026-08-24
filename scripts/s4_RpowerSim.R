

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
library("readr")
library("dplyr")
library("tidyr")
library("stringr")
library("scales")
source(here("functions", "metaPoweR.R"))

df <- read_excel(here::here("data", "studyInformation", "studyInformation_publication.xlsx"))

# Load necessary package
library(Matrix)



################################
# test power function
metaPower(effectSize=0.10, tau=0.06, sampleSizes=df$`sample size quest`)


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


dat20000 <- read.csv(here("results", "power", "data20000voxels.csv"))
corr20000 <- fastCor(dat20000)

t1 <- Sys.time()
effectiveN_fullNet <- meff(corr20000, method = effMethod)
t2 <- Sys.time()
t2-t1
# effective N is 53!


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






#########################################
# plot cluster-insertion simulation


# ------------------------------------------------------------------
# file names
# ------------------------------------------------------------------
file_t_006 <- here("results", "power", "simulation_results_r0.10_t0.06.csv")
file_t_008 <- here("results", "power", "simulation_results_r0.10_t0.08.csv")

# ------------------------------------------------------------------
# helper
# ------------------------------------------------------------------
read_power_run <- function(file, tau_label) {
  
  dat <- read_csv(file, show_col_types = FALSE)
  
  tibble(
    scenario = factor(c("1", "2", "3", "4"), levels = c("1", "2", "3", "4")),
    voxels_total = c("800", "1300", "1600", "1800"),
    ROI = c(mean(dat$pow1_reg), mean(dat$pow2_reg), mean(dat$pow3_reg), mean(dat$pow4_reg)),
    Whole_brain = c(mean(dat$pow1_all), mean(dat$pow2_all), mean(dat$pow3_all), mean(dat$pow4_all))
  ) %>%
    pivot_longer(
      cols = c(ROI, Whole_brain),
      names_to = "scope",
      values_to = "power"
    ) %>%
    mutate(
      tau = tau_label,
      x_lab = factor(
        paste0("Clusters: ", scenario, "\nVoxels: ", voxels_total),
        levels = paste0(
          "Clusters: ", c("1", "2", "3", "4"),
          "\nVoxels: ", c("800", "1300", "1600", "1800")
        )
      )
    )
}

plot_dat <- bind_rows(
  read_power_run(file_t_006, "τ = .06"),
  read_power_run(file_t_008, "τ = .08")
)

# ------------------------------------------------------------------
# plot
# ------------------------------------------------------------------
p <- ggplot(
  plot_dat,
  aes(
    x = x_lab,
    y = power,
    color = scope,
    linetype = tau,
    group = interaction(scope, tau)
  )
) +
  geom_hline(yintercept = .80, linewidth = 0.5, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = label_percent(accuracy = 1)
  ) +
  scale_color_manual(
    values = c("ROI" = "#1b9e77", "Whole_brain" = "#d95f02"),
    labels = c("Network of Interest", "Whole-brain")
  ) +
  scale_linetype_manual(
    values = c("τ = .06" = "solid", "τ = .08" = "longdash")
  ) +
  labs(
    x = NULL,
    y = "Power",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(size = 10, lineheight = 0.95),
    legend.position = "right"
  ) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.08), legend.justification = c(1, 0)) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2,
                            override.aes = list(color = "black", linewidth = 1.4),
                            keywidth = unit(1.8, "cm"))
  )

p

ggsave(
  filename = here("results", "figures", "power_simulation_cluster.png"),
  plot = p,
  width = 7,
  height = 4.5,
  dpi = 300
)
