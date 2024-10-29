library("metafor")
library("MASS")
library("meta")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load example data and write to folder
dat <- dat.mcdaniel1994[complete.cases(dat.mcdaniel1994), ]
write.csv(dat, "mcdanielsTestData.csv", row.names = FALSE)

# run meta-analysis on first 5 rows (to assess decimal precision)
fit1 <- metacor(dat$ri[11:15], dat$ni[11:15], method.tau = "DL")
fit1

# download images used in Bossier Meta-Analysis


           






