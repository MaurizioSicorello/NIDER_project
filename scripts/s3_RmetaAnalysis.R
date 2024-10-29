##############################################################################
#### Meta-Analysis of congruence between emotion regulation outcomes and conduct power analyses
##############################################################################
#To install the correct packages and reproduce the code, run:
# install.packages("renv")
# library("renv")
# renv::restore()


#######################################
# load packages, data and custom function

library("renv")
library("here")
library("readxl")
library("stringr")
library("meta")
library("WebPower")
library("eulerr")
library("flextable")


df <- read_excel(here::here("data", "studyInformation", "studyInformation_publication.xlsx"))
dfcorrs <- df[,c("image name stem", "sample size quest", "sample size rating", "sample size amygdala", "amyQuestCorr", "questSelfRating", "amySelfRating", "questionnaire name")]

extractMAinfo <- function(model){
  
  # Extract the effect size (r) and its confidence interval
  r <- model$TE.random
  r_ci_lower <- model$lower.random
  r_ci_upper <- model$upper.random
  p_value <- model$pval.random
  k <- model$k
  N <- sum(model$data$.n)
  
  # Extract tau and its confidence interval
  tau <- model$tau
  tau_ci_lower <- model$lower.tau
  tau_ci_upper <- model$upper.tau
  
  # Extract Q statistic, degrees of freedom, and p-value for heterogeneity
  Q <- model$Q
  Q_df <- model$df.Q
  Q_p <- model$pval.Q
  
  #out <- c(round(r, 2), round(r_ci_lower, 2), round(r_ci_upper, 2), round(p_value, 3), round(k,1), round(N,1), round(tau, 2), round(tau_ci_lower, 2), round(tau_ci_upper, 2), round(Q_df,1), round(Q, 2), round(Q_p, 3))
  # Create a vector with specified formatting, including confidence intervals in [lower, upper] format
  out <- c(
    sprintf("%.2f", r),
    sprintf("[%.2f, %.2f]", r_ci_lower, r_ci_upper),
    sprintf("%.3f", p_value),
    as.character(k),
    as.character(N),
    sprintf("%.2f", tau),
    sprintf("[%.2f, %.2f]", tau_ci_lower, tau_ci_upper),
    sprintf("%.2f", Q),
    as.character(Q_df),
    sprintf("%.3f", Q_p)
  )
  
  return(out)
}


#######################################
# average morawetz 2016b

morawetzRows <- dfcorrs[dfcorrs$`image name stem` == "Morawetz2016b", ][1:2, 5:7]
dfcorrs[which(dfcorrs$`image name stem` == "Morawetz2016b")[1], 5:7] <- as.list(apply(morawetzRows, 2, function(x) {as.character(mean(as.numeric(x)))}))
dfcorrs <- dfcorrs[-which(dfcorrs$`image name stem` == "Morawetz2016b")[2], ]

onlineRatings <- !dfcorrs$`image name stem` %in% c("Diers2021", "GaertnerDiersUnpublished", "HofhanselUnpublished")

#######################################
# Questionnaires and self reports

# subset df
dfQuestSelf <- dfcorrs[!is.na(dfcorrs$questSelfRating) & onlineRatings, ]
dfQuestSelf$`image name stem`

# perform MA
meta_questSelf <- metacor(as.numeric(dfQuestSelf$questSelfRating), as.numeric(dfQuestSelf$'sample size quest'), method.tau = "DL", studlab = dfQuestSelf$`image name stem`)
meta_questSelf

# plot results
meta::forest(meta_questSelf)
png(file=here("results", "figures", "forestplot_questSelf.png"), width = 800, height = 600)
meta::forest(meta_questSelf)
dev.off() 




#######################################
# Questionnaires and amygdala

# subset df
dfAmyQuest <- dfcorrs[!is.na(dfcorrs$amyQuestCorr), ]

# perform MA
meta_amyQuest <- metacor(as.numeric(dfAmyQuest$amyQuestCorr), dfAmyQuest$'sample size quest', method.tau = "DL", studlab = dfAmyQuest$`image name stem`)
meta_amyQuest

# plot results
meta::forest(meta_amyQuest)
png(file=here("results", "figures", "forestplot_amyQuest.png"), width = 800, height = 700)
meta::forest(meta_amyQuest)
dev.off() 



#######################################
# amygdala and self reports

# subset df
dfAmySelf <- dfcorrs[!is.na(dfcorrs$amySelfRating) & onlineRatings, ]

# perform MA
meta_amySelf <- metacor(as.numeric(dfAmySelf$amySelfRating), dfAmySelf$'sample size rating', method.tau = "DL", studlab = dfAmySelf$`image name stem`)
meta_amySelf

# plot results
meta::forest(meta_amySelf)
png(file=here("results", "figures", "forestplot_meta_amySelf.png"), width = 800, height = 600)
meta::forest(meta_amySelf)
dev.off() 

# calculate average power
wp.correlation(r = meta_amySelf$TE.random, power = 0.8)
mean(sapply(meta_amySelf$n, function(x){wp.correlation(n = x, meta_amySelf$TE.random)$power}))
sd(sapply(meta_amySelf$n, function(x){wp.correlation(n = x, meta_amySelf$TE.random)$power}))
median(meta_amySelf$n)



#######################################
# results table

tableColNames <- c("Model", "r", "95% CIa", "pa", "k", "N", "tau", "95% CIb", "Q", "df", "pb")

outTable <- data.frame(rbind(extractMAinfo(meta_amySelf), extractMAinfo(meta_questSelf), extractMAinfo(meta_amyQuest)))
outTable <- cbind(c("Rating-Amygdala", "Rating-Questionnaires", "Amygdala-Questionnaires"), outTable)
names(outTable) <- tableColNames

# flextable::save_as_docx(
#   flextable(
#     outTable
#   ),
#   path = here("results", "tables", "sharedMAtable.docx")
# )


#######################################
# overlap plot (Y = quest, X1 = selfreport, X2 = amygdala)
# https://www.andrewheiss.com/blog/2021/08/21/r2-euler/
# https://stats.stackexchange.com/questions/314926/can-you-calculate-r2-from-correlation-coefficents-in-multiple-linear-regressi/364630#364630

ADEG = 1
BDFG = 1
CEFG = 1
DG = meta_questSelf$TE.random^2
EG = meta_amyQuest$TE.random^2
FG = meta_amySelf$TE.random^2
A = 1 - t(c(meta_questSelf$TE.random, meta_amyQuest$TE.random)) %*% solve(matrix(c(1,meta_amySelf$TE.random, meta_amySelf$TE.random, 1), nrow=2)) %*% c(meta_questSelf$TE.random, meta_amyQuest$TE.random)
B = 1 - t(c(meta_questSelf$TE.random, meta_amySelf$TE.random)) %*% solve(matrix(c(1,meta_amyQuest$TE.random, meta_amyQuest$TE.random, 1), nrow=2)) %*% c(meta_questSelf$TE.random, meta_amySelf$TE.random)
C = 1 - t(c(meta_amyQuest$TE.random, meta_amySelf$TE.random)) %*% solve(matrix(c(1,meta_questSelf$TE.random, meta_questSelf$TE.random, 1), nrow=2)) %*% c(meta_amyQuest$TE.random, meta_amySelf$TE.random)
D = ADEG - A - EG
E = CEFG - C - FG
F = BDFG - B - DG
G = ifelse((DG - D) < 0, 0, (DF - D))

defineOverlap = c(
  "Questionnaire" = 100,  # Total variance for Questionnaire
  "Self-Rating" = 100,    # Total variance for Self-Rating
  "Amygdala" = 100,       # Total variance for Amygdala
  "Questionnaire&Self-Rating" = round(DG*100,2),
  "Questionnaire&Amygdala" = round(EG*100,3),# Overlap between Questionnaire and Self-Rating     # Overlap between Questionnaire and Amygdala
  "Self-Rating&Amygdala" = round(FG*100,2) # Overlap between Self-Rating and Amygdala
)

fit <- euler(defineOverlap)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

pdf(here("results", "figures", "venn_diagram.pdf"))
plot(fit, fills = myCol, labels = list(fontsize=16))
dev.off()




r = 0.5
ut = 1/1.2
up = 1/1.2

r*ut*up + sqrt((1-ut^2)*(1-up^2))

test = data.frame(r = c(0.5, 0.3, 0.1), n = c(50, 70, 90), this = c(1.5, 1.5, 1.5))

summary(ma_r(rxyi = r, n = n, data = test, ))

