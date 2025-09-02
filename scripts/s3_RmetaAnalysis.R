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
library("dplyr")
library("meta")
library("metafor")
library("WebPower")
library("eulerr")
library("flextable")
library("metaBMA")
library("ggplot2")
source(here("functions", "metaPoweR.R"))


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

dfcorrs <- dfcorrs %>%
  dplyr::mutate(formatted_author = if_else(
    str_detect(`image name stem`, "\\d"),
    str_replace(`image name stem`, "^([A-Za-z\\+]+)(\\d.*)", "\\1 et al. (\\2)"),
    `image name stem`  
  ))

# Create a named vector for mapping
author_mapping <- c(
  "Benzait et al. (2023)" = "Benzait et al. (2023)",
  "Berboth et al. (2021)" = "Berboth et al. (2021)",
  "Brehl et al. (2020)" = "Brehl et al. (2021; partially unpublished)",
  "Burghart+SchmidtUnpublished" = "Burghart et al. (unpublished)",
  "Diers et al. (2021)" = "Diers, Gärtner et al. (2023)",
  "Doerfel et al. (2014a)" = "Dörfel et al. (2014; distancing)",
  "Doerfel et al. (2014b)" = "Dörfel et al. (2014; reinterpretation)",
  "Gaebler et al. (2014)" = "Gaebler et al. (2014)",
  "GaertnerDiersUnpublished" = "Scheffel et al. (2019)",
  "Gianaros et al. (2020-AHAB-II)" = "Gianaros et al. (2020; AHAB-II)",
  "Gianaros et al. (2020-PIP)" = "Gianaros et al. (2020; PIP)",
  "Glosemeyer et al. (2020)" = "Glosemeyer et al. (2020)",
  "HofhanselUnpublished" = "Hofhansel et al. (2023)",
  "Jentsch et al. (2019)" = "Jentsch et al. (2019)",
  "PowersUnpublished" = "Powers et al. (2022)",
  "WessaUnpublished" = "Wessa et al. (unpublished)",
  "KimUnpublished" = "Kim et al. (unpublished)",
  "LaBarUnpublished" = "LaBar et al. (unpublished)",
  "MarinMorales et al. (2021)" = "Marín-Morales et al. (2022)",
  "MinUnpublished" = "Min et al. (2022)",
  "Morawetz et al. (2016a)" = "Morawetz et al. (2016)",
  "Morawetz et al. (2016b)" = "Morawetz et al. (2016; pictures)",
  "Morawetz et al. (2019)" = "Morawetz et al. (2016; videos)",
  "Morawetz et al. (2020)" = "Morawetz et al. (2019)",
  "Morawetz et al. (2021)" = "Morawetz et al. (2020)",
  "MulejBratec et al. (2015)" = "Morawetz et al. (2021)",
  "MuellerPinzlerUnpublished" = "Mulej Bratec et al. (2015)",
  "HunekeUnpublished" = "Müller-Pinzler et al. (unpublished)",
  "NetaUnpublished" = "Huneke et al. (unpublished)",
  "Paschke et al. (2016)" = "Pierce et al. (2022)",
  "Rehbein et al. (2021a)" = "Rehbein et al. (2021; sample 1)",
  "Rehbein et al. (2021b)" = "Rehbein et al. (2021; sample 2)",
  "Sandner et al. (2021)" = "Sandner et al. (2021)",
  "SidorenkoUnpublished" = "Doren et al. (unpublished)",
  "Guendelman&DziobekUnpublished" = "Guendelman et al. (2022)",
  "SokolowskiUnpublished" = "Sokolowski et al. (2022)",
  "Steward et al. (2021)" = "Steward et al. (2021)",
  "VanReekum et al. (2021)" = "Lloyd et al. (2021)",
  "VanReekumUnpublished" = "Tupitsa et al. (2023)"
)

dfcorrs <- dfcorrs %>%
  mutate(
    cleaned_author = ifelse(
      `formatted_author` %in% names(author_mapping),
      author_mapping[formatted_author],
      formatted_author
    )
  )

dfcorrs$`image name stem` <- dfcorrs$cleaned_author
onlineRatings <- !dfcorrs$`image name stem` %in% c("Diers, Gärtner et al. (2023)", "Scheffel et al. (2019)" , "Hofhansel et al. (2023)")
#onlineRatings <- !dfcorrs$`image name stem` %in% c("Diers2021", "GaertnerDiersUnpublished", "HofhanselUnpublished")


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


# Bayes factors
# tau prior based on https://openpsychologydata.metajnl.com/articles/10.5334/jopd.33
# r prior based on: Effect size guidelines for individual differences researchers &
# Efficient alternatives for Bayesian hypothesis tests in psychology

dfQuestSelf$questSelfRating <- as.numeric(dfQuestSelf$questSelfRating) 
dfQuestSelf$'sample size quest' <- as.numeric(dfQuestSelf$'sample size quest')

dfQuestSelf_conv <- escalc(measure="ZCOR", ri=questSelfRating, ni=dfQuestSelf$'sample size quest', data=dfQuestSelf)
dfQuestSelf_conv$vi_sqrt <- sqrt(dfQuestSelf_conv$vi)

tau2 <- 0.2^2/2
effSize_prior_oneSide <- prior("custom", function(x) x^2/tau2 * dnorm(x, sd=sqrt(tau2)), 0, 1)
effSize_prior_twoSide <- prior("custom", function(x) x^2/tau2 * dnorm(x, sd=sqrt(tau2)), -1, 1)
tau_prior <- prior("t", c(location=0, scale=.15, nu=1), lower=0, upper=1)


Bmeta_questSelf <- meta_bma(y=yi, SE=vi_sqrt, data=dfQuestSelf_conv,
                            d = effSize_prior_twoSide,
                            tau = tau_prior,
                            iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")
Bmeta_questSelf



plot(effSize_prior_twoSide, xaxt = "n")
axis(1, at = seq(-1, 1, by = 0.1))

plot(tau_prior, xaxt = "n")
axis(1, at = seq(0, 1, by = 0.1))


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

# Bayes factors
dfAmyQuest$amyQuestCorr <- as.numeric(dfAmyQuest$amyQuestCorr) 

dfAmyQuest_conv <- escalc(measure="ZCOR", ri=amyQuestCorr, ni=dfAmyQuest$'sample size quest', data=dfAmyQuest)
dfAmyQuest_conv$vi_sqrt <- sqrt(dfAmyQuest_conv$vi)

Bmeta_amyQuest <- meta_bma(y=yi, SE=vi_sqrt, data=dfAmyQuest_conv,
                            d = effSize_prior_twoSide,
                            tau = tau_prior,
                            iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")

Bmeta_amyQuest


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


# Bayes factors
dfAmySelf$amySelfRating <- as.numeric(dfAmySelf$amySelfRating) 

dfAmySelf_conv <- escalc(measure="ZCOR", ri=amySelfRating, ni=dfAmySelf$'sample size quest', data=dfAmySelf)
dfAmySelf_conv$vi_sqrt <- sqrt(dfAmySelf_conv$vi)

Bmeta_amySelf <- meta_bma(y=yi, SE=vi_sqrt, data=dfAmySelf_conv,
                           d = effSize_prior_twoSide,
                           tau = tau_prior,
                           iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")

Bmeta_amySelf



#######################################
# plot power curve

metaPower(effectSize = 0.2, sampleSizes = as.numeric(meta_questSelf$n), tau = meta_amyQuest$tau)


# grid of effect sizes
eff <- seq(0, 1, length.out = 201)

# helper to build one curve; coerces n to numeric just in case
curve_df <- function(model, label) {
  n_vec <- as.numeric(model$n)
  pow <- vapply(eff, function(e) {
    res <- metaPower(effectSize = e, sampleSizes = n_vec, tau = model$tau)
    # Prefer $randomPower; fall back to $random or first numeric if needed
    if (!is.null(res$randomPower)) {
      as.numeric(res$randomPower)
    } else if (!is.null(res$random)) {
      as.numeric(res$random)
    } else {
      as.numeric(res[[1]])
    }
  }, numeric(1))
  data.frame(effect = eff, power = pow, model = label)
}



# build curves
df <- rbind(
  curve_df(meta_questSelf, "questSelf"),
  curve_df(meta_amyQuest,  "amyQuest"),
  curve_df(meta_amySelf,   "amySelf")
)

# plot
ggplot(df, aes(effect, power, color = model)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.80, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0.10, color = "red", linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_color_discrete(
    breaks = c("amyQuest", "amySelf", "questSelf"),
    labels = c(
      "amyQuest"  = "Amygdala/Questionnaires",
      "amySelf"   = "Amygdala/Ratings",
      "questSelf" = "Questionnaires/Ratings"
    )
  ) +
  labs(x = "Effect size", y = "Power", color = "Model",
       title = "Meta-analytic power curves") +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.02),   # inside, bottom-right
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "grey80")
  )

ggsave(here("results", "figures", "overlap_powerCurved.svg"), width = 7, height = 4.5, units = "in", dpi = 300)




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

defineOverlap = c(
  "Trait\nquestionnaires" = 100,
  "Task-based\naffective ratings" = 100,
  "Amygdala\ndown-regulation" = 100,
  "Trait\nquestionnaires&Task-based\naffective ratings" = round(DG*100,2),
  "Trait\nquestionnaires&Amygdala\ndown-regulation" = round(EG*100,3),
  "Task-based\naffective ratings&Amygdala\ndown-regulation" = round(FG*100,2)
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

