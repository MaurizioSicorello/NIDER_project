

library("metaBMA")
library("here")
library("metafor")

library("future")
library("future.apply")
library("progressr")


########################
# define priors
tau2 <- 0.2^2/2
effSize_prior_oneSide <- prior("custom", function(x) x^2/tau2 * dnorm(x, sd=sqrt(tau2)), 0, 1)
effSize_prior_twoSide <- prior("custom", function(x) x^2/tau2 * dnorm(x, sd=sqrt(tau2)), -1, 1)
tau_prior <- prior("t", c(location=0, scale=.15, nu=1), lower=0, upper=1)

# function to convert t-values to r

convert_t <- function(t, N){
  r <- t/sqrt(t^2 + N-2)
  return(r)
}




########################
# average network activity

############
# questionnaire

df_NOI_quest <- read.csv(here("results", "bayesfactors", "NOI_quest.csv"))

# Network 1
df_NOI_quest$networkOne <- convert_t(df_NOI_quest$networkOne, df_NOI_quest$sampleSize)
df_NOI_quest_conv_n1 <- escalc(measure="ZCOR", ri=networkOne, ni=sampleSize, data=df_NOI_quest)
df_NOI_quest_conv_n1$vi_sqrt <- sqrt(df_NOI_quest_conv_n1$vi)

Bmeta_NOI_quest <- meta_bma(y=yi, SE=vi_sqrt, data=df_NOI_quest_conv_n1,
                            d = effSize_prior_twoSide,
                            tau = tau_prior,
                            iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")
Bmeta_NOI_quest

# Network 2
df_NOI_quest$networkTwo <- convert_t(df_NOI_quest$networkTwo, df_NOI_quest$sampleSize)
df_NOI_quest_conv_n1 <- escalc(measure="ZCOR", ri=networkTwo, ni=sampleSize, data=df_NOI_quest)
df_NOI_quest_conv_n1$vi_sqrt <- sqrt(df_NOI_quest_conv_n1$vi)

Bmeta_NOI_quest <- meta_bma(y=yi, SE=vi_sqrt, data=df_NOI_quest_conv_n1,
                            d = effSize_prior_twoSide,
                            tau = tau_prior,
                            iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")
Bmeta_NOI_quest


# ability questionnaires for revision

df_NOI_quest <- read.csv(here("results", "bayesfactors", "NOI_quest_ability.csv"))

# Network 1
df_NOI_quest$networkOne <- convert_t(df_NOI_quest$networkOne, df_NOI_quest$sampleSize)
df_NOI_quest_conv_n1 <- escalc(measure="ZCOR", ri=networkOne, ni=sampleSize, data=df_NOI_quest)
df_NOI_quest_conv_n1$vi_sqrt <- sqrt(df_NOI_quest_conv_n1$vi)

Bmeta_NOI_quest <- meta_bma(y=yi, SE=vi_sqrt, data=df_NOI_quest_conv_n1,
                            d = effSize_prior_twoSide,
                            tau = tau_prior,
                            iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")
Bmeta_NOI_quest

# Network 2
df_NOI_quest$networkTwo <- convert_t(df_NOI_quest$networkTwo, df_NOI_quest$sampleSize)
df_NOI_quest_conv_n1 <- escalc(measure="ZCOR", ri=networkTwo, ni=sampleSize, data=df_NOI_quest)
df_NOI_quest_conv_n1$vi_sqrt <- sqrt(df_NOI_quest_conv_n1$vi)

Bmeta_NOI_quest <- meta_bma(y=yi, SE=vi_sqrt, data=df_NOI_quest_conv_n1,
                            d = effSize_prior_twoSide,
                            tau = tau_prior,
                            iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")
Bmeta_NOI_quest




############
# rating

df_NOI_rating <- read.csv(here("results", "bayesfactors", "NOI_rating.csv"))

# Network 1
df_NOI_rating$networkOne <- convert_t(df_NOI_rating$networkOne, df_NOI_rating$sampleSize)
df_NOI_rating_conv_n1 <- escalc(measure="ZCOR", ri=networkOne, ni=sampleSize, data=df_NOI_rating)
df_NOI_rating_conv_n1$vi_sqrt <- sqrt(df_NOI_rating_conv_n1$vi)

Bmeta_NOI_rating <- meta_bma(y=yi, SE=vi_sqrt, data=df_NOI_rating_conv_n1,
                            d = effSize_prior_twoSide,
                            tau = tau_prior,
                            iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")
Bmeta_NOI_rating

# Network 2
df_NOI_rating$networkTwo <- convert_t(df_NOI_rating$networkTwo, df_NOI_rating$sampleSize)
df_NOI_rating_conv_n1 <- escalc(measure="ZCOR", ri=networkTwo, ni=sampleSize, data=df_NOI_rating)
df_NOI_rating_conv_n1$vi_sqrt <- sqrt(df_NOI_rating_conv_n1$vi)

Bmeta_NOI_rating <- meta_bma(y=yi, SE=vi_sqrt, data=df_NOI_rating_conv_n1,
                            d = effSize_prior_twoSide,
                            tau = tau_prior,
                            iter = 8000, logml_iter = 5000, rel.tol = .1, summarize="integrate")
Bmeta_NOI_rating



########################
# voxel-wise analyses

# --- Wrapper function for one voxel ---
run_meta <- function(voxel_data, sample_size) {
  # Drop studies with missing values for this voxel
  keep <- !is.na(voxel_data)
  voxel_data <- voxel_data[keep]
  sample_size <- sample_size[keep]
  
  # If too few studies remain, skip
  if (length(voxel_data) < 3) return(NA_real_)
  
  # Preprocessing
  voxel_conv <- convert_t(voxel_data, sample_size)
  voxel_es <- escalc(measure = "ZCOR", ri = voxel_conv, ni = sample_size)
  voxel_es$vi_sqrt <- sqrt(voxel_es$vi)
  
  # Meta-analysis
  Bmeta <- meta_bma(
    y   = voxel_es$yi,
    SE  = voxel_es$vi_sqrt,
    #data = voxel_es,
    d   = effSize_prior_twoSide,
    tau = tau_prior,
    iter = 8000, logml_iter = 5000, rel.tol = .1,
    summarize = "integrate"
  )
  
  # Return just the inclusion Bayes Factor
  Bmeta$inclusion$incl.BF
}

safe_run_meta <- function(voxel_data, sample_size) {
  tryCatch(
    run_meta(voxel_data, sample_size),
    error = function(e) NA_real_   # return NA instead of stopping
  )
}


###########
# questionnaires

# prepare df
#df <- read.csv(here("results", "bayesfactors", "WB_quest.csv"))
#df <- read.csv(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/WB_quest.csv"))
# df <- read.csv(here("results", "bayesfactors", "WB_quest_ability.csv"))
# 
# sampleSize <- df$sampleSize
# datTable <- as.matrix(df[, -(1:2)])
# 
# #datTable <- datTable[,1:14]
# 
# plan(multisession, workers = availableCores() - 1)
# 
# t1 <- Sys.time()
# 
# BF_vector_quest <- with_progress({
#   
#   # define a progressor inside with_progress()
#   p <- progressor(along = seq_len(ncol(datTable)))
#   
#   future_sapply(
#     seq_len(ncol(datTable)),
#     function(j) {
#       p()  # increment progress bar
#       safe_run_meta(datTable[, j], sampleSize)
#     },
#     future.seed = TRUE
#   )
# })
# 
# 
# t2 <- Sys.time()
# t2-t1
# 
# plan(sequential)
# 
# #write.csv(BF_vector_quest, here("results", "bayesfactors", "WB_quest_results.csv"), row.names = FALSE)
# #write.csv(BF_vector_quest, paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/WB_quest_results.csv"), row.names = FALSE)
# write.csv(BF_vector_quest, here("results", "bayesfactors", "WB_quest_results_ability.csv"), row.names = FALSE)
# 
# psych::describe(BF_vector_quest)
# 
# sum(is.na(datTable[1,]))

###########
# ratings

# # prepare df
# #df <- read.csv(here("results", "bayesfactors", "WB_quest.csv"))
# df <- read.csv(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/WB_rating.csv"))
# 
# sampleSize <- df$sampleSize
# datTable <- as.matrix(df[, -(1:2)])
# 
# #datTable <- datTable[,1:14]
# 
# plan(multisession, workers = availableCores() - 1)
# 
# t1 <- Sys.time()
# 
# BF_vector_quest <- with_progress({
#   
#   # define a progressor inside with_progress()
#   p <- progressor(along = seq_len(ncol(datTable)))
#   
#   future_sapply(
#     seq_len(ncol(datTable)),
#     function(j) {
#       p()  # increment progress bar
#       safe_run_meta(datTable[, j], sampleSize)
#     },
#     future.seed = TRUE
#   )
# })
# 
# 
# t2 <- Sys.time()
# t2-t1
# 
# plan(sequential)
# 
# #write.csv(BF_vector_quest, here("results", "bayesfactors", "WB_quest_results.csv"), row.names = FALSE)
# write.csv(BF_vector_quest, paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/WB_rating_results.csv"), row.names = FALSE)


########################################################
# evaluate whole-brain BFs

bf_quests <- read.csv(here("results", "bayesfactors", "WB_quest_results.csv"))$x
bf_ratings <- read.csv(here("results", "bayesfactors", "WB_rating_results.csv"))$x
bf_quests2 <- read.csv(here("results", "bayesfactors", "WB_quest_results_ability.csv"))$x


plot(bf_quests, bf_ratings)

# prior prob based on effective number of tests (=53, derived in s4_RpowerSim.R)
H0_prior_prob <- 0.5^(1/53)
H1_prior_prob <- 1-H0_prior_prob
prior_odds <- H1_prior_prob/H0_prior_prob

# inspect questionnaire BFs
sum(bf_quests < 0.33, na.rm =TRUE)
sum(bf_quests > 3, na.rm =TRUE)

bf_quests2 <- bf_quests2[!is.na(bf_quests2)]
table(ifelse(bf_quests2>3, "h1", ifelse(bf_quests2 < 1/3, "h0", "inconcl.")))
table(ifelse((bf_quests2*prior_odds)>3, "h1", ifelse((bf_quests2*prior_odds) < 1/3, "h0", "inconcl.")))
psych::describe(bf_quests2*prior_odds)
hist(bf_quests2[bf_quests2>3])

max(bf_quests2*prior_odds)

bf_quests_c <- bf_quests * prior_odds
sum(bf_quests_c < 0.33, na.rm =TRUE)/length(bf_quests_c)
sum(bf_quests_c < 0.1, na.rm =TRUE)/length(bf_quests_c)
sum(bf_quests_c < 0.01, na.rm =TRUE)/length(bf_quests_c)
sum(bf_quests_c > 0.33, na.rm =TRUE)

# inspect ratings BFs
sum(bf_ratings < 0.33, na.rm =TRUE)
sum(bf_ratings > 3, na.rm =TRUE)

hist(bf_ratings[bf_ratings > 3 & !is.na(bf_ratings)])

bf_ratings_c <- bf_ratings * prior_odds
sum(bf_ratings_c < 0.01, na.rm =TRUE)/length(bf_ratings_c)
sum(bf_ratings_c > 3, na.rm =TRUE)
psych::describe(bf_ratings_c[bf_ratings_c < 0.33], na.rm =TRUE)

sum(bf_ratings_c < 3 & bf_ratings_c > 1/3, na.rm =TRUE)
psych::describe(bf_ratings_c[bf_ratings_c < 3 & bf_ratings_c > 1/3], na.rm =TRUE)

bf_ratings_c[bf_ratings_c > 3 & !is.na(bf_ratings)]
psych::describe(bf_ratings_c[bf_ratings_c > 3], na.rm =TRUE)

sum(bf_ratings > 3 & !is.na(bf_ratings))
psych::describe(bf_ratings[bf_ratings > 3], na.rm =TRUE)
