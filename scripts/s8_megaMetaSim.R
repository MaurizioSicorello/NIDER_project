##################
# load packages and sample sizes

library("here")
library("metafor")
library("ggplot2")
library("psych")
library("tidyverse")
library(dplyr)
library(lme4)
library(meta)

Nemp <- read.csv(here("results", "bayesfactors", "NOI_quest.csv"))$sampleSize


#################################### FUNCTIONS ####################################

##################
# function to simulate data
simulate_mega_cor <- function(n_vec, rho, tau, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  k <- length(n_vec)
  
  # draw study-specific correlations
  rho_i <- rnorm(k, mean = rho, sd = tau)
  
  # truncate to valid correlation range
  rho_i[rho_i >= 0.999] <- 0.999
  rho_i[rho_i <= -0.999] <- -0.999
  
  dat_list <- vector("list", k)
  
  for (i in seq_len(k)) {
    
    n <- n_vec[i]
    r <- rho_i[i]
    
    Sigma <- matrix(c(1, r,
                      r, 1), 2, 2)
    
    xy <- MASS::mvrnorm(n = n, mu = c(0,0), Sigma = Sigma)
    
    dat_list[[i]] <- data.frame(
      study = i,
      x = xy[,1],
      y = xy[,2]
    )
  }
  
  dat <- do.call(rbind, dat_list)
  
  list(
    data = dat,
    true_r = rho_i,
    n = n_vec,
    k = k
  )
}


##################
# function to randomly remove studies
drop_studies <- function(sim_obj, drop_prop = 0, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if (drop_prop <= 0) return(sim_obj)
  if (drop_prop >= 1) stop("drop_prop must be smaller than 1.")
  
  k <- sim_obj$k
  n_drop <- floor(k * drop_prop)
  
  # keep at least 2 studies
  if ((k - n_drop) < 2) {
    n_drop <- k - 2
  }
  
  if (n_drop <= 0) return(sim_obj)
  
  dropped_studies <- sample(seq_len(k), size = n_drop, replace = FALSE)
  kept_studies <- setdiff(seq_len(k), dropped_studies)
  
  dat_new <- sim_obj$data %>%
    filter(study %in% kept_studies)
  
  # renumber studies consecutively
  study_map <- data.frame(
    old_study = kept_studies,
    study = seq_along(kept_studies)
  )
  
  dat_new <- dat_new %>%
    left_join(study_map, by = c("study" = "old_study")) %>%
    select(study, x, y)
  
  list(
    data = dat_new,
    true_r = sim_obj$true_r[kept_studies],
    n = sim_obj$n[kept_studies],
    k = length(kept_studies),
    dropped_studies = dropped_studies,
    kept_studies = kept_studies
  )
}


##################
# function to analyze data
analyze_mega_meta <- function(dat) {
  
  # 1) z-standardize within studies
  dat_std <- dat %>%
    group_by(study) %>%
    mutate(
      x = as.numeric(scale(x)),
      y = as.numeric(scale(y))
    ) %>%
    ungroup()
  
  # 2) mega-analysis: no fixed or random intercepts
  mega_mod <- lmer(
    y ~ 0 + x + (0 + x | study),
    data = dat_std,
    REML = TRUE
  )
  
  mega_est <- unname(fixef(mega_mod)["x"])
  mega_se  <- sqrt(vcov(mega_mod)["x", "x"])
  mega_ci  <- mega_est + c(-1, 1) * qnorm(0.975) * mega_se
  
  mega_res <- data.frame(
    estimate = mega_est,
    ci_lb = mega_ci[1],
    ci_ub = mega_ci[2]
  )
  
  # 3) study-level correlations for meta-analysis
  cor_dat <- dat_std %>%
    group_by(study) %>%
    summarise(
      r = cor(x, y),
      n = n(),
      .groups = "drop"
    )
  
  # 4) random-effects meta-analysis via metacor
  meta_mod <- metacor(
    cor = r,
    n = n,
    studlab = study,
    data = cor_dat,
    method.tau = "REML",
    method.random.ci = "classic",
    sm = "ZCOR",
    backtransf = TRUE,
    control = list(maxiter = 1000, stepadj = 0.5)
  )
  
  meta_res <- data.frame(
    estimate = tanh(meta_mod$TE.random),
    ci_lb = tanh(meta_mod$lower.random),
    ci_ub = tanh(meta_mod$upper.random)
  )
  
  list(
    mega = mega_res,
    meta = meta_res,
    mega_model = mega_mod,
    meta_model = meta_mod,
    study_effects = cor_dat
  )
}


#################################### SIMULATION ####################################

t1 <- Sys.time()

iterations <- 2000
drop_prop <- 0.80   

compOut <- matrix(NA, nrow = iterations, ncol = 6)

for(i in 1:iterations){
  
  simTemp <- simulate_mega_cor(
    n_vec = Nemp,
    rho = 0.20,
    tau = 0.10
  )
  
  simTemp <- drop_studies(
    sim_obj = simTemp,
    drop_prop = drop_prop
  )
  
  resTemp <- analyze_mega_meta(simTemp$data)
  
  compOut[i, 1:3] <- as.numeric(resTemp$mega)
  compOut[i, 4:6] <- as.numeric(resTemp$meta)
}

t2 <- Sys.time()
t2-t1

compOut <- as.data.frame(compOut)
names(compOut) <- c("mega_mean", "mega_CI_l", "mega_CI_u", "meta_mean", "meta_CI_l", "meta_CI_u")

write.csv(compOut, here::here("results", "power", paste0("megaMetaSim_drop", drop_prop)))



############################
# evaluate results

df0 <- read.csv(here::here("results", "power", "megaMetaSim_drop0"))[,-1]
#df0.05 <- read.csv(here::here("results", "power", "megaMetaSim_drop0.05"))[,-1]
df0.1 <- read.csv(here::here("results", "power", "megaMetaSim_drop0.1"))[,-1]
df0.2 <- read.csv(here::here("results", "power", "megaMetaSim_drop0.2"))[,-1]

dfplot <- as.data.frame(rbind(as.matrix(df0[,4:6]), as.matrix(df0[,1:3]), as.matrix(df0.1[,1:3]), as.matrix(df0.2[,1:3])))
names(dfplot) <- c("mean", "CI_l", "CI_u")
modelNames <- c("Meta", "Mega", "Mega - 10% missing", "Mega - 20% missing")
dfplot$model <- factor(rep(modelNames, each = iterations), levels = modelNames)


ggplot(dfplot, aes(x = model, y = mean)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.12
  ) +
  theme_classic() +
  labs(
    x = "Model",
    y = "Mean"
  ) +
  geom_hline(yintercept = 0.20, linetype=2) +
  ylim(-0.1, 0.5)




dfplot_prec <- dfplot %>%
  mutate(
    precision = CI_u - CI_l
  )

describeBy(dfplot_prec, group="model")


ggplot(dfplot_prec, aes(x = model, y = precision)) +
  geom_violin(trim = FALSE, alpha = 0.3) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.12
  ) +
  theme_classic() +
  labs(
    x = "Model",
    y = "CI width"
  ) +
  ylim(0, 0.3)






# adapt these column selections if needed:
# columns 1:3 = mega, columns 4:6 = meta
dfplot <- bind_rows(
  data.frame(df0[, 4:6], model = "Meta"),
  data.frame(df0[, 1:3], model = "Mega")
)

# rename columns if necessary
names(dfplot)[1:3] <- c("mean", "CI_l", "CI_u")

# create precision and stack into long format
dfplot_long <- dfplot %>%
  mutate(precision = CI_u - CI_l) %>%
  select(model, mean, precision) %>%
  pivot_longer(
    cols = c(mean, precision),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    mean = "Mean estimate",
                    precision = "Precision (CI width)")
  )



df_meta <- as.data.frame(df0[, 4:6])
names(df_meta) <- c("mean", "CI_l", "CI_u")
df_meta$model <- "Meta"

df_mega <- as.data.frame(df0[, 1:3])
names(df_mega) <- c("mean", "CI_l", "CI_u")
df_mega$model <- "Mega"

dfplot <- bind_rows(df_meta, df_mega)

dfplot_long <- dfplot %>%
  mutate(precision = CI_u - CI_l) %>%
  select(model, mean, precision) %>%
  pivot_longer(
    cols = c(mean, precision),
    names_to = "metric",
    values_to = "value"
  )

dfplot_long <- dfplot %>%
  mutate(precision = CI_u - CI_l) %>%
  select(model, mean, precision) %>%
  pivot_longer(
    cols = c(mean, precision),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    mean = "Mean estimate",
                    precision = "Precision (CI width)")
  )


ggplot(dfplot_long, aes(x = model, y = value, fill = model, colour = model)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.12
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Model",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )




library(dplyr)
library(ggplot2)
library(patchwork)

df0 <- read.csv(here::here("results", "power", "megaMetaSim_drop0"))[, -1]

df_meta <- as.data.frame(df0[, 4:6])
names(df_meta) <- c("mean", "CI_l", "CI_u")
df_meta$model <- "Meta-analysis"

df_mega <- as.data.frame(df0[, 1:3])
names(df_mega) <- c("mean", "CI_l", "CI_u")
df_mega$model <- "Mega-analysis"

dfplot <- bind_rows(df_meta, df_mega) %>%
  mutate(
    ci_width = CI_u - CI_l,
    model = factor(model, levels = c("Meta-analysis", "Mega-analysis"))
  )

p_mean <- ggplot(dfplot, aes(x = model, y = mean, fill = model, colour = model)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.12
  ) +
  geom_hline(yintercept = 0.20, linetype = 2) +
  scale_y_continuous(
    limits = c(-0.10, 0.50),
    breaks = seq(-0.10, 0.50, by = 0.10)
  ) +
  theme_classic() +
  labs(
    title = "Mean estimate",
    x = "Model",
    y = "Mean estimate"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_prec <- ggplot(dfplot, aes(x = model, y = ci_width, fill = model, colour = model)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.12
  ) +
  scale_y_continuous(
    limits = c(0, 0.30),
    breaks = seq(0, 0.30, by = 0.05)
  ) +
  theme_classic() +
  labs(
    title = "Precision",
    x = "Model",
    y = "CI width"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_mean + p_prec

ggsave(here::here("results", "figures", "mega_meta_two_column.png"), plot = p_mean + p_prec, width = 7, height = 3.5, units = "in", dpi = 300)



mean(df0[,1] - df0[,4])
mean(df0[,4]-0.20)
mean(df0[,1]-0.20)
sd(df0[,1] - df0[,4])

mean((df0[,3] - df0[,2]) - (df0[,6] - df0[,5]))
sd((df0[,3] - df0[,2]) - (df0[,6] - df0[,5]))


#####################################################
# simulate power by sample size

source(here("functions", "metaPoweR.R"))

# parameters
effVal <- 0.10
tauVec <- c(0, 0.15, 0.30)
numTruePos <- 45
voxelFDR <- 10000
perStudyN <- 57

# vary number of studies, starting at 5
kVec <- 5:37

# create parameter grid
outTable <- expand.grid(
  tau = tauVec,
  k = kVec
)

# sample size vectors and total N
outTable$sampleSizes <- lapply(outTable$k, function(k) rep(perStudyN, k))
outTable$totalN <- outTable$k * perStudyN

# compute power
outTable$power <- mapply(
  function(tau, sampleSizes) {
    metaPower(
      effectSize = effVal,
      sampleSizes = sampleSizes,
      tau = tau,
      alpha = 0.05 / voxelFDR,
      numTruePositives = numTruePos
    )$random
  },
  tau = outTable$tau,
  sampleSizes = outTable$sampleSizes
)

# labels
outTable$tau <- factor(
  outTable$tau,
  levels = c(0, 0.15, 0.30),
  labels = c("tau = 0", "tau = 0.15", "tau = 0.30")
)

# axis breaks
kBreaks <- c(5, 10, 15, 20, 25, 30, 35)
nBreaks <- kBreaks * perStudyN

# plot
ggplot(outTable, aes(x = k, y = power, color = tau)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  scale_x_continuous(
    breaks = kBreaks,
    labels = kBreaks,
    name = "Number of studies",
    sec.axis = sec_axis(
      ~ . * perStudyN,
      breaks = nBreaks,
      labels = nBreaks,
      name = "Total sample size"
    )
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    y = "Power",
    color = "Between-study heterogeneity"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x.top = element_text(color = "black"),
    axis.title.x.top = element_text(color = "black"),
    axis.text.x.bottom = element_text(color = "black"),
    axis.title.x.bottom = element_text(color = "black")
  )

# save if wanted
ggsave(here("results", "figures", "power_by_k_tau_lines.png"), width = 8, height = 6)
