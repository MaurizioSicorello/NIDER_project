library("ggplot2")
library("here")
library("cowplot")

# load an combine data
effSizes_quest <- read.csv(here::here("results", "tables", "quest_effectSizes.csv"))
effSizes_quest$tau_quest <- rep(0, nrow(effSizes_quest))
effSizes_rating <- read.csv(here::here("results", "tables", "rating_effectSizes_withoutBrehl.csv"))
effSizes_rating$sig_rating_withoutBrehl <- factor(effSizes_rating$sig_rating_withoutBrehl, levels=c("0", "1"))
effSizes_amy <- read.csv(here::here("results", "tables", "amygdala_effectSizes.csv"))

questPlot <- ggplot(data=effSizes_quest, aes(x=effectSizes_quest)) +
  geom_histogram(binwidth=0.001) +
  xlim(-0.2, 0.2) +
  theme_classic() +
  
  xlab("Correlations") +
  ylab("Density")


ratingPlot <- ggplot(data=effSizes_rating, aes(x=effectSizes_rating_withoutBrehl)) +
  geom_histogram(binwidth=0.001) +
  #geom_density(aes(y = after_stat(density * n/nrow(effSizes_rating)))) +
  xlim(-0.2, 0.2) +
  theme_classic() +
  
  xlab("Correlations") +
  ylab("Density")


amyPlot <- ggplot(data=effSizes_amy, aes(x=effectSizes_quest)) +
  geom_histogram(binwidth=0.001) +
  xlim(-0.8, 0.8) +
  theme_classic() +
  
  xlab("Correlations") +
  ylab("Density")



cowplot::plot_grid(questPlot, ratingPlot, amyPlot, nrow=1)




# plot jackknife rating results
jackknife_rating <- read.csv(here::here("results", "tables", "jackknife_ratings.csv"))

ggplot(data=jackknife_rating, aes(x=numSigEffects)) +
  geom_histogram(colour="black") +
  theme_classic() +
  ylab("Absolute Study Frequency") +
  xlab("# Significant Voxels after Study Inclusion")


ggsave(here::here("results", "figures", "jackknife_ratings.png"))
ggsave(here::here("results", "figures", "jackknife_ratings.svg"))



# plot amygdala effect sizes
effSize_amy <- read.csv(here::here("results", "tables", "amygdala_effectSizes.csv"))

ggplot(data=effSize_amy, aes(x=effectSizes_quest)) +
  geom_histogram(colour="black", binwidth=0.01) +
  theme_classic() +
  ylab("Absolute Frequency") +
  xlab("Correlations of Amygdala and remaining brain (between-person)") +
  xlim(c(-1, 1))


ggsave(here::here("results", "figures", "amygdala_effectSizes.png"))
ggsave(here::here("results", "figures", "amygdala_effectSizes.svg"))

round(mean(effSize_amy$effectSizes_quest), 2)
round(range(effSize_amy$effectSizes_quest), 2)
