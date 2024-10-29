library("ggplot2")
library("here")
library("cowplot")



##########################################################################
# PLOT EFFECT SIZE DENSITY DISTRIBUTION

# load an combine data
effSizes_quest <- read.csv(here::here("results", "tables", "quest_effectSizes.csv"))
effSizes_quest$tau_quest <- rep(0, nrow(effSizes_quest))
effSizes_rating <- read.csv(here::here("results", "tables", "rating_effectSizes_withoutBrehl.csv"))
effSizes_rating$sig_rating_withoutBrehl <- factor(effSizes_rating$sig_rating_withoutBrehl, levels=c("", "significant"))
effSizes_amy <- read.csv(here::here("results", "tables", "amygdala_effectSizes.csv"))

questPlot <- ggplot(data=effSizes_quest, aes(x=effectSizes_quest)) +
  geom_histogram(binwidth=0.001) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, by = 0.15)) +
  theme_classic() +
  
  xlab("Questionnaires") +
  ylab(NULL)


ratingPlot <- ggplot(data=effSizes_rating, aes(x=effectSizes_rating_withoutBrehl)) +
  geom_histogram(binwidth=0.001) +
  #geom_density() +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, by = 0.15)) +
  theme_classic() +
  
  xlab("Ratings") +
  ylab(NULL)


amyPlot <- ggplot(data=effSizes_amy, aes(x=effectSizes_quest)) +
  geom_histogram(binwidth=0.001) +
  scale_x_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, by = 0.15)) +
  theme_classic() +
  
  xlab("Amygdala") +
  ylab(NULL)


plot_grid(questPlot, ratingPlot, amyPlot, nrow=1)




##########################################################################
# PLOT BRAIN-BRAIN CORRELATIONS ACROSS MODALITIES

dfcrossMod <- read.csv(here::here("results", "tables", "kragelAnalyses.csv"))

mean(dfcrossMod$corrMean)

# Custom labels for the x-axis (Subdomain) and the facets (Domain)
custom_x_labels <- c("Inhibition" = "Inhibitory Control", 
                     "ResponseSelect" = "Response Selection", 
                     "WorkingMem" = "Working Memory", 
                     "Images" = "Aversive Images", 
                     "Social" = "Social Rejection", 
                     "Sounds" = "Aversive Sounds", 
                     "Mechanical" = "Mechanical Pain", 
                     "Thermal" = "Thermal Pain", 
                     "Visceral" = "Visceral Pain")

custom_facet_labels <- c("Cog_control" = "Cognitive Control", 
                         "Neg_Emotion" = "Negative Emotion", 
                         "Pain" = "Pain Processing")

ggplot(dfcrossMod, aes(x = Subdomain, y = corrMean)) +
  
  geom_hline(yintercept=0, colour="lightgrey") +
  
  geom_point(aes(group = Studynumber), position = position_dodge(width = 0.7), size = 3) +  # Plot means as dots
  geom_errorbar(aes(ymin = corrMean - corrSD, ymax = corrMean + corrSD, group = Studynumber), 
                width = 0.2, position = position_dodge(width = 0.7)) +  # Add error bars for SD
  
  scale_x_discrete(labels = custom_x_labels) +
  facet_wrap(~Domain, scales = "free_x", nrow = 1, strip.position = "bottom", labeller = as_labeller(custom_facet_labels)) +  # Move domain indicators to bottom
  theme_classic() +
  theme(strip.placement = "outside",  # Move facet labels outside the plot
        strip.background = element_blank(),  # Remove the box around facet labels
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(x = NULL, y = "Average inter-region correlation") 


ggsave(here::here("results", "figures", "kragelAnalyses.svg"), device="svg", height = 6.75/2, width=6.75)

