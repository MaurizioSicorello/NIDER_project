library(tidyverse)
library(here)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

dat <- read.csv(
  here("results", "tables", "kragelAnalyses.csv")
)


# ------------------------------------------------------------
# Prepare labels
# ------------------------------------------------------------

dat <- dat %>%
  mutate(
    domain_label = recode(
      Domain,
      "Cog_control" = "Cognitive Control",
      "Neg_Emotion" = "Negative Emotion",
      "Pain"        = "Pain Processing"
    ),
    
    subdomain_label = recode(
      Subdomain,
      "Inhibition"     = "Inhibitory Control",
      "ResponseSelect" = "Response Selection",
      "WorkingMem"     = "Working Memory",
      "Images"         = "Aversive Images",
      "Social"         = "Social Rejection",
      "Sounds"         = "Aversive Sounds",
      "Mechanical"     = "Mechanical Pain",
      "Thermal"        = "Thermal Pain",
      "Visceral"       = "Visceral Pain"
    )
  )


# ------------------------------------------------------------
# Set order
# ------------------------------------------------------------

dat$domain_label <- factor(
  dat$domain_label,
  levels = c(
    "Cognitive Control",
    "Negative Emotion",
    "Pain Processing"
  )
)

subdomain_order <- c(
  "Inhibitory Control",
  "Response Selection",
  "Working Memory",
  "Aversive Images",
  "Social Rejection",
  "Aversive Sounds",
  "Mechanical Pain",
  "Thermal Pain",
  "Visceral Pain"
)

dat$subdomain_label <- factor(
  dat$subdomain_label,
  levels = subdomain_order
)


# ------------------------------------------------------------
# Create ID for the two observations within each subdomain
# ------------------------------------------------------------

dat <- dat %>%
  group_by(domain_label, subdomain_label) %>%
  mutate(point_id = factor(row_number())) %>%
  ungroup()


# ------------------------------------------------------------
# Horizontal separation of observations
# ------------------------------------------------------------

pd <- position_dodge(width = 0.30)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(
  dat,
  aes(
    x = subdomain_label,
    y = corrMean,
    group = point_id
  )
) +
  
  # grey reference line at zero
  geom_hline(
    yintercept = 0,
    colour = "grey85",
    linewidth = 0.5
  ) +
  
  # error bars
  geom_errorbar(
    aes(
      ymin = corrMean - corrSD,
      ymax = corrMean + corrSD
    ),
    width = 0.08,
    linewidth = 0.5,
    position = pd
  ) +
  
  # points
  geom_point(
    size = 2.2,
    position = pd
  ) +
  
  # three separate x-axis sections
  facet_grid(
    . ~ domain_label,
    scales = "free_x",
    space = "free_x",
    switch = "x"
  ) +
  
  # y-axis
  scale_y_continuous(
    limits = c(-0.35, 0.95),
    breaks = c(-0.3, 0, 0.3, 0.6, 0.9),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = "Inter-region correlation"
  ) +
  
  # ----------------------------------------------------------
# Formatting
# ----------------------------------------------------------

theme_classic(
  base_family = "Arial"
) +
  
  theme(
    
    # subdomain labels
    axis.text.x = element_text(
      family = "Arial",
      size = 6,
      angle = 45,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 2),
      colour = "black"
    ),
    
    # y-axis tick labels
    axis.text.y = element_text(
      family = "Arial",
      size = 6,
      colour = "black"
    ),
    
    # y-axis title
    axis.title.y = element_text(
      family = "Arial",
      size = 7,
      margin = margin(r = 5)
    ),
    
    # domain labels below x-axis
    strip.placement = "outside",
    
    strip.background = element_blank(),
    
    strip.text.x.bottom = element_text(
      family = "Arial",
      size = 6,
      face = "bold",
      margin = margin(t = 5, b = 1)
    ),
    
    # visible gaps between domains
    panel.spacing.x = unit(2, "mm"),
    
    # axes
    axis.line = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    
    # extra bottom space prevents descenders such as "g"
    # from being clipped
    plot.margin = margin(
      t = 2,
      r = 2,
      b = 5,
      l = 2,
      unit = "mm"
    )
  )


# ------------------------------------------------------------
# Show plot
# ------------------------------------------------------------

p


# ------------------------------------------------------------
# Save SVG in final size
# ------------------------------------------------------------

ggsave(
  filename = here("results", "figures", "Figure3a.svg"),
  plot = p,
  width = 17,
  height = 7,
  units = "cm"
)