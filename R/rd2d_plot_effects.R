# rd2d: plot point estimation and confidence intervals / bands
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu
# Last update: 2025/5/14

rm(list=ls(all=TRUE))

install.packages('rd2d')
library(rd2d)
library(ggplot2)
library(latex2exp)

################################ Load data #####################################

data_rd2d <- read.csv("Data/data_rd2d.csv")
eval <- read.csv("Data/eval.csv")
D <- as.matrix(read.table("Data/D.csv", sep = ",", header = FALSE))

y <- data_rd2d$y
X <- data_rd2d[,c("x.1", "x.2")]
t <- data_rd2d$t

result.rd2d <- rd2d(y, X, t, eval)

CI.lower.biv <- result.rd2d$results$CI.lower
CI.upper.biv <- result.rd2d$results$CI.upper
CB.lower.biv <- result.rd2d$results$CB.lower
CB.upper.biv <- result.rd2d$results$CB.upper

############## Point Estimation and Confidence Interval / Bands ################

neval <- 40
indx <- c(1:neval)
temp_plot <- ggplot() + theme_bw()

# Point estimation
df <- data.frame(indx = indx, y = result.rd2d$results$Est.p, label = rep(c("Bivariate"), each = length(indx)))
temp_plot <- temp_plot + geom_point(data = df,
                                    aes(x = indx, y = y, color = label, shape = label, fill = label, linetype = label))

# Confidence intervals
df_errorbar <- data.frame(indx = indx, ymin = CI.lower.biv, ymax = CI.upper.biv, label = "Pointwise CI")
temp_plot <- temp_plot + geom_errorbar(data = df_errorbar, aes(x = indx, ymin = ymin, ymax = ymax,
                                                               color = label, shape = label, fill = label, linetype = label))

# Confidence bands
df_ribbon <- data.frame(indx = indx, ymin = CB.lower.biv, ymax = CB.upper.biv, label = "Uniform CB")
temp_plot <- temp_plot + geom_ribbon(data = df_ribbon, aes(x = indx, ymin = ymin, ymax = ymax,
                                                           color = label, shape = label, fill = label, linetype = label), alpha = 0.1)

################################# Format Plots #################################

temp_plot <- temp_plot + xlab("Cutoffs on the Boundary") + ylab("Treatment Effect")

legend_order <- c("Bivariate", "Pointwise CI", "Uniform CB")

temp_plot <- temp_plot + scale_color_manual(
  values = c("Bivariate" = "blue", "Pointwise CI" = "blue","Uniform CB" = "dodgerblue4"),
  name = NULL, breaks = legend_order
) + scale_shape_manual(
  values = c("Bivariate" = 16, "Pointwise CI" = 124,"Uniform CB" = 0),  # 16 = filled circle, 17 = triangle
  name = NULL, breaks = legend_order
) + scale_fill_manual(
  values = c("Bivariate" = NA, "Pointwise CI" = NA, "Uniform CB" = "dodgerblue4"),
  name = NULL, breaks = legend_order
) + scale_linetype_manual(
  values = c("Bivariate" = 0, "Pointwise CI" = 5,"Uniform CB" = 0),
  name = NULL, breaks = legend_order) +
  guides(
    fill = guide_legend(order = 1),  # Keep "Uniform CB" first
    color = guide_legend(order = 1),  # Ensure everything is in one group
    shape = guide_legend(order = 1),   # Align shapes with color
    linetype = guide_legend(order = 1)
  )

temp_plot <- temp_plot + geom_vline(xintercept = c(21), color = "lightgrey", size = 1, linetype = "dotted")

temp_plot <- temp_plot +
  annotate(
    "text",
    x = 19,                # Same x as the vertical line
    y = 0.45,               # Adjust this so it appears where you want
    label = "kink",
    color = "dimgrey",
    size = 4,              # Text size
    vjust = -0.5,          # Vertical justification (pulls text below this y-value)
    fontface = "bold"
  )

# Place legend inside and adjust text sizes
temp_plot <- temp_plot + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold",
                               size = 15),
    axis.text.y = element_text(face = "bold",
                               size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

temp_plot <- temp_plot + scale_x_continuous(
  breaks = c(1,5,10,15,21,25,30,35,40),
  labels = c(TeX("$\\textbf{b}_{1}$"), TeX("$\\textbf{b}_{5}$"), TeX("$\\textbf{b}_{10}$"), 
             TeX("$\\textbf{b}_{15}$"),TeX("$\\textbf{b}_{21}$"), TeX("$\\textbf{b}_{25}$"), 
             TeX("$\\textbf{b}_{30}$"), TeX("$\\textbf{b}_{35}$"), TeX("$\\textbf{b}_{40}$"))
)

temp_plot <- temp_plot + coord_cartesian(xlim = c(1, 40), ylim = c(0.45, 1.05))

# Print the plot
print(temp_plot)

ggsave("Results/ci-and-cb.png", temp_plot, width = 6, height = 5)
