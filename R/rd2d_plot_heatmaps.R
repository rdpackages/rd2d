# rd2d: plot heatmaps
# Authors: M. D. Cattaneo, R. Titiunik, R. R. Yu
# Last update: 2025/5/14

rm(list=ls(all=TRUE))

install.packages('rd2d')

library(rd2d)
library(ggplot2)
library(latex2exp)
library(dplyr)

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

############################### Plot heatmaps ##################################

# heat map for treatment effect
data.plot <- cbind(eval, result.rd2d$results$Est.p)
colnames(data.plot) <- c("x.1", "x.2", "tau.hat")

# Function to interpolate points and color values between two consecutive points
interpolate_points <- function(df, n=ninter){
  do.call(rbind, lapply(1:(nrow(df)-1), function(i){
    xseq <- seq(df$x.1[i], df$x.1[i+1], length.out = n+2)[2:(n+1)]
    yseq <- seq(df$x.2[i], df$x.2[i+1], length.out = n+2)[2:(n+1)]
    colorseq <- seq(df$tau.hat[i], df$tau.hat[i+1], length.out = n+2)[2:(n+1)]
    data.frame(x.1 = xseq, x.2 = yseq, tau.hat = colorseq)
  }))
}

# Generating interpolated points
ninter <- 10
interpolated_data <- interpolate_points(data.plot)

# Plotting
heat_wd <- 6.5
heatcol_low <- "blue"
heatcol_high <- "red"
heat_lab <- TeX("$\\tau(\\textbf{b})$")
heat_title <- NULL
xlabel <- TeX("$X_{1i}$")
ylabel <- TeX("$X_{2i}$")

augmented_data <-rbind(data.plot, interpolated_data)
ord <- order(augmented_data[,1], augmented_data[,2])
augmented_data <- augmented_data[ord,]

if (is.null(heat_title)) heat_title <- "Heat Map"
plot_heat <- ggplot(data.plot, aes(x=x.1, y=x.2)) +
  geom_segment(data = augmented_data, aes(xend = lead(x.1, order_by=x.1), yend = lead(x.2, order_by=x.2), color = tau.hat),
               size = heat_wd, lineend = "round") +
  scale_color_gradient(low=heatcol_low, high=heatcol_high) +
  labs(color = heat_lab) +
  xlab(xlabel) +
  ylab(ylabel)

plot_heat <- plot_heat + geom_text(data=data.plot, aes(x=x.1, y=x.2, label=sprintf("%02d", 1:nrow(data.plot))), color="black", size=3)

plot_heat <- plot_heat + theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15,face = "bold"),  # X-axis label size
    axis.title.y = element_text(size = 15,face = "bold"),  # Y-axis label size
    plot.title = element_text(size = 20, hjust = 0.5),  # Title size and centering
    text=element_text(family="Times New Roman", face="bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),  # White background, no border
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + coord_cartesian(xlim = c(-10, 55), ylim = c(-10, 55))

print(plot_heat)

ggsave("Results/heat-treatment_effect.png", plot_heat, width = 6, height = 5)

# heatmap for p-value

data.plot$p.value <- result.rd2d$results$`P>|z|`
data.plot$p.sig <- cut(data.plot$p.value,
                       breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       labels = c("p < 0.001", "0.001 ≤ p < 0.01", "0.01 ≤ p < 0.05", "0.05 ≤ p < 0.1", "p ≥ 0.1"))
sig_colors <- c("p < 0.001" = "#d73027",       # red
                "0.001 ≤ p < 0.01" = "#fc8d59", # orange
                "0.01 ≤ p < 0.05" = "#fee08b",  # yellow
                "0.05 ≤ p < 0.1" = "#d9ef8b",   # light green
                "p ≥ 0.1" = "#91cf60")          # green
library(ggplot2)

plot_heat_pvalue <- ggplot(data.plot, aes(x = x.1, y = x.2, fill = p.sig)) +
  geom_tile(color = "white",show.legend = TRUE) +
  scale_fill_manual(values = sig_colors, name = "P-value",drop = FALSE) +
  labs(x = TeX("$X_{1i}$"), y = TeX("$X_{2i}$")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5),
    text = element_text(family = "Times New Roman", face = "bold"),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.position = c(0.8, 1),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

plot_heat_pvalue <- plot_heat_pvalue +  coord_cartesian(xlim = c(-10, 55), ylim = c(-10, 55)) 
plot_heat_pvalue <- plot_heat_pvalue + geom_text(data=data.plot, aes(x=x.1, y=x.2, label=sprintf("%02d", 1:nrow(data.plot))), color="black", size=3)

print(plot_heat_pvalue)

ggsave("Results/heat-pvalue.png", plot_heat_pvalue, width = 6, height = 5)

