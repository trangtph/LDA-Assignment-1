# import sas dataset --------------------------------------------------------
library(haven)
library(tidyverse)

hemo <- read_sas("hemodialysis.sas7bdat")

# remove NA hb, then remove patients with only 1 hb
hemo <- hemo %>% 
  filter(!is.na(hb)) %>% 
  group_by(ID) %>% 
  filter(n() > 2)

# count distinct patients
n_distinct(hemo$ID) # 2933 patients having more than 1 hb
# 2089 patients having more than 2 hb

# test linear trend ---------------------------------------------------------

## fit linear model for each patient ----------------------------------------
# Initialize lists to store statistics for each patient
r2 <- list()
SSE <- list()
SSR <- list()
ni <- list()
intercept <- list()
slope <- list()
SEintercept <- list()
SEslope <- list()

# Loop through each unique patient ID
for (i in unique(hemo$ID)) {
  # Subset data for the current patient
  patient_data <- hemo[hemo$ID == i, ]
  
  # Fit the linear model for the current patient
  Model1 <- lm(hb ~ month, data = patient_data, x = TRUE)
  
  # Calculate and store the statistics
  r2[[i]] <- summary(Model1)$r.squared
  SSE[[i]] <- sum(Model1$residuals^2)
  SSR[[i]] <- anova(Model1)$`Sum Sq`[1]
  ni[[i]] <- length(Model1$x[, 1])
  intercept[[i]] <- Model1$coefficients[1]
  slope[[i]] <- Model1$coefficients[2]
  SEintercept[[i]] <- summary(Model1)$coefficients[1, 2]
  SEslope[[i]] <- summary(Model1)$coefficients[2, 2]
}

# Calculate the meta R-squared statistic
R2_meta <- sum(unlist(SSR)) / sum(unlist(SSR) + unlist(SSE))

# Display the meta R-squared result
R2_meta

## plotting ------------------------------------------------------------------

# Convert lists to numeric vectors for plotting
r2_values <- unlist(r2)
ni_values <- unlist(ni)

# Create a data frame for ggplot
plot_data <- data.frame(ni = ni_values, r2 = r2_values)

# Plot using ggplot2
ggplot(plot_data, aes(x = ni, y = r2)) +
  geom_point(color = "blue", size = 2, alpha = 0.7) +       # Scatter plot points
  geom_hline(yintercept = R2_meta, color = "red", linetype = "dashed", size = 1) + # Horizontal line for R2_meta
  labs(title = "The linear function with time seems to describe the data reasonably well",
       x = expression(Number~n[i]~of~measurements),
       y = expression(Coefficient~R[i]^2)) +
  theme_minimal() +                                          # Minimal theme for clean look
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Center and style title
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = max(plot_data$ni), y = R2_meta,       # Label for R2_meta line
           label = expression(R[meta]^2),
           hjust = 1.1, vjust = -0.5, color = "red", size = 4)

# save plot
ggsave("Results/r-squared-linear.png", width = 8, height = 6, dpi = 300)

# test quadratic trend ------------------------------------------------------
# Initialize lists to store statistics for each patient
r2 <- list()
SSE <- list()
SSR <- list()
ni <- list()
intercept <- list()
slope <- list()
SEintercept <- list()
SEslope <- list()

# Loop through each unique patient ID
for (i in unique(hemo$ID)) {
  # Subset data for the current patient
  patient_data <- hemo[hemo$ID == i, ]
  
  # Fit the quadratic model (including month squared term) for the current patient
  Model1 <- lm(hb ~ month + I(month^2), data = patient_data, x = TRUE)
  
  # Calculate and store the statistics
  r2[[i]] <- summary(Model1)$r.squared
  SSE[[i]] <- sum(Model1$residuals^2)
  SSR[[i]] <- anova(Model1)$`Sum Sq`[1]
  ni[[i]] <- length(Model1$x[, 1])
  intercept[[i]] <- Model1$coefficients[1]
  slope[[i]] <- Model1$coefficients[2]
  SEintercept[[i]] <- summary(Model1)$coefficients[1, 2]
  SEslope[[i]] <- summary(Model1)$coefficients[2, 2]
}

# Calculate the meta R-squared statistic
R2_meta <- sum(unlist(SSR)) / sum(unlist(SSR) + unlist(SSE))

# Display the meta R-squared result
R2_meta

# plotting ------------------------------------------------------------------

# Convert lists to numeric vectors for plotting
r2_values <- unlist(r2)
ni_values <- unlist(ni)

# Create a data frame for ggplot
plot_data <- data.frame(ni = ni_values, r2 = r2_values)

# Plot using ggplot2
ggplot(plot_data, aes(x = ni, y = r2)) +
  geom_point(color = "blue", size = 2, alpha = 0.7) +       # Scatter plot points
  geom_hline(yintercept = R2_meta, color = "red", linetype = "dashed", size = 1) + # Horizontal line for R2_meta
  labs(
    title = "Quadratic Model of hb vs. Month per Patient",
    x = expression(Number~n[i]~of~measurements),
    y = expression(Coefficient~R[i]^2)
  ) +
  theme_minimal() +                                          # Minimal theme for clean look
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), # Center and style title
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = max(plot_data$ni), y = R2_meta,       # Label for R2_meta line
           label = expression(R[meta]^2),
           hjust = 1.1, vjust = -0.5, color = "red", size = 4)

# Save plot
ggsave("Results/r-squared-quadratic.png", width = 8, height = 6, dpi = 300)


