# Function to install and load packages
install_and_load_packages <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

# List of packages to be installed and loaded
packages <- c("tidyverse", "stringr", "gridExtra")

# Use the function to install and load packages
install_and_load_packages(packages)

# Load data
data <- read.csv("/home/users/allstaff/lmcintosh/ev_likelihood/summary2.csv")
print("Data loaded")

# Extract lam and alpha values from the filename
data <- data %>%
  mutate(lam = str_extract(filename, "(?<=lam)[0-9.]+"),
         alpha = str_extract(filename, "(?<=alpha)[0-9._]+")) %>%
  separate(alpha, into = c("alpha1", "alpha2", "alpha3"), sep = "_", extra = "drop")

# Convert lam to numeric
data$lam <- as.numeric(data$lam)

print("Data prepared")
print(head(data))

# Determine overall y-axis limits
y_min <- 0
y_max <- max(c(data$proportion_estimated, data$best_estimate_accuracy, data$arbitrary_estimate_accuracy), na.rm = TRUE)

# Read alpha_values.csv to get the order of alphas
alpha_order <- read.csv("/home/users/allstaff/lmcintosh/ev_likelihood/alpha_values.csv", stringsAsFactors = FALSE)
alpha_order <- strsplit(alpha_order[[1]], ",")[[1]]
# Reshape alpha_order into a matrix, fill it by column, then turn it back into a vector by concatenating the rows
nrows <- length(multipliers)  # assuming multipliers has been defined before this script
alpha_order_matrix <- matrix(alpha_order, nrow = nrows, byrow = FALSE)
alpha_order <- as.vector(t(alpha_order_matrix))

# Generate separate plots for each alpha triplet
alphas <- unique(paste(data$alpha1, data$alpha2, data$alpha3, sep = "_"))
alphas <- alphas[alphas != "NA_NA_NA"]  # Exclude "NA_NA_NA" alpha values

# Order alphas according to the order in alpha_values.csv
alphas <- alphas[match(alpha_order, alphas)]

print("Ordered alphas found")
print(alphas)

# Create two empty lists to store plots
plot_list_best <- list()
plot_list_arbitrary <- list()

for(i in seq_along(alphas)) {
  alpha <- alphas[i]
  
  plot_data <- data %>% filter(paste(alpha1, alpha2, alpha3, sep = "_") == alpha)
  print(paste("Plotting data for alpha: ", alpha))
  print(head(plot_data))
  
  # Reshape data for bar plots
  plot_data_long_best <- plot_data %>% select(lam, proportion_estimated, best_estimate_accuracy) %>%
    gather(key = "measure", value = "value", proportion_estimated, best_estimate_accuracy)
  
  plot_data_long_arbitrary <- plot_data %>% select(lam, proportion_estimated, arbitrary_estimate_accuracy) %>%
    gather(key = "measure", value = "value", proportion_estimated, arbitrary_estimate_accuracy) %>%
    mutate(measure = if_else(measure == "arbitrary_estimate_accuracy", "PCAWG_cuttoff_accuracy", measure))
  
  # Plot for best_estimate_accuracy
  plot_best <- ggplot(plot_data_long_best, aes(x = lam, y = value, fill = measure)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = alpha,
         x = "Avg. CN events",
         y = "Accuracy") +
    theme_bw() +
    coord_cartesian(ylim = c(y_min, y_max)) +  # Add this line to set the y-axis limits
    theme(legend.position = "none")  # Hide legend in individual plots
  
  # Plot for arbitrary_estimate_accuracy
  plot_arbitrary <- ggplot(plot_data_long_arbitrary, aes(x = lam, y = value, fill = measure)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = alpha,
         x = "Avg. CN events",
         y = "Accuracy") +
    theme_bw() +
    coord_cartesian(ylim = c(y_min, y_max)) +  # Add this line to set the y-axis limits
    theme(legend.position = "none")  # Hide legend in individual plots
  
  # Add the plots to the lists
  plot_list_best[[i]] <- plot_best
  plot_list_arbitrary[[i]] <- plot_arbitrary
}

# Create a separate legend plot
legend_plot <- ggplot(data = data.frame(x = factor(1:3), y = 0, fill = c("proportion_estimated", "best_estimate_accuracy", "PCAWG_cuttoff_accuracy")), aes(x = x, y = y, fill = fill)) +
  geom_bar(stat = "identity") +
  theme_void() +  # No axes or grid
  theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 6))  # Legend at the bottom, with smaller text

# Save the legend plot
ggsave("/home/users/allstaff/lmcintosh/ev_likelihood/pictures/legend.png", legend_plot, width = 10, height = 2)  # Adjusted width and height for 5:1 aspect ratio
print(legend_plot)
print("Legend plot saved")

# Arrange plots in two grids
grid_plot_best <- gridExtra::grid.arrange(grobs = plot_list_best, nrow = 3, top = "Accuracy under different dirichelet priors: best_estimate_accuracy")
grid_plot_arbitrary <- gridExtra::grid.arrange(grobs = plot_list_arbitrary, nrow = 3, top = "Accuracy under different dirichelet priors: PCAWG_cuttoff_accuracy")

# Save the grid plots
ggsave("/home/users/allstaff/lmcintosh/ev_likelihood/pictures/summary_bar_plot_best.png", grid_plot_best, width = 16, height = 10)  # Adjusted width and height for 16:10 aspect ratio
ggsave("/home/users/allstaff/lmcintosh/ev_likelihood/pictures/summary_bar_plot_arbitrary.png", grid_plot_arbitrary, width = 16, height = 10)  # Adjusted width and height for 16:10 aspect ratio
print(grid_plot_best)
print(grid_plot_arbitrary)
print("Done")
