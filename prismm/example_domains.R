# Required Libraries
library(ggplot2)
library(gridExtra)

# Required Libraries
library(combinat)

# Generate permutations of 1, 5, 2
base_alphas <- permn(c(1, 5, 2))

# Append c(1, 1, 1) to the list
base_alphas <- c(base_alphas, list(c(1, 1, 1)))

# Values to multiply
multipliers <- c(1, 10, 100)

# Function to sample from Dirichlet Distribution
sample_dirichlet <- function(alpha, n=1000) {
  matrix(sapply(1:n, function(x) {
    x <- rgamma(length(alpha), alpha)
    x / sum(x)
  }), nrow=n, byrow=TRUE)
}

# Initialize a matrix for storing plots
plots <- matrix(list(), nrow = length(multipliers), ncol = length(base_alphas), byrow = TRUE)

# Initialize a matrix for storing alpha values
alpha_values <- matrix("", nrow = length(multipliers), ncol = length(base_alphas), byrow = TRUE)

# Iterate over base_alphas (columns) and multipliers (rows) to create plots
for (i in seq_along(base_alphas)) {
  for (j in seq_along(multipliers)) {
    alpha <- base_alphas[[i]] * multipliers[j]
    df <- as.data.frame(sample_dirichlet(alpha))
    colnames(df) <- c("prob_down", "prob_up", "Var3") # Update column names
    p <- ggplot(df, aes(x = prob_down, y = prob_up)) + # Update aes_string arguments
      geom_point(aes(colour = Var3), alpha = 0.5) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = paste("Alpha =", paste(alpha, collapse = ", ")),
           x = "prob_down", # Update x label
           y = "prob_up", # Update y label
           color = "p") +
      theme_minimal() +
      xlim(0,1) +
      ylim(0,1)
    # Add the plot to the list
    plots[j, i] <- list(p)
    
    # Add the alpha value to the alpha_values matrix
    alpha_values[j, i] <- paste(alpha, collapse = "_")
  }
}

# Convert the matrix to a list for use with grid.arrange
plots <- as.list(t(plots))

# Arrange the plots
arranged_plots <- do.call(gridExtra::grid.arrange, c(plots, ncol = length(base_alphas)))

# Save the plot
ggsave("/home/users/allstaff/lmcintosh/ev_likelihood/pictures/my_plot.png",
       arranged_plots,
       width = 16, height = 10)


# Save the alpha_values matrix to a CSV file
# Convert the data frame to a vector
alpha_values_vector <- as.vector(unlist(alpha_values))

# Concatenate all the values into a single string
alpha_values_single_row <- data.frame(t(paste(alpha_values_vector, collapse = ",")))

# Write the result to a CSV
write.csv(alpha_values_single_row, "/home/users/allstaff/lmcintosh/ev_likelihood/alpha_values.csv", row.names = FALSE)
