# Required Libraries
library(ggplot2)
library(gridExtra)

# Function to sample from Dirichlet Distribution
sample_dirichlet <- function(alpha, n=1000) {
  matrix(sapply(1:n, function(x) {
    x <- rgamma(length(alpha), alpha)
    x / sum(x)
  }), nrow=n, byrow=TRUE)
}

# Loop for asking user input and generating the plot
while(TRUE) {
  # Ask user for alpha values
  alpha <- readline(prompt="Enter three alpha values (separated by space) or 'quit' to exit: ")
  
  # If user enters 'quit', break the loop and exit
  if(tolower(alpha) == "quit"){
    break
  }
  
  alpha <- as.numeric(unlist(strsplit(alpha, " ")))
  
  if(length(alpha) != 3){
    print("Please enter exactly three alpha values. Try again.")
  } else {
    # Generate the data frame from user defined alpha
    df <- as.data.frame(sample_dirichlet(alpha))
    colnames(df) <- paste0("Var", 1:length(alpha))
    
    # Generate the plot for the user defined alpha
    plot <- ggplot(df, aes_string(x = names(df)[1], y = names(df)[2])) +
      geom_point(aes_string(colour = names(df)[3]), alpha = 0.5) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = paste("Alpha =", paste(alpha, collapse = ", ")),
           x = names(df)[1],
           y = names(df)[2],
           color = names(df)[3]) +
      xlim(0, 1) + 
      ylim(0, 1) +
      theme_minimal()
    
    print(plot)
  }
}
