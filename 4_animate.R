## Produce animated example of how a multiplier adjusts a set of species frequencies relative to some target frequency-weighted mean frequency (Phi in Frescalo)
## Oli Pescott
## 04th March 2024
#rm(list=ls())
library(spatstat)
library(ggplot2)
library(magick) ## gganimate didn't work without ggplotly, which I couldn't install due to some version conflict
set.seed(92)

# Parameters
x_range <- c(0, 20)
y_range <- c(0, 20)
grid_size <- 10
num_species <- 100
max_presence <- 0.98
# Generate species frequencies using a sigmoidal curve
x <- 1:num_species
intensities <- max_presence / (1 + exp(-0.075 * (x - num_species / 2)))
calculateFWMF <- function(intensities) {
  sum(intensities^2) / sum(intensities)
}
RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}

calculateFWMF(intensities) # "pure" Phi based on underlying Poisson intensities

simulateSpeciesOccurrences <- function(intensities = intensities, x_range = x_range, y_range = y_range, num_species = num_species) {
  W <- owin(c(x_range[1], x_range[2]), c(y_range[1], y_range[2]))
  all_species <- lapply(1:num_species, function(i) {
    pattern <- rpoispp(intensities[i], win = W)
    marks(pattern) <- as.factor(rep(i, pattern$n))
    return(pattern)
  })
  combined_species <- do.call(superimpose, all_species)
  return(combined_species)
}

combined_species <- simulateSpeciesOccurrences(intensities = intensities, x_range = x_range, y_range = y_range, num_species = num_species)

plotGriddedSpeciesFrequency <- function(combined_species = combined_species, x_range = x_range, y_range = y_range, grid_size = grid_size) {
  # Create breaks for the grid cells
  x_breaks <- seq(x_range[1], x_range[2], length.out = grid_size + 1)
  y_breaks <- seq(y_range[1], y_range[2], length.out = grid_size + 1)
  
  # Count over all occurrences of every species in each grid cell
  counts <- quadratcount(combined_species, tess = tess(xgrid = x_breaks, ygrid = y_breaks))
  
  # Now assuming 'combined_species' is a ppp object with marks indicating species
  # Split the combined point pattern into a list of point patterns for each species
  species_ppps <- split.ppp(combined_species)
  
  # Use lapply to apply quadratcount to each species' point pattern
  species_counts <- lapply(species_ppps, function(ppp) {
    quadratcount(ppp, tess = tess(xgrid = x_breaks, ygrid = y_breaks))
  })
  
  # Extract total counts per species correctly
  total_counts_per_species <- sapply(species_counts, function(counts) {
    # If counts is directly an atomic vector, sum it directly
    if (is.atomic(counts)) {
      out <- sapply(counts, function(x) ifelse(x>=1,1,0))
      return(sum(out))
    } else if (is.list(counts) && !is.null(counts$counts)) {
      # If counts is a list and has a $counts component, sum that
      out2 <- sapply(counts$counts, function(x) ifelse(x>=1,1,0))
      #return(sum(counts$counts))
      return(sum(out2))
    } else {
      # Otherwise, return NA or handle as appropriate for your data structure
      return(NA)
    }
  })
  
  # Calculate species grid/quadrat frequencies
  species_frequenciesQ <- total_counts_per_species / grid_size^2
  
  # Prepare data for plotting
  plot_data <- data.frame(Species = 1:length(species_frequenciesQ), Frequency = species_frequenciesQ)
  
  # Assuming 'plot_data' is your data frame with 'Species' and 'Frequency'
  nls_model <- nls(Frequency ~ A / (1 + exp(-B * (Species - C))), 
                   data = plot_data, 
                   start = list(A = max(plot_data$Frequency), B = 0.1, C = median(plot_data$Species)))
  # Create a new data frame for predictions
  species_seq <- seq(min(plot_data$Species), max(plot_data$Species), length.out = 100)
  pred_data <- data.frame(Species = species_seq)
  pred_data$Predicted <- predict(nls_model, newdata = pred_data)
  # Original data plot with qplot
  p <- ggplot2::qplot(data = plot_data, x = Species, y = Frequency, geom = "point", ylim = c(0,2.1),
             main = "Grid-based species frequency curve with logistic curve fit: Target Phi = 0.84", xlab = "Species ID", ylab = "Frequency")
  
  # Add NLS fitted line
  p <- p + geom_line(data = pred_data, aes(x = Species, y = Predicted), color = "blue")
  
  # Empirical FWMF based on grid-based summary at given scale
  fwmf <- calculateFWMF(species_frequenciesQ)
  
  median_species_id <- median(plot_data$Species)
  
  p <- p + geom_point(aes(x = median_species_id, y = fwmf), color = "red", size = 3) +
    annotate("text", x = 1, y = 1.8, label = "Perfect!", 
             hjust = 0, vjust = 0, color = "red", size = 5)
  
    # Save the plot as an image
  #ggsave(filename = paste0("outputs/plot_", as.character(round(fwmf,2)), ".png"), plot = p, width = 10, height = 8)
  ggsave(filename = paste0("outputs/plot_", "009", ".png"), plot = p, width = 8, height = 6)
  # also return items
  return(list(plot = p, fwmf = fwmf, gridFreqs = species_frequenciesQ))
}

out <- plotGriddedSpeciesFrequency(combined_species = combined_species, x_range = x_range, y_range = y_range, grid_size = grid_size)

# Assuming original grid frequencies are stored in `original_intensities`
# and the original grid FWMF is calculated and stored in `target_fwmf`

findMultiplier <- function(original_intensities = out$gridFreqs/3, target_fwmf = 0.84, tolerance = 0.01, max_iterations = 10) {
  low <- 0
  high <- 20
  current_multiplier <- 1
  iterList <- list()
  for (i in 1:max_iterations) {
    current_intensities <- original_intensities * current_multiplier
    iterList[[i]] <- current_intensities
    current_fwmf <- calculateFWMF(current_intensities)
    if (abs(current_fwmf - target_fwmf) < tolerance) {
      break
    } else if (current_fwmf < target_fwmf) {
      low <- current_multiplier
    } else {
      high <- current_multiplier
    }
    current_multiplier <- (low + high) / 2
  }
  current_multiplier <- round(current_multiplier,2)
  return(c(final_multiplier = current_multiplier, iterList = iterList))
}

out2 <- findMultiplier(original_intensities = out$gridFreqs/5, target_fwmf = 0.84, tolerance = 0.01, max_iterations = 10)

# We now have a list of interated grid frequencies that we can plot

#generatePlot <- function(grid_freqs, iter) {
generatePlot <- function(grid_freqs) {
  # Prepare data for plotting
  plot_data <- data.frame(Species = 1:length(grid_freqs), Frequency = grid_freqs)
  
  # Assuming 'plot_data' is your data frame with 'Species' and 'Frequency'
  nls_model <- nls(Frequency ~ A / (1 + exp(-B * (Species - C))), 
                   data = plot_data, 
                   start = list(A = max(plot_data$Frequency), B = 0.1, C = median(plot_data$Species)))
  # Create a new data frame for predictions
  species_seq <- seq(min(plot_data$Species), max(plot_data$Species), length.out = 100)
  pred_data <- data.frame(Species = species_seq)
  pred_data$Predicted <- predict(nls_model, newdata = pred_data)
  # Original data plot with qplot
  p <- ggplot2::qplot(data = plot_data, x = Species, y = Frequency, geom = "point", ylim = c(0,2.1),
                      main = "Grid-based species frequency curve with logistic curve fit: Target Phi = 0.84", xlab = "Species ID", ylab = "Frequency")
  
  # Add NLS fitted line
  p <- p + geom_line(data = pred_data, aes(x = Species, y = Predicted), color = "blue")
  
  # Empirical FWMF based on grid-based summary at given scale
  fwmf <- calculateFWMF(grid_freqs)
  
  median_species_id <- median(plot_data$Species)
  
  # p <- p + geom_point(aes(x = median_species_id, y = fwmf), color = "red", size = 3) +
  #   annotate("text", x = 1, y = 1.8, label = bquote(phi == .(round(fwmf, 2))), 
  #            hjust = 0, vjust = 0, color = "red", size = 5)
  
    # Find the corresponding frequency for the median species for a more accurate y-coordinate
  median_species_freq <- plot_data$Frequency[which.min(abs(plot_data$Species - median_species_id))]
  
  p <- p + geom_point(aes(x = median_species_id, y = median_species_freq), color = "red", size = 3) +
           annotate("text", x = 1, y = 1.8, label = bquote(phi == .(round(fwmf, 2))), 
                  hjust = 0, vjust = 0, color = "red", size = 5)
  
  # Add iteration text
  # number <- RIGHT(iter, 1)
  # p <- p + annotate("text", x = 1, y = 2.05, label = paste0("Iteration: ", as.character(number)), 
  #                   hjust = 0, vjust = 0, color = "black", size = 5)
    
   # Save the plot as an image
  #ggsave(filename = paste0("outputs/plot_", as.character(round(fwmf,3)), ".png"), plot = p, width = 10, height = 8)
  #ggsave(filename = sprintf("outputs/plot_%03d.png", plot_id), plot = p, width = 8, height = 6)
  
  # also
  return(p)
}

#out3 <- lapply(out2[2:9], function(x) generatePlot(grid_freqs = x, iter = names(x)))
out3 <- lapply(out2[2:9], function(x) generatePlot(grid_freqs = x))

# Iterate through the list of plots and save each
for (i in seq_along(out3)) {
  # Define a unique filename for each plot
  plotname <- sprintf("outputs/plot_%03d.png", i) # This names files as plot_001.png, plot_002.png, etc.
  
  # Save the current plot in the loop
  ggsave(filename = plotname, plot = out3[[i]], width = 8, height = 6)
}

# Now animate using magick
# Read the images
image_files <- list.files(path = "outputs/", pattern = "plot_\\d+\\.png", full.names = TRUE)
image_list <- lapply(image_files, image_read)
image_combined <- image_join(image_list)

# Animate and save
animation <- image_animate(image_combined, fps = 1)
image_write(animation, "outputs/species_frequency_animation.gif")


