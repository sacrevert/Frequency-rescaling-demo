library(spatstat)
library(ggplot2)
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

# Calculate the frequency-weighted mean frequency (FWMF)
fwmf <- sum(intensities^2) / sum(intensities)
fwmf

# Visualize species frequencies
qplot(x, intensities, geom="line") + 
  ggtitle("True species frequencies") + 
  xlab("Species") + 
  ylab("Frequency")

# Window of the simulation
W <- owin(c(x_range[1], x_range[2]), c(y_range[1], y_range[2]))

# Simulate species occurrences
simulateSpecies <- function(intensity, W, species_id, num_species) {
  # Simulate Poisson process for a species
  pattern <- rpoispp(intensity, win = W)
  
  # Assign species ID as mark
  marks(pattern) <- as.factor(rep(species_id, pattern$n))
  return(pattern)
}

# Simulate for each species and combine
all_species <- lapply(1:num_species, function(i) simulateSpecies(intensities[i], W, i, num_species))
combined_species <- do.call(superimpose, all_species)

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

# First, ensure species_counts is structured as expected
# Assuming it's a list where each element is the result of quadratcount() for each species

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

# Use qplot to plot empirical species frequency curve
qplot(data = plot_data, x = Species, y = Frequency, geom = "line", 
      main = "Grid-based species frequency curve", xlab = "Species ID", ylab = "Frequency") +
  geom_point() # Adds points to the line plot

# Assuming 'plot_data' is your data frame with 'Species' and 'Frequency'
nls_model <- nls(Frequency ~ A / (1 + exp(-B * (Species - C))), 
                 data = plot_data, 
                 start = list(A = max(plot_data$Frequency), B = 0.1, C = median(plot_data$Species)))
# Create a new data frame for predictions
species_seq <- seq(min(plot_data$Species), max(plot_data$Species), length.out = 100)
pred_data <- data.frame(Species = species_seq)
pred_data$Predicted <- predict(nls_model, newdata = pred_data)
# Original data plot with qplot
p <- qplot(data = plot_data, x = Species, y = Frequency, geom = "point", 
           main = "Grid-based species frequency curve with logistic curve fit", xlab = "Species ID", ylab = "Frequency")

# Add NLS fitted line
p + geom_line(data = pred_data, aes(x = Species, y = Predicted), color = "blue")


########################################################################################################
## Better plot for "counts" object
# Assuming 'counts' is from quadratcount and you have 'x_breaks' and 'y_breaks'
# Calculate midpoints for x and y breaks to use as coordinates
x_midpoints <- (x_breaks[-length(x_breaks)] + x_breaks[-1]) / 2
y_midpoints <- (y_breaks[-length(y_breaks)] + y_breaks[-1]) / 2

# Expand grid to create a data frame of all possible combinations of x and y midpoints
grid_df <- expand.grid(x = x_midpoints, y = y_midpoints)

# Flatten counts to a vector (assuming counts is a matrix with species as columns if multiple species were counted separately)
sppPaPerGridCell <- list()
sppPaPerGridCell <- lapply(species_counts, function(counts) {
  # If counts is directly an atomic vector, sum it directly
  if (is.atomic(counts)) {
    out <- sapply(counts, function(x) ifelse(x>=1,1,0))
    return(out)
  } else if (is.list(counts) && !is.null(counts$counts)) {
    # If counts is a list and has a $counts component, sum that
    out2 <- sapply(counts$counts, function(x) ifelse(x>=1,1,0))
    #return(sum(counts$counts))
    return(out2)
  } else {
    # Otherwise, return NA or handle as appropriate for your data structure
    return(NA)
  }
})

counts_vector <- colSums(do.call(rbind, sppPaPerGridCell))

# Add counts to the grid data frame
grid_df$counts <- counts_vector

# viridis option
library(viridis)
ggplot(grid_df, aes(x = x, y = y, fill = counts)) + 
  geom_tile() +
  scale_fill_viridis_c(name = "Counts", option = "C") +
  labs(title = "Heatmap of species richness", x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal()
