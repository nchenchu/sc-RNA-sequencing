Load required libraries
library(data.table)
library(ineq)
library(ggplot2)

# Define Output Directory
output_dir <- "/Users/bigley/Library/CloudStorage/Box-Box/Bigley Lab/Navyasree Chenchu/Chenchu Projects/Results/Combined_Data/"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the combined dataset
combined_data_file <- file.path(output_dir, "Combined_Data.csv")
combined_df <- fread(combined_data_file)

# Debug: Print column names in the dataset
print("Available columns in dataset:")
print(names(combined_df))

# Function to compute the Gini coefficient
gini_coeff <- function(column) {
  values <- na.omit(column)  # Remove NA values
  if (length(values) > 1 && is.numeric(values)) {
    return(Gini(values))
  } else {
    return(NA)
  }
}

# Filter TCR (T Cells) and BCR (B Cells) based on 'cell_type' column
tcr_data <- combined_df[cell_type %like% "T Cells"]
bcr_data <- combined_df[cell_type %like% "B Cells"]

# Select numeric columns
numeric_cols <- names(combined_df)[sapply(combined_df, is.numeric)]

# Compute Gini coefficients for TCR and BCR
tcr_gini <- sapply(tcr_data[, ..numeric_cols, with=FALSE], gini_coeff)
bcr_gini <- sapply(bcr_data[, ..numeric_cols, with=FALSE], gini_coeff)

# Convert results to dataframe for visualization
gini_df <- data.table(
  Column = c(names(tcr_gini), names(bcr_gini)),
  Gini = c(tcr_gini, bcr_gini),
  Type = rep(c("TCR (T Cells)", "BCR (B Cells)"), times = c(length(tcr_gini), length(bcr_gini)))
)

# Remove NA values from the dataframe to avoid plotting issues
gini_df <- na.omit(gini_df)

# Ensure the dataframe is not empty
if (nrow(gini_df) == 0) {
  stop("No valid numeric columns found for TCR/BCR Gini computation.")
}

# Plot Gini coefficient comparison for TCR and BCR
gini_plot <- ggplot(gini_df, aes(x = reorder(Column, Gini), y = Gini, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Gini Coefficient Analysis of TCR (T Cells) and BCR (B Cells) Diversity",
       x = "Sample Type",
       y = "Gini Coefficient") +
  scale_fill_manual(values = c("blue", "red"))

# Save the plot as PNG
plot_output_file <- file.path(output_dir, "Gini_Coefficient_Analysis.png")
ggsave(filename = plot_output_file, plot = gini_plot, width = 10, height = 6, dpi = 300)

# Save Gini coefficient results as CSV
gini_output_file <- file.path(output_dir, "Gini_Coefficient_Results.csv")
fwrite(gini_df, gini_output_file)

# Print confirmation messages
message("Gini coefficient results saved to: ", gini_output_file)
message("Gini coefficient plot saved to: ", plot_output_file)








