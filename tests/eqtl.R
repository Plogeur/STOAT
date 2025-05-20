# Load command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript linear_regression_from_files.R X_file.txt Y_file.txt")
}

x_file <- args[1]
y_file <- args[2]

# Read the data
X <- scan(x_file, what = numeric(), quiet = TRUE)
Y <- scan(y_file, what = numeric(), quiet = TRUE)

# Check length
if (length(X) != length(Y)) {
  stop("Error: X and Y must have the same number of entries")
}

# Fit linear model
model <- lm(Y ~ X)

# Print summary
summary(model)

# Rscript eqtl.R X_eqtl.tsv Y_eqtl.tsv