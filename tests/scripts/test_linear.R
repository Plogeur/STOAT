# Load necessary library
library(car)     # For VIF
library(ggplot2) # For plots

# Generate a single valid row that sums to 1 with elements in {0, 0.5, 1}
generate_valid_row <- function() {
  repeat {
    row <- sample(c(0, 0.5, 1), 4, replace = TRUE)
    if (sum(row) == 1) return(row)
  }
}

# Generate dataset: 20 rows of valid samples
generate_dataset <- function(n_rows = 20) {
  X <- t(replicate(n_rows, generate_valid_row()))
  colnames(X) <- c("Intercept", "X1", "X2", "X3")
  return(X)
}

# Run a simulation comparing C vs R implementation
run_simulation <- function(n_simulations = 100) {
  results <- data.frame()

  for (i in 1:n_simulations) {
    X <- generate_dataset()
    beta_true <- c(2, -1, 0.5, 3)  # True coefficients
    noise <- rnorm(nrow(X), mean = 0, sd = 0.1)
    y <- X %*% beta_true + noise

    # Fit with R's linear model
    df <- as.data.frame(X)
    df$y <- y
    model_r <- lm(y ~ X1 + X2 + X3, data = df)

    # Simulate "C model" output (pretend it has small noise/errors)
    # Replace this with your actual C outputs
    coef_c <- coef(model_r) + rnorm(length(coef(model_r)), sd = 0.05)

    # Multicollinearity via VIF
    vif_values <- vif(model_r)

    # Store results
    results <- rbind(results, data.frame(
      sim = i,
      R_Intercept = coef(model_r)[1],
      R_X1 = coef(model_r)[2],
      R_X2 = coef(model_r)[3],
      R_X3 = coef(model_r)[4],
      C_Intercept = coef_c[1],
      C_X1 = coef_c[2],
      C_X2 = coef_c[3],
      C_X3 = coef_c[4],
      VIF_X1 = vif_values[1],
      VIF_X2 = vif_values[2],
      VIF_X3 = vif_values[3]
    ))
  }

  return(results)
}

# Run it!
set.seed(42)
sim_results <- run_simulation(10000)

# -------------------------------
# 📊 Plotting Comparisons
# -------------------------------

# Melt data for plotting
library(reshape2)
coefs <- melt(sim_results[, c("sim", "R_X1", "R_X2", "R_X3", "C_X1", "C_X2", "C_X3")], id.vars = "sim")
coefs$impl <- ifelse(grepl("^R_", coefs$variable), "R", "C")
coefs$feature <- gsub("^(R_|C_)", "", coefs$variable)

# Plot: Coefficient comparison by implementation
ggplot(coefs, aes(x = feature, y = value, fill = impl)) +
  geom_boxplot(position = "dodge") +
  labs(title = "Coefficient Estimates: R vs C Implementation",
       y = "Estimated Coefficient", x = "Feature") +
  theme_minimal()

# Plot: VIF values
vifs <- melt(sim_results[, c("sim", "VIF_X1", "VIF_X2", "VIF_X3")], id.vars = "sim")
ggplot(vifs, aes(x = variable, y = value)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "Variance Inflation Factor (VIF) per Feature",
       x = "Feature", y = "VIF") +
  theme_minimal()

# Plot: Coefficient Differences (R - C)
diffs <- data.frame(
  sim = sim_results$sim,
  Diff_X1 = sim_results$R_X1 - sim_results$C_X1,
  Diff_X2 = sim_results$R_X2 - sim_results$C_X2,
  Diff_X3 = sim_results$R_X3 - sim_results$C_X3
)
diffs_melted <- melt(diffs, id.vars = "sim")
ggplot(diffs_melted, aes(x = variable, y = value)) +
  geom_boxplot(fill = "orange") +
  labs(title = "Coefficient Differences (R - C)", y = "Difference", x = "Feature") +
  theme_minimal()
