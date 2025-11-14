############################################
############################################
# load genotype data
genoGE <- read.csv("~/Downloads/GenotypesGxE.csv", row.names=1)

#load data, ran one year at a time
GxE23forModelPerms <- read.csv("~/Downloads/GxE23forModelPerms.csv")
all<-GxE23forModelPerms

GxE22forModelPerms <- read.csv("~/GxE22forModelPerms.csv")
all<- GxE22forModelPerms


############################################
### Define response variables 
# 2022: 6:49
# 2023: 6:47
############################################

#change this to correspond to the correct year data
response_vars <- colnames(all)[6:47]

# Keep only numeric columns
numeric_mask  <- sapply(all[response_vars], is.numeric)
response_vars <- response_vars[numeric_mask]

if (length(response_vars) == 0) {
  stop("No numeric response variables found in columns 6:49 of 'all'")
}

############################################
### Match genotypes between all and genoGE
############################################

# Get genotype names from genoGE (columns) and all (Genotype column)
genoGE_genotypes <- colnames(genoGE)
all_genotypes    <- all$Genotype

# Find common genotypes
common_genotypes <- intersect(genoGE_genotypes, all_genotypes)
if (length(common_genotypes) == 0) {
  stop("No matching genotypes found between 'all' and 'genoGE'")
}

message("Found ", length(common_genotypes), " common genotypes")

# Subset and align data
all_matched      <- all[all$Genotype %in% common_genotypes, ]
genoGE_matched   <- genoGE[, common_genotypes, drop = FALSE]

# Order both data frames consistently
all_matched    <- all_matched[match(common_genotypes, all_matched$Genotype), ]
genoGE_matched <- genoGE_matched[, match(common_genotypes, colnames(genoGE_matched)), drop = FALSE]

# Verify dimensions
message("Matching genotypes: ", paste(common_genotypes, collapse = ", "))
message("'all' dimensions: ", nrow(all_matched), " rows")
message("'genoGE' dimensions: ", ncol(genoGE_matched), " columns")

############################################
### Initialize containers
############################################

set.seed(6666)  # optional, for reproducibility
n_iter <- 100

# For full per-test results across all traits
all_results_list <- vector("list", length(response_vars))
names(all_results_list) <- response_vars

# For per-trait min p and empirical threshold
trait_summary <- data.frame(
  ResponseVar = response_vars,
  MinPvalue   = NA_real_,
  EmpThresh95 = NA_real_,
  stringsAsFactors = FALSE
)

############################################
### Main loop over response variables
############################################

for (rv in response_vars) {
  message("\n=== Working on response variable: ", rv, " ===")
  
  # Initialize results for this response variable
  results <- data.frame(
    Variable = character(n_iter),
    Pvalue   = numeric(n_iter),
    stringsAsFactors = FALSE
  )
  
  # Perform n_iter random regressions
  for (i in seq_len(n_iter)) {
    # Randomly sample one row from genoGE (each row is a variable)
    sampled_row      <- sample(nrow(genoGE_matched), 1)
    sampled_var_name <- rownames(genoGE_matched)[sampled_row]
    sampled_values   <- unlist(genoGE_matched[sampled_row, ])
    
    # Check if the sampled variable has variation
    if (length(unique(sampled_values)) < 2) {
      message("Skipping ", sampled_var_name, " - no variation in values")
      results$Variable[i] <- sampled_var_name
      results$Pvalue[i]   <- NA
      next
    }
    
    # Perform linear regression
    model <- tryCatch(
      {
        lm(all_matched[[rv]] ~ sampled_values)
      },
      error = function(e) {
        message("Error with variable ", sampled_var_name, ": ", e$message)
        return(NULL)
      }
    )
    
    # Store results
    results$Variable[i] <- sampled_var_name
    if (!is.null(model)) {
      results$Pvalue[i] <- summary(model)$coefficients[2, 4]
    } else {
      results$Pvalue[i] <- NA
    }
    
    # Print progress
    if (i %% 10 == 0) {
      message("Completed ", i, "/", n_iter, " iterations for ", rv)
    }
  }
  
  # Clean and tag with response variable name
  results <- na.omit(results)
  results$ResponseVar <- rv
  
  # Store in list
  all_results_list[[rv]] <- results
  
  ##########################################
  ### Per-trait min p + empirical threshold
  ##########################################
  
  # Minimum p for this trait
  min_p_rv <- min(results$Pvalue, na.rm = TRUE)
  
  # Empirical distribution of min p-values for this trait
  n_draw <- min(100, nrow(results))  # in case fewer than 100 rows
  min_p_distribution_rv <- replicate(100, {
    min(sample(results$Pvalue, n_draw, replace = TRUE), na.rm = TRUE)
  })
  
  threshold_rv <- quantile(min_p_distribution_rv, 0.95, na.rm = TRUE)
  
  # Store in trait_summary
  trait_summary[trait_summary$ResponseVar == rv, "MinPvalue"]   <- min_p_rv
  trait_summary[trait_summary$ResponseVar == rv, "EmpThresh95"] <- threshold_rv
  
  # Optional: print quick info
  message("Response ", rv,
          " | min p = ", signif(min_p_rv, 4),
          " | empirical 95% threshold = ", signif(threshold_rv, 4))
}

############################################
### Combine all results into a single data frame
############################################

all_results <- do.call(rbind, all_results_list)
row.names(all_results) <- NULL

############################################
### Optional: plot for one trait
############################################

# You can pick any trait to visualize
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

rv_to_plot <- response_vars[1]  # change to any trait name you want
plot_data  <- subset(all_results, ResponseVar == rv_to_plot)

ggplot(plot_data, aes(x = seq_along(Pvalue), y = Pvalue,
                      color = Pvalue < 0.05)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                     labels = c(
                       paste("Not significant (", sum(plot_data$Pvalue >= 0.05), ")"),
                       paste("Significant (", sum(plot_data$Pvalue < 0.05), ")")
                     )) +
  labs(title = paste("P-values from", length(common_genotypes),
                     "matched genotypes\nResponse:", rv_to_plot),
       x = "Iteration", y = "P-value",
       color = "Significance") +
  theme_minimal() +
  ylim(0, 1) +
  annotate("text", x = 50, y = 0.9,
           label = paste("Minimum p-value:", round(min(plot_data$Pvalue), 5)))

############################################
### Global summary over all traits
############################################

message("\n=== Global Results Summary ===")
message("Total tests performed (all traits): ", nrow(all_results))
message("Significant at p<0.05: ", sum(all_results$Pvalue < 0.05))
message("Minimum p-value observed: ", min(all_results$Pvalue))
message("Median p-value: ", median(all_results$Pvalue))
