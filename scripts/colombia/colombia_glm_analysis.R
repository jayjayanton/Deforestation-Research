#!/usr/bin/env Rscript
# Colombia deforestation GLM analysis — portable research workflow
# Based on my historical WWF GLM script. Requires an authorized input CSV.
# Usage: Rscript scripts/colombia/colombia_glm_analysis.R --input data/colombia.csv --output-dir outputs/colombia
#
# The model uses binomial logistic regression, not Random Forest. The historical
# Random Forest risk map in this portfolio was produced by a separate analysis.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  pos <- match(flag, args)
  if (is.na(pos)) return(default)
  if (pos == length(args)) stop("Missing value after ", flag)
  args[[pos + 1L]]
}
input <- get_arg("--input")
out_dir <- get_arg("--output-dir", "outputs/colombia")
target <- get_arg("--target", "Deforestation")
predictor_arg <- get_arg("--predictors")
seed <- as.integer(get_arg("--seed", "42"))
if (is.null(input) || !file.exists(input)) {
  stop("Supply an existing CSV with --input. Example: --input data/colombia.csv")
}
if (!is.finite(seed)) stop("--seed must be an integer")
set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
df <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
if (!target %in% names(df)) stop("Missing target column: ", target)
if (!"Elevation" %in% names(df)) stop("Missing Elevation column for the historical baseline")

# Explicit predictor selection avoids fragile column-number exclusions in the
# original script. Supply comma-separated predictor names to --predictors.
id_cols <- c("ID", "id_cell", "OBJECTID", "Lattitude", "Latitude", "Longitude",
             "Clust_150km", "Col_ac")
if (is.null(predictor_arg)) {
  predictors <- setdiff(names(df), c(target, id_cols))
  predictors <- predictors[vapply(df[predictors], is.numeric, logical(1))]
} else {
  predictors <- trimws(strsplit(predictor_arg, ",", fixed = TRUE)[[1]])
}
if (anyDuplicated(predictors) || !length(predictors)) stop("No valid predictors selected")
if (!all(predictors %in% names(df))) {
  stop("Missing predictors: ", paste(setdiff(predictors, names(df)), collapse = ", "))
}
if (!"Elevation" %in% predictors) predictors <- c(predictors, "Elevation")
if (!all(vapply(df[predictors], is.numeric, logical(1)))) {
  stop("All selected predictors must be numeric")
}
# Historical distance variables were converted from metres to kilometres.
# This is optional because other datasets may already use kilometres.
if (identical(get_arg("--distances-in-metres", "false"), "true")) {
  distance_cols <- grep("_Dist$", predictors, value = TRUE)
  df[distance_cols] <- lapply(df[distance_cols], function(x) x / 1000)
}
# Require documented 0/1 labels. Do not silently infer factor ordering.
raw_y <- df[[target]]
if (is.factor(raw_y)) raw_y <- as.character(raw_y)
y <- suppressWarnings(as.numeric(raw_y))
if (anyNA(y) || !all(unique(y) %in% c(0, 1)) || length(unique(y)) != 2) {
  stop("The target must be explicitly coded as 0 and 1, with both classes present")
}
df[[target]] <- y
needed <- unique(c(target, predictors))
complete <- complete.cases(df[, needed, drop = FALSE])
message("Excluding ", sum(!complete), " rows with missing target/predictor values")
df <- df[complete, , drop = FALSE]
if (nrow(df) < 30) stop("At least 30 complete records are required")
if (min(table(df[[target]])) < 5) stop("At least five observations per class are required")

formula_full <- reformulate(predictors, response = target)
formula_naive <- reformulate("Elevation", response = target)
fit_safe <- function(formula, data) {
  suppressWarnings(glm(formula, data = data, family = binomial(),
                       control = glm.control(maxit = 100)))
}
full <- fit_safe(formula_full, df)
stepwise <- suppressWarnings(step(full, direction = "both", trace = 0))
naive <- fit_safe(formula_naive, df)

models <- list(Full = full, Stepwise = stepwise, Naive = naive)
coefs <- do.call(rbind, lapply(names(models), function(name) {
  sm <- summary(models[[name]])$coefficients
  data.frame(model = name, term = rownames(sm), estimate = sm[, 1],
             std_error = sm[, 2], z = sm[, 3], p_value = sm[, 4],
             row.names = NULL)
}))
write.csv(coefs, file.path(out_dir, "glm_coefficients.csv"), row.names = FALSE)

# Holdout split; fit feature selection on training data only.
idx0 <- which(df[[target]] == 0)
idx1 <- which(df[[target]] == 1)
train_idx <- sort(c(sample(idx0, floor(length(idx0) * .7)),
                    sample(idx1, floor(length(idx1) * .7))))
train <- df[train_idx, , drop = FALSE]
test <- df[-train_idx, , drop = FALSE]
train_full <- fit_safe(formula_full, train)
train_step <- suppressWarnings(step(train_full, direction = "both", trace = 0))
train_naive <- fit_safe(formula_naive, train)
train_models <- list(Full = train_full, Stepwise = train_step, Naive = train_naive)
pred <- data.frame(observed = test[[target]])
for (name in names(train_models)) {
  pred[[name]] <- predict(train_models[[name]], newdata = test, type = "response")
}
write.csv(pred, file.path(out_dir, "holdout_predictions.csv"), row.names = FALSE)

# Evaluation at a fixed 0.5 threshold; do not optimize using the test set.
score <- function(actual, probability, threshold = .5) {
  if (anyNA(probability)) return(c(accuracy = NA, sensitivity = NA,
                                  specificity = NA, tss = NA, auc = NA))
  predicted <- as.integer(probability >= threshold)
  tp <- sum(predicted == 1 & actual == 1)
  tn <- sum(predicted == 0 & actual == 0)
  fp <- sum(predicted == 1 & actual == 0)
  fn <- sum(predicted == 0 & actual == 1)
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
  # Rank-based AUC with tie handling, requires both observed classes.
  auc <- if (length(unique(actual)) == 2) {
    r <- rank(probability, ties.method = "average")
    (sum(r[actual == 1]) - sum(actual == 1) * (sum(actual == 1) + 1) / 2) /
      (sum(actual == 1) * sum(actual == 0))
  } else NA_real_
  c(accuracy = mean(predicted == actual), sensitivity = sens,
    specificity = spec, tss = sens + spec - 1, auc = auc)
}
metrics <- do.call(rbind, lapply(names(train_models), function(name) {
  data.frame(model = name, t(score(pred$observed, pred[[name]])),
             row.names = NULL)
}))
write.csv(metrics, file.path(out_dir, "holdout_metrics.csv"), row.names = FALSE)

# Optional spatial holdout, using geographic coordinates to form broad folds.
# This is a simplified geographic-group CV, not a substitute for the original
# 150-km spatial clustering or a spatial autocorrelation assessment.
lat_col <- if ("Lattitude" %in% names(df)) "Lattitude" else
           if ("Latitude" %in% names(df)) "Latitude" else NULL
if (!is.null(lat_col) && "Longitude" %in% names(df)) {
  coords <- df[, c(lat_col, "Longitude"), drop = FALSE]
  if (all(complete.cases(coords)) && all(vapply(coords, is.numeric, logical(1)))) {
    lat0 <- mean(coords[[lat_col]])
    xy <- cbind(coords$Longitude * 111.32 * cos(lat0 * pi / 180),
                coords[[lat_col]] * 111.32)
    n_folds <- min(5L, nrow(df) %/% 20L)
    if (n_folds >= 2L) {
      clusters <- kmeans(xy, centers = n_folds, nstart = 20)$cluster
      spatial_rows <- list()
      for (fold in sort(unique(clusters))) {
        tr <- df[clusters != fold, , drop = FALSE]
        te <- df[clusters == fold, , drop = FALSE]
        if (length(unique(tr[[target]])) < 2) next
        fitted <- tryCatch(fit_safe(formula_full, tr), error = function(e) NULL)
        if (is.null(fitted)) next
        probs <- tryCatch(predict(fitted, newdata = te, type = "response"),
                          error = function(e) rep(NA_real_, nrow(te)))
        spatial_rows[[length(spatial_rows) + 1L]] <-
          data.frame(fold = fold, n = nrow(te),
                     t(score(te[[target]], probs)))
      }
      if (length(spatial_rows)) {
        write.csv(do.call(rbind, spatial_rows),
                  file.path(out_dir, "spatial_fold_metrics.csv"), row.names = FALSE)
      }
    }
  }
}
message("Wrote coefficients and validation outputs to ", normalizePath(out_dir))
