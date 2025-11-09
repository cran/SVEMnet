## ------------------------------------------------------------------
## Comparison: SVEMnet_glm (binomial) vs glmnet (binomial)
## ------------------------------------------------------------------
## Assumes:
##   - library(SVEMnet) is available and contains SVEMnet_glm + predict.svem_binomial
##   - glmnet is installed
## ------------------------------------------------------------------

library(SVEMnet)
library(glmnet)

set.seed(2026)

## ---------- Helper: fast AUC (no extra packages) -------------------
auc_fast <- function(prob, y) {
  # y in {0,1}, prob = predicted probability for class 1
  o <- order(prob, decreasing = TRUE)
  prob <- prob[o]
  y    <- y[o]

  pos <- sum(y == 1)
  neg <- sum(y == 0)
  if (pos == 0 || neg == 0) return(NA_real_)

  # Rank-based AUC
  r <- rank(prob, ties.method = "average")
  (sum(r[y == 1]) - pos * (pos + 1) / 2) / (pos * neg)
}

## ---------- Helper: ROC curve (for plotting) -----------------------
roc_curve <- function(prob, y) {
  o <- order(prob, decreasing = TRUE)
  prob <- prob[o]
  y    <- y[o]

  P <- sum(y == 1)
  N <- sum(y == 0)
  if (P == 0 || N == 0) {
    return(data.frame(fpr = c(0, 1), tpr = c(0, 1)))
  }

  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)

  tpr <- tp / P
  fpr <- fp / N

  # pad with (0,0) and (1,1) for nicer plotting
  data.frame(
    fpr = c(0, fpr, 1),
    tpr = c(0, tpr, 1)
  )
}

## ---------- Simulation settings ------------------------------------
n_train   <- 35
n_test    <- 500
p_signal  <- 3      # number of true predictors
p_noise   <- 7      # pure noise predictors
p_total   <- p_signal + p_noise

B         <- 20     # number of simulation replicates (bump to 50+ if you like)
nBoot_svem <- 60    # bootstrap reps for SVEMnet_glm (moderate for speed)

## Storage for results
auc_svem  <- numeric(B)
auc_glm   <- numeric(B)
acc_svem  <- numeric(B)
acc_glm   <- numeric(B)

## To store last replicate ROC for plotting
roc_last_svem <- NULL
roc_last_glm  <- NULL

## ---------- Main simulation loop -----------------------------------
for (b in seq_len(B)) {
  cat("Replicate", b, "of", B, "...\n")

  ## ---- Generate data ----
  X_tr  <- matrix(rnorm(n_train * p_total), n_train, p_total)
  X_te  <- matrix(rnorm(n_test  * p_total), n_test,  p_total)

  colnames(X_tr) <- paste0("x", seq_len(p_total))
  colnames(X_te) <- colnames(X_tr)

  # True linear predictor uses only first p_signal
  beta_true <- c(0.8, -1.0, 0.6, rep(0, p_noise))  # length p_total
  intercept <- -0.4

  lp_tr <- intercept + X_tr %*% beta_true
  lp_te <- intercept + X_te %*% beta_true

  p_tr  <- plogis(lp_tr)
  p_te  <- plogis(lp_te)

  y_tr  <- rbinom(n_train, size = 1, prob = p_tr)
  y_te  <- rbinom(n_test,  size = 1, prob = p_te)

  dat_tr <- data.frame(y = y_tr, X_tr)
  dat_te <- data.frame(y = y_te, X_te)

  ## ---- Fit SVEMnet_glm (binomial) ----
  fit_svem <- SVEMnet_glm(
    y ~ .,
    data    = dat_tr,
    family  = "binomial",
    nBoot   = nBoot_svem,
    relaxed = TRUE
  )

  # SVEM probabilities on test data (ensemble mean)
  p_hat_svem <- as.numeric(
    predict(
      fit_svem,
      newdata = dat_te,
      type    = "response",
      agg     = "mean"
    )
  )

  ## ---- Fit glmnet (binomial) with CV ----
  X_tr_mm <- model.matrix(y ~ ., dat_tr)[, -1, drop = FALSE]
  X_te_mm <- model.matrix(y ~ ., dat_te)[, -1, drop = FALSE]

  cv_glm <- cv.glmnet(
    x            = X_tr_mm,
    y            = y_tr,
    family       = "binomial",
    alpha        = 1,
    nfolds       = 5,
    type.measure = "auc"
  )

  p_hat_glm <- as.numeric(
    predict(
      cv_glm,
      newx = X_te_mm,
      s    = "lambda.min",
      type = "response"
    )
  )

  ## ---- Metrics on test set ----
  # AUC
  auc_svem[b] <- auc_fast(p_hat_svem, y_te)
  auc_glm[b]  <- auc_fast(p_hat_glm,  y_te)

  # Accuracy at threshold 0.5
  y_pred_svem <- as.numeric(p_hat_svem >= 0.5)
  y_pred_glm  <- as.numeric(p_hat_glm  >= 0.5)

  acc_svem[b] <- mean(y_pred_svem == y_te)
  acc_glm[b]  <- mean(y_pred_glm  == y_te)

  ## ---- Store ROC for last replicate ----
  if (b == B) {
    roc_last_svem <- roc_curve(p_hat_svem, y_te)
    roc_last_glm  <- roc_curve(p_hat_glm,  y_te)
  }
}

cat("\n------------------------------------------------------------\n")
cat("Summary over", B, "replicates\n")
cat("------------------------------------------------------------\n\n")

## ---------- Summaries ---------------------------------------------
cat("AUC (SVEMnet_glm):   mean =", mean(auc_svem, na.rm = TRUE),
    ", sd =", sd(auc_svem, na.rm = TRUE), "\n")
cat("AUC (glmnet):        mean =", mean(auc_glm,  na.rm = TRUE),
    ", sd =", sd(auc_glm,  na.rm = TRUE), "\n\n")

cat("Accuracy (SVEMnet_glm): mean =", mean(acc_svem, na.rm = TRUE),
    ", sd =", sd(acc_svem, na.rm = TRUE), "\n")
cat("Accuracy (glmnet):      mean =", mean(acc_glm,  na.rm = TRUE),
    ", sd =", sd(acc_glm,  na.rm = TRUE), "\n\n")

cat("In how many replicates is SVEMnet_glm AUC > glmnet AUC?\n")
cat("   count =", sum(auc_svem > auc_glm, na.rm = TRUE),
    "out of", B, "\n")
cat("   proportion =", mean(auc_svem > auc_glm, na.rm = TRUE), "\n\n")

cat("In how many replicates is SVEMnet_glm accuracy > glmnet accuracy?\n")
cat("   count =", sum(acc_svem > acc_glm, na.rm = TRUE),
    "out of", B, "\n")
cat("   proportion =", mean(acc_svem > acc_glm, na.rm = TRUE), "\n\n")

## ---------- ROC curve for last replicate ---------------------------
op <- par(no.readonly = TRUE)
par(mar = c(5, 5, 4, 2) + 0.1)

plot(roc_last_glm$fpr, roc_last_glm$tpr,
     type = "l", lwd = 2,
     xlab = "False Positive Rate",
     ylab = "True Positive Rate",
     main = "ROC curve (last replicate)\nSVEMnet_glm vs glmnet")
lines(roc_last_svem$fpr, roc_last_svem$tpr,
      lwd = 2, col = 2)
abline(0, 1, lty = 3)

legend("bottomright",
       legend = c(
         paste0("glmnet (AUC = ", round(auc_glm[B],  3), ")"),
         paste0("SVEMnet_glm (AUC = ", round(auc_svem[B], 3), ")")
       ),
       col = c(1, 2), lwd = 2, bty = "n")

par(op)

## Also print the AUCs for the last replicate:
cat("Last replicate AUCs:\n")
cat("  glmnet:       ", round(auc_glm[B],  4), "\n")
cat("  SVEMnet_glm:  ", round(auc_svem[B], 4), "\n")


## Paired deltas
d_auc <- auc_svem - auc_glm
d_acc <- acc_svem - acc_glm

## Summary
cat("AUC delta (SVEM - glmnet): mean =", mean(d_auc),
    "sd =", sd(d_auc), "median =", median(d_auc), "\n")
cat("ACC delta (SVEM - glmnet): mean =", mean(d_acc),
    "sd =", sd(d_acc), "median =", median(d_acc), "\n\n")

## Paired t-tests
cat("Paired t-test (AUC):\n"); print(t.test(auc_svem, auc_glm, paired = TRUE))
cat("\nPaired t-test (Accuracy):\n"); print(t.test(acc_svem, acc_glm, paired = TRUE))

## Wilcoxon signed-rank (robust)
cat("\nWilcoxon (AUC):\n"); print(wilcox.test(auc_svem, auc_glm, paired = TRUE))
cat("\nWilcoxon (Accuracy):\n"); print(wilcox.test(acc_svem, acc_glm, paired = TRUE))

## 95% bootstrap CIs for mean delta
boot_mean_ci <- function(x, R = 2000) {
  m <- replicate(R, mean(sample(x, replace = TRUE)))
  c(mean = mean(x), lwr = quantile(m, 0.025), upr = quantile(m, 0.975))
}
cat("\nBootstrap 95% CI for mean AUC delta:\n"); print(boot_mean_ci(d_auc))
cat("\nBootstrap 95% CI for mean ACC delta:\n"); print(boot_mean_ci(d_acc))

avg_precision <- function(prob, y) {
  o <- order(prob, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  P  <- sum(y == 1); N <- sum(y == 0)
  if (P == 0 || N == 0) return(NA_real_)
  precision <- tp / (tp + fp)
  recall    <- tp / P
  # step-wise AP via trapezoids on recall (monotone increasing)
  sum(diff(c(0, recall)) * precision)
}

## You already have p_hat_svem, p_hat_glm, and y_te in the loop for the last replicate.
## If you don't, recreate them for the last replicate similarly to before.

ap_svem <- avg_precision(p_hat_svem, y_te)
ap_glm  <- avg_precision(p_hat_glm,  y_te)
cat("\nAverage precision (last replicate):\n")
cat("  SVEMnet_glm:", round(ap_svem, 3), "\n")
cat("  glmnet     :", round(ap_glm,  3), "\n")

## Plot PR curve for last replicate
pr_curve <- function(prob, y) {
  o <- order(prob, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  P  <- sum(y == 1); N <- sum(y == 0)
  precision <- tp / pmax(1, tp + fp)
  recall    <- if (P == 0) rep(0, length(tp)) else tp / P
  data.frame(recall = c(0, recall), precision = c(1, precision))
}

pr_svem <- pr_curve(p_hat_svem, y_te)
pr_glm  <- pr_curve(p_hat_glm,  y_te)

op <- par(no.readonly = TRUE)
par(mar = c(5,5,4,2)+0.1)
plot(pr_glm$recall, pr_glm$precision, type = "l", lwd = 2,
     xlab = "Recall", ylab = "Precision",
     main = "Precision–Recall (last replicate)")
lines(pr_svem$recall, pr_svem$precision, lwd = 2, col = 2)
legend("topright",
       legend = c(paste0("glmnet (AP=", round(ap_glm,3), ")"),
                  paste0("SVEMnet_glm (AP=", round(ap_svem,3), ")")),
       col = c(1,2), lwd = 2, bty = "n")
par(op)


best_thresh <- function(prob, y, metric = c("youden","f1")) {
  metric <- match.arg(metric)
  th <- sort(unique(prob))
  if (length(th) > 500) th <- quantile(prob, probs = seq(0,1,length.out=500))
  score <- function(t) {
    pred <- as.integer(prob >= t)
    tp <- sum(pred==1 & y==1); tn <- sum(pred==0 & y==0)
    fp <- sum(pred==1 & y==0); fn <- sum(pred==0 & y==1)
    if (metric=="youden") {
      tpr <- ifelse(tp+fn>0, tp/(tp+fn), 0)
      fpr <- ifelse(fp+tn>0, fp/(fp+tn), 0)
      return(tpr - fpr)
    } else {
      prec <- ifelse(tp+fp>0, tp/(tp+fp), 0)
      rec  <- ifelse(tp+fn>0, tp/(tp+fn), 0)
      if (prec+rec==0) return(0)
      return(2*prec*rec/(prec+rec))
    }
  }
  s <- vapply(th, score, numeric(1))
  th[which.max(s)]
}

t_svem_youden <- best_thresh(p_hat_svem, y_te, "youden")
t_glm_youden  <- best_thresh(p_hat_glm,  y_te, "youden")

t_svem_f1 <- best_thresh(p_hat_svem, y_te, "f1")
t_glm_f1  <- best_thresh(p_hat_glm,  y_te, "f1")

cat("\nBest thresholds (last replicate):\n")
cat("  SVEM (Youden):", round(t_svem_youden,3),
    "  glmnet (Youden):", round(t_glm_youden,3), "\n")
cat("  SVEM (F1):    ", round(t_svem_f1,3),
    "  glmnet (F1):    ", round(t_glm_f1,3), "\n")








