#' Cross-fitted U-statistic for the embedded centred CATE
#'
#' Source u_statistic_crossfit.R and nuisance_func.R before calling this
#' function. Each outer block trains only on block$train. Within that sample,
#' nuisance_func() fits the regressions and marginal on disjoint inner folds.
#'
#' @param data Data frame with Y, binary A, and numeric covariates.
#' @param pair_folds The object returned by make_pair_folds() for these data,
#'   in the same row order.
#' @param ... Arguments to nuisance_func(), e.g. splits, K, SL.lib, randomized,
#'   pi_hat, hal, or kernel. K here controls the learner's internal CV.
#' @param keep_pair_scores Save the symmetric n-by-n matrix of cross-fitted
#'   pair scores with zero diagonal for bootstrap_embedded_u(). Default FALSE.
#' @return A list with estimate, n_pairs, and a blocks data frame. block_sum
#'   is the inner-fold-weighted sum of the scores for the block. The estimate
#'   is sum(block_sum) / choose(nrow(data), 2). When requested, pair_scores[i,j]
#'   holds the complete inner-averaged score for original observations i and j.
#' @export
#' @examples
#' \dontrun{
#' set.seed(3)
#' n <- 125L
#' dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), A = rbinom(n, 1, 0.5))
#' dat$Y <- with(dat, x2 + A * (1 + x1)) + rnorm(n)
#' plan <- make_pair_folds(n, K = 2L, seed = 17L)
#' set.seed(18)
#' ans <- crossfit_embedded_u(dat, plan, splits = 2L, K = 3L,
#'                           SL.lib = "SL.glm", randomized = TRUE, pi_hat = 0.5,
#'                           keep_pair_scores = TRUE)
#' ans$estimate
#' stopifnot(sum(ans$blocks$n_pairs) == choose(n, 2),
#'           isTRUE(all.equal(sum(ans$pair_scores[upper.tri(ans$pair_scores)]) /
#'                            ans$n_pairs, ans$estimate)),
#'           isTRUE(all.equal(ans$estimate,
#'                            weighted.mean(ans$blocks$block_mean, ans$blocks$n_pairs))))
#' set.seed(19)
#' test <- bootstrap_embedded_u(ans, B = 1999L, alpha = 0.05)
#' test[c("statistic", "critical_value", "p_value", "reject")]
#' # Reproduce the first Gaussian draw directly from the unordered pairs.
#' set.seed(19)
#' xi <- rnorm(n)
#' upper <- upper.tri(ans$pair_scores)
#' direct <- 2 / (n - 1) * sum((ans$pair_scores * tcrossprod(xi))[upper])
#' stopifnot(isTRUE(all.equal(direct, test$draws[1])))
#' # For this small example, enumerate the unordered pair labels once.
#' pairs <- do.call(rbind, lapply(plan$blocks, function(b) {
#'   if (b$type == "triangle") return(t(utils::combn(b$left, 2L)))
#'   unname(as.matrix(expand.grid(b$left, b$right)))
#' }))
#' keys <- paste(pmin(pairs[, 1], pairs[, 2]), pmax(pairs[, 1], pairs[, 2]), sep = "/")
#' stopifnot(length(keys) == choose(n, 2), all(pairs[, 1] != pairs[, 2]),
#'           anyDuplicated(keys) == 0L)
#' }
crossfit_embedded_u <- function(data, pair_folds, ..., keep_pair_scores = FALSE) {
  
  
  #-----------------------------------------------------------------------------
  # 1.0 Verify the inputs
  #-----------------------------------------------------------------------------
  
  if (!is.data.frame(data) || !all(c("Y", "A") %in% names(data))) {
    stop("data must be a data frame with Y, A, and the covariates.")
  }
  
  data <- as.data.frame(data)
  n <- nrow(data)
  
  if (!is.logical(keep_pair_scores) || length(keep_pair_scores) != 1L ||
      is.na(keep_pair_scores)) {
    stop("keep_pair_scores must be TRUE or FALSE.")
  }
  
  if (!inherits(pair_folds, "pair_folds") || !isTRUE(pair_folds$n == n)) {
    stop("Supply the result of make_pair_folds() for the rows of data.")
  }
  
  blocks <- pair_folds$blocks
  n_pairs <- choose(n, 2L)
  
  counts <- vapply(blocks, .eu_block_count, numeric(1), rows = seq_len(n))
  if (!length(blocks) || sum(counts) != n_pairs) {
    stop("The pair counts in the supplied blocks do not sum to choose(n, 2).")
  }
  
  block_sums <- numeric(length(blocks))
  if (keep_pair_scores) {
    H <- matrix(NA_real_, nrow = n, ncol = n)
    diag(H) <- 0
  }
  
  
  #-----------------------------------------------------------------------------
  # 2.0 Loop over the outer folds and average the inner-fold scores for each block. 
  # The outer block trains on the complement of its evaluation rows. 
  # The inner folds train on disjoint subsets of that training complement.
  #-----------------------------------------------------------------------------
  
  for (l in seq_along(blocks)) {
    
    b <- blocks[[l]]
    
    contribution <- tryCatch({
      
      # Train the nuisance functions on the training complement of this outer block.
      fit <- nuisance_func(train_data = data[b$train, , drop = FALSE], ...)
      
      # Local positions in z are the left rows followed by the right rows.
      z <- data[c(b$left, b$right), , drop = FALSE]
      x <- z[, fit$x_names, drop = FALSE]
      
      # Average the inner-fold kernels for this outer block
      fit$average(function(s) .eu_block_sum(s, z, x, b, keep_pair_scores))
      
    }, error = function(err) {
      stop(sprintf("Outer block %d (%s): %s", l, b$type, conditionMessage(err)),
           call. = FALSE)
    })
    
    block_sums[l] <- sum(contribution)
    if (keep_pair_scores) {
      
      if (identical(b$type, "triangle")) {
        
        H[b$left, b$left] <- contribution + t(contribution)
        
      } else {
        
        H[b$left, b$right] <- contribution
        H[b$right, b$left] <- t(contribution)
        
      }
    }
  }
  
  
  by_block <- data.frame(
    block = seq_along(blocks),
    type = vapply(blocks, function(b) b$type, character(1)),
    n_train = vapply(blocks, function(b) length(b$train), integer(1)),
    n_eval = vapply(blocks, function(b) length(b$eval), integer(1)),
    n_pairs = counts, block_sum = block_sums, block_mean = block_sums / counts
  )
  
  result <- list(estimate = sum(block_sums) / n_pairs, n_pairs = n_pairs, blocks = by_block)
  
  if (keep_pair_scores) {
    if (anyNA(H)) stop("The supplied pair plan did not cover every unordered pair once.")
    result$pair_scores <- H
  }
  
  result
  
}


# For one outer block and one inner fit, return the sum or the individual
# scores. Triangles store only the strict upper triangle in local row order.
.eu_block_sum <- function(s, z, x, b, keep_pair_scores = FALSE) {
  
  # Retrieve the R term from the inner fit.
  r <- s$residual(a = z$A, y = z$Y, x = x)
  
  # Calculate the first-order influence function term
  h <- (s$embedded_cate(x) - s$mean_embedding) * (r + s$centred_cate(x))
  
  # 
  triangle <- identical(b$type, "triangle")
  left <- seq_along(b$left)
  right <- if (triangle) left else length(left) + seq_along(b$right)
  
  # Restore the matrix shape when either side of a rectangle has one row.
  kc <- s$centred_kernel(x[left, , drop = FALSE], x[right, , drop = FALSE])
  if (length(kc) != as.double(length(left)) * length(right)) {
    stop("centred_kernel returned the wrong number of pair values.")
  }
  kc <- matrix(kc, nrow = length(left), ncol = length(right))
  
  
  # In a triangle use i < j in the local group order. In a rectangle every
  # left/right combination is one unordered pair, regardless of row labels.
  sums <- numeric(length(right))
  for (col in seq_along(right)) {
    i <- if (triangle) seq_len(col - 1L) else left
    j <- right[col]
    
    # Calculate the pairwise kernels
    values <- h[i] + h[j] + r[i] * r[j] * kc[i, col] - s$embedding_norm_sq
    if (any(!is.finite(values))) stop("A pair score is non-finite.")
    if (keep_pair_scores) {
      # Reuse the kernel column's storage. Unused triangle entries stay zero.
      kc[, col] <- 0
      kc[i, col] <- values
    } else {
      sums[col] <- sum(values)
    }
  }
  if (keep_pair_scores) kc else sum(sums)
}





#' I.i.d. multiplier bootstrap for the constant-CATE null
#'
#' For each replicate draw n independent, mean-zero, unit-variance multipliers,
#' one for each original observation. Given the saved pair scores H, the draw
#' is xi' H xi / (n-1), on the scale of n times the observed estimate.
#' All folds, fitted models, kernel parameters, and pair scores are held fixed.
#'
#' Under the homogeneity null the oracle pair kernel is degenerate. A
#' sufficient justification with estimated scores is (i) oracle equivalence
#' of the observed statistic on the n scale and (ii) vanishing empirical
#' mean-square error of the off-diagonal scores relative to the oracle kernel.
#' Reference: Dehling and Mikosch (1994), equation (1.6) and Theorem 3.1.
#'
#' @param object Result of crossfit_embedded_u(..., keep_pair_scores = TRUE).
#' @param B Positive number of Monte Carlo replicates.
#' @param alpha Upper-tail test level between zero and one.
#' @param multiplier "normal" uses standard Gaussian weights. "rademacher"
#'   uses independent equiprobable signs. Use set.seed() before calling.
#' @return A list with statistic (n * estimate), draws, critical_value, the
#'   upper-tail p_value with the plus-one correction, reject, alpha, B and
#'   multiplier. The order-statistic cutoff agrees with p_value <= alpha,
#'   with strict exceedance required. It is Inf if B is too small to reject.
#' @export
bootstrap_embedded_u <- function(object, B = 1999L, alpha = 0.05,
                                 multiplier = c("normal", "rademacher")) {
  
  H <- object$pair_scores
  
  if (is.null(H)) {
    stop("Refit with crossfit_embedded_u(..., keep_pair_scores = TRUE) to save individual pair scores.")
  }
  
  if (!is.matrix(H) || !is.numeric(H) || nrow(H) < 2L || nrow(H) != ncol(H) ||
      any(!is.finite(H)) || any(diag(H) != 0) || !isSymmetric(H)) {
    stop("pair_scores must be a finite symmetric matrix with zero diagonal.")
  }
  
  if (!is.numeric(object$estimate) || length(object$estimate) != 1L ||
      !is.finite(object$estimate)) stop("The fitted estimate must be a finite number.")
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B) || B > .Machine$integer.max) {
    stop("B must be a positive integer.")
  }
  
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) stop("alpha must be strictly between zero and one.")
  multiplier <- match.arg(multiplier)
  
  
  n <- nrow(H)
  draws <- numeric(B)
  
  # Each column of xi is a full replicate. Batching limits working memory.
  for (first in seq.int(1, B, by = 64L)) {
    cols <- seq.int(first, min(first + 63, B))
    size <- as.double(n) * length(cols)
    weights <- if (multiplier == "normal") {
      stats::rnorm(size)
    } else {
      2 * sample.int(2L, size, replace = TRUE) - 3
    }
    xi <- matrix(weights, nrow = n, ncol = length(cols))
    draws[cols] <- colSums(xi * (H %*% xi)) / (n - 1)
  }
  statistic <- n * object$estimate
  p_value <- (1 + sum(draws >= statistic)) / (B + 1)
  
  
  # Use the same floating-point p-value grid for the cutoff and the decision.
  rank <- B + 1 - findInterval(alpha, seq_len(B) / (B + 1))
  critical_value <- if (rank > B) Inf else sort(draws, partial = rank)[rank]
  
  list(statistic = statistic, draws = draws, critical_value = critical_value,
       p_value = p_value, reject = p_value <= alpha, alpha = alpha, B = B,
       multiplier = multiplier)
}




# Verify training complements and pair counts before running any learners.
.eu_block_count <- function(b, rows) {
  triangle <- identical(b$type, "triangle")
  if (!triangle && !identical(b$type, "rectangle")) stop("Unknown pair block type.")
  if ((triangle && (length(b$left) < 2L || length(b$right) != 0L)) ||
      (!triangle && (length(b$left) == 0L || length(b$right) == 0L))) {
    stop("An outer block has invalid sides.")
  }
  held_out <- c(b$left, b$right)
  if (!is.numeric(held_out) || anyDuplicated(held_out) || !all(held_out %in% rows) ||
      length(b$eval) != length(held_out) || !setequal(b$eval, held_out) ||
      anyDuplicated(b$train) || !is.numeric(b$train) ||
      !setequal(b$train, setdiff(rows, held_out))) {
    stop("Each block must train on exactly the complement of its evaluation rows.")
  }
  count <- if (triangle) {
    choose(length(b$left), 2L)
  } else {
    as.double(length(b$left)) * length(b$right)
  }
  if (!isTRUE(b$n_pairs == count)) {
    stop("The pair count of an outer block is inconsistent with its sides.")
  }
  count
}
