################################################################################
#                          U-statistic cross-fitting                           #
################################################################################
#                                                                              #
#         This script contains functions to perform the U-statistic            #
#         cross-fitting procedure for the estimation                           #
#                                                                              #
################################################################################
#                                                                              #
#                               Dependencies                                   #
#                                                                              #
# No scripts are needed for this script                                        #
#                                                                              #
#                                                                              #
################################################################################


# Escanciano-Terschuur Algorithm 1: K triangles and four subrectangles for
# each pair of groups; L = K(2K-1). For n divisible by 2K all training sets
# have size n(K-1)/K. Otherwise balanced sizes differ by at most two.
# The random partition depends only on row indices.
make_pair_folds <- function(n, K = 2L, seed = NULL) {
  stopifnot(length(K) == 1L, 
            K == as.integer(K), 
            K >= 2, 
            n >= 2 * K,
            n == as.integer(n))
  
  #browser()
  set.seed(seed)
  
  #-----------------------------------------------------------------------------
  # 0. Permute the data
  #-----------------------------------------------------------------------------
  perm <- sample.int(n)
  
  #-----------------------------------------------------------------------------
  # 1. Partition the n observations into K parts G_1,...,G_K
  #-----------------------------------------------------------------------------
  groups <- unname(split(perm, rep(seq_len(K), length.out = n)))
  
  
  #-----------------------------------------------------------------------------
  # 2. For each part, construct blocks of pairs (Within-group blocks)
  #-----------------------------------------------------------------------------
  blocks <- lapply(groups, function(g) list(type = "triangle", left = g))
  
  
  #-----------------------------------------------------------------------------
  # 3. Partition into two equal disjoint subsets 
  #-----------------------------------------------------------------------------
  halves <- lapply(groups, function(g) {
    cut <- ceiling(length(g) / 2)
    list(g[seq_len(cut)], g[seq.int(cut + 1L, length(g))])
  })
  
  #-----------------------------------------------------------------------------
  # 4. Construct blocks of pairs for each pair of parts (Between-group blocks)
  #-----------------------------------------------------------------------------
  for (u in seq_len(K - 1L)) for (v in seq.int(u + 1L, K))
    for (a in 1:2) for (b in 1:2)
      blocks[[length(blocks) + 1L]] <- list(type = "rectangle",
                                            left = halves[[u]][[a]], 
                                            right = halves[[v]][[b]])
  
  #-----------------------------------------------------------------------------
  # 5. Compute the evaluation and training sets for each block
  #-----------------------------------------------------------------------------
  blocks <- lapply(blocks, function(b) {
    b$eval <- c(b$left, b$right)
    b$train <- setdiff(seq_len(n), b$eval)
    b$n_pairs <- if (b$type == "triangle") choose(length(b$left), 2)
    else length(b$left) * length(b$right)
    b
  })
  structure(list(n = n, K = K, groups = groups, blocks = blocks),
            class = "pair_folds")
}




