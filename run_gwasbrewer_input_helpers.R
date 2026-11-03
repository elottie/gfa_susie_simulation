#source('gwasbrewer_input_helpers.R')
source('build_gwasbrewer_input.R')
options(warn=1)

# 1. Load credible sets and their relationships
ref <- readRDS("/nfs/turbo/sph-jvmorr/GFA_susie/simulation/metab_susie/ld_ref/metsim_ld_ref_1M.RDS")
pool <- load_cs_pool(
  dir = "/nfs/turbo/sph-jvmorr/GFA_susie/simulation/metab_susie/credible_set_info/cs_objects/",
  prefix = "metsim_2022.susie.all_in.pval5e-8",
  chrs = 1:22,
  r2_min = 0.05,
  ref_info = ref$snp_info
)

#print(str(pool))

# For each credible set, extract the trait name (should only be one, why we have cs_key to be diff for diff traits w/ same cs)
traits <- vapply(
  pool$cs,
  function(x) {
    value <- unique(as.character(x$metabolite))
    value <- value[!is.na(value)]

    if (length(value) != 1) {
      stop("Each credible set must have exactly one non-missing metabolite")
    }

    value
  },
  character(1)
)

# Count credible sets by trait
n_cs_by_trait <- setNames(
  as.integer(table(traits)),
  names(table(traits))
)

print(head(n_cs_by_trait))

# Deterministic split
# we assume rows of myF (traits) in same order as rows of n_cs_by_trait, can change later
assign_cs_to_factors <- function(
    myF,
    cs_per_trait,
    factor_fraction = 0.5,
    detection_buffer = 1.5
) {
  K <- ncol(myF)
  T <- nrow(myF)

  if (anyNA(cs_per_trait)) {
    stop("cs_per_trait must contain non-negative integer counts")
  }

  if (is.na(factor_fraction)) {
    stop("factor_fraction must be one number between 0 and 1")
  }

  # multiply each element of cs_by_trait by detection buffer, simulating we did not catch all credible sets
  cs_per_trait <- floor(cs_per_trait * detection_buffer + 0.5)

  # number of credible sets for the trait, length of traits
  # Round-half-up rather than R's round-to-even
  n_cs_factor <- as.integer(floor(cs_per_trait * factor_fraction + 0.5))
  n_cs_direct <- cs_per_trait - n_cs_factor

  # Rows = traits; columns = factors
  n_cs_mat <- matrix(0, nrow = T, ncol = K)

  for (i in seq_len(T)) {
    # Allocate the trait cs destined for factors evenly across factors
    n_cs_mat[i, ] <- n_cs_factor[i] %/% K

    # Assign leftover factor-destined cs to the first factor
    n_cs_factor_remain <- n_cs_factor[i] %% K
    if (n_cs_factor_remain > 0) {
      n_cs_mat[i, seq_len(n_cs_factor_remain)] <- n_cs_mat[i, seq_len(n_cs_factor_remain)] + 1
    }
  }

  list(
    n_cs_factor = n_cs_factor,
    n_cs_direct = n_cs_direct,
    n_cs_mat = n_cs_mat
  )
}

myF <- matrix(
  c(
    0, -1,
    1,  1,
    2,  1
  ),
  nrow = 3,
  byrow = TRUE
)

myF <- apply(
  myF,
  2,
  function(x) x / sqrt(sum(x^2))
)

cs_per_trait <- n_cs_by_trait[1:nrow(myF)]

factor_spec <- assign_cs_to_factors(
  myF = myF,
  cs_per_trait = cs_per_trait,
  factor_fraction = 0.5
)

factor_spec

stop('told to stop here')

#factor_fraction <- 0.3

#n_factor <- floor(n_cs_by_trait * factor_fraction + 0.5)
#n_direct <- n_cs_by_trait - n_factor

#print(n_factor)
#print(n_direct)

#count_summary <- data.table(
#  trait = names(n_cs_by_trait),
#  n_total = as.integer(n_cs_by_trait),
#  n_factor = as.integer(n_factor),
#  n_direct = as.integer(n_direct)
#)

#count_summary

# Allocate direct effects evenly across traits
T <- nrow(myF)

n_causal_direct <- rep.int(
  factor_spec$n_direct_total %/% T,
  T
)

remainder <- factor_spec$n_direct_total %% T

if (remainder > 0L) {
  n_causal_direct[seq_len(remainder)] <-
    n_causal_direct[seq_len(remainder)] + 1L
}

# Factor-mediated inputs
factor_inputs <- make_gwasbrewer_inputs(
  pool = pool,
  n_col = factor_spec$n_col,
  n_causal_fn = factor_spec$n_causal_fn,
  n_close = 2,
  seed = 123
)

identical(factor_inputs$selected$pool_cs_name,factor_inputs$selected$cs_key)

# Direct-effect inputs
direct_inputs <- make_gwasbrewer_inputs(
  pool = pool,
  n_col = T,
  n_causal_fn = as.list(n_causal_direct),
  n_close = 0,
  exclude_cs = unique(
    factor_inputs$selected$pool_cs_name
  ),
  seed = 456
)


str(factor_inputs)
str(direct_inputs)

expand_pi <- function(inputs, full_snp_ids) {
  stopifnot(
    is.matrix(inputs$pi),
    length(inputs$snp_ids) == nrow(inputs$pi)
  )

  full_snp_ids <- as.character(full_snp_ids)
  small_snp_ids <- as.character(inputs$snp_ids)

  if (anyNA(full_snp_ids) || anyDuplicated(full_snp_ids)) {
    stop("Full SNP IDs must be nonmissing and unique")
  }

  if (anyNA(small_snp_ids) || anyDuplicated(small_snp_ids)) {
    stop("Input SNP IDs must be nonmissing and unique")
  }

  idx <- match(small_snp_ids, full_snp_ids)

  if (anyNA(idx)) {
    missing <- small_snp_ids[is.na(idx)]
    stop(
      length(missing),
      " selected SNP(s) are missing from the LD SNP universe: ",
      paste(head(missing, 10), collapse = ", ")
    )
  }

  out <- matrix(
    0,
    nrow = length(full_snp_ids),
    ncol = ncol(inputs$pi),
    dimnames = list(
      full_snp_ids,
      colnames(inputs$pi)
    )
  )

  out[idx, ] <- inputs$pi
  out
}

ld_snp_ids <- ref$snp_info$SNP

pi_L <- expand_pi(
  inputs = factor_inputs,
  full_snp_ids = ld_snp_ids
)

pi_theta <- expand_pi(
  inputs = direct_inputs,
  full_snp_ids = ld_snp_ids
)

dim(pi_L)
#str(pi_L)
colSums(pi_L)
sum(rowSums(pi_L) != 0)

save(factor_inputs, pi_L, pi_theta, file='gwasbrewer_inputs.RDS')
