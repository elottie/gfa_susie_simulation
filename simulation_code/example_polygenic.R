source("build_gwasbrewer_input_polygenic.R")
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
#traits <- vapply(
#  pool$cs,
#  function(x) {
#    value <- unique(as.character(x$metabolite))
#    value <- value[!is.na(value)]
#
#    if (length(value) != 1) {
#      stop("Each credible set must have exactly one non-missing metabolite")
#    }
#
#    value
#  },
#  character(1)
##)

# Count credible sets by trait
#n_cs_by_trait <- setNames(
#  as.integer(table(traits)),
#  names(table(traits))
#)
#n_cs <- unname(n_cs_by_trait)
#saveRDS(n_cs, file = "credible_set_counts.RDS")
n_cs <- readRDS("../credible_set_counts.RDS")



buffer <- 1.5
factor_fraction <- 0.5
n_causal_fn <- function(n_col){
    x <- sample(n_cs, size = n_col, replace = TRUE)
    x <- x*buffer*factor_fraction
    return(ceiling(x))
}


# let K = 3 factors for example
K <- 3

L <- make_gwasbrewer_inputs(
  pool = pool,
  n_col = K,
  n_causal_fn = n_causal_fn,
  n_close = 4,
  snp_ids = ref$snp_info$SNP,
  seed = 123
)

L_poly <- make_gwasbrewer_inputs(
  pool = pool,
  n_col = K,
  n_causal_fn = n_causal_fn,
  n_close = 4,
  snp_ids = ref$snp_info$SNP,
  seed = 123,
  polygenic_pi = 1000/nrow(ref$snp_inf),
  polygenic_proportion = 0.5,
)

str(L)
str(L_poly)

saveRDS(L,'example_L_direct_gb_inputs.RDS')
saveRDS(L_poly,'example_L_poly_gb_inputs.RDS')
