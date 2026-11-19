source("build_gwasbrewer_input_polygenic.R")
options(warn=1)

# --- options ----
seed <- 123

# assume we didn't detect all cs
buffer <- 1.5
# how many cs you want to act thru factors vs traits
# default 0.5, make 1 to make all cs thru facs or traits
factor_fraction <- 1

# set number of traits or facts
# let K = 3 factors for example
#K <- 3
# let M = 20 traits for example
M <- 20

n_close <- 4

# 0.5 too high
polygen_prop <- 0.7

# --- load data ---
# load credible sets and their relationships
ref <- readRDS("/nfs/turbo/sph-jvmorr/GFA_susie/simulation/metab_susie/ld_ref/metsim_ld_ref_1M.RDS")
pool <- load_cs_pool(
  dir = "/nfs/turbo/sph-jvmorr/GFA_susie/simulation/metab_susie/credible_set_info/cs_objects/",
  prefix = "metsim_2022.susie.all_in.pval5e-8",
  chrs = 1:22,
  r2_min = 0.05,
  ref_info = ref$snp_info
)

n_cs <- readRDS("../credible_set_counts.RDS")

# --- make causal eff func ---
n_causal_fn <- function(n_col){
    x <- sample(n_cs, size = n_col, replace = TRUE)
    x <- x*buffer*factor_fraction
    return(ceiling(x))
}

# --- make gwasbrewer inputs ---
gb_inp <- make_gwasbrewer_inputs(
  pool = pool,
  n_col = M,
  n_causal_fn = n_causal_fn,
  n_close = n_close,
  snp_ids = ref$snp_info$SNP,
  seed = seed
)

gb_inp_poly <- make_gwasbrewer_inputs(
  pool = pool,
  n_col = M,
  n_causal_fn = n_causal_fn,
  n_close = n_close,
  snp_ids = ref$snp_info$SNP,
  seed = seed,
  polygenic_pi = 1000/nrow(ref$snp_inf),
  polygenic_proportion = polygen_prop,
)

str(gb_inp)
str(gb_inp_poly)

# --- save inputs ---
saveRDS(gb_inp,'GBInp_20Trt_FacProp0_NClose4_Poly0.RDS')
saveRDS(gb_inp_poly,paste0('GBInp_20Trt_FacProp0_NClose4_Poly',polygen_prop,'.RDS'))
