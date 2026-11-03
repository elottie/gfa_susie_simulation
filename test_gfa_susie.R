library(GWASBrewer)
library(dplyr)
library(ggplot2)
library(viridis)
library(GFA)
library(susieR)
library(Matrix)


set.seed(100)

pdf("Rplots_3traits_2closepeak_testing.pdf", width = 8, height = 10)

#use factors from vignette

# True factor structure
myF <- matrix(c(0, -1, 1, 1, 2, 1), nrow = 3, byrow = TRUE)

myF <- apply(myF, 2, function(x){x/sqrt(sum(x^2))})

plot_factors(myF, row_names = c("T3", "T2", "T1"), col_names = c("F1", "F2")) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2], name = "Normalized\nEffect") +
  xlab("Factor") +
  ggtitle("True Factor Network") +
  scale_x_discrete(position = "top")+
  theme(axis.text.x = element_text(angle = 0),
        panel.background = element_rect(fill = "white"),
        axis.ticks = element_blank())
#ggsave('TrueFacNetwork.png')


# simulate data
load('gwasbrewer_inputs.RDS')

#str(factor_inputs)
#str(direct_inputs)

ref <- readRDS("/nfs/turbo/sph-jvmorr/GFA_susie/simulation/metab_susie/ld_ref/metsim_ld_ref_1M.RDS") 
print(head(ref$snp_info$AF))
print(nrow(ref$snp_info$AF))
dat <- GWASBrewer::sim_lf(F_mat = myF,
                          N = 50000,
                          J = nrow(pi_L),
                          h2_trait = rep(0.3, 3),
                          omega = rep(0.8, 3),
                          pi_L = pi_L,
                          pi_theta = colMeans(pi_theta),
                          R_LD = ref$R_LD,  # jean made from metsim
                          af = ref$snp_info$AF,
                          est_s = TRUE)
dat <- GWASBrewer:::calc_ld_scores(dat, ref$R_LD)

str(dat)

# trait manhat plots to see if looks like how we want ---
library(dplyr)
library(tidyr)
library(ggplot2)

stopifnot(
  nrow(dat$beta_hat) == nrow(dat$snp_info),
  all(dim(dat$beta_hat) == dim(dat$se_beta_hat))
)

trait_names <- colnames(dat$beta_hat)
if (is.null(trait_names)) {
  trait_names <- paste0("Trait ", seq_len(ncol(dat$beta_hat)))
}

plot_dat <- dat$snp_info %>%
  mutate(row_id = row_number()) %>%
  bind_cols(
    as.data.frame(dat$beta_hat) %>%
      setNames(paste0("beta_", seq_len(ncol(dat$beta_hat)))),
    as.data.frame(dat$se_beta_hat) %>%
      setNames(paste0("se_", seq_len(ncol(dat$se_beta_hat))))
  ) %>%
  pivot_longer(
    cols = matches("^(beta|se)_"),
    names_to = c(".value", "trait_number"),
    names_pattern = "(beta|se)_(\\d+)"
  ) %>%
  mutate(
    trait = factor(
      paste0("Trait ", trait_number),
      levels = paste0("Trait ", seq_along(trait_names)),
      labels = trait_names
    ),
    z = beta / se,
    p = 2 * pnorm(abs(z), lower.tail = FALSE),
    logp = -log10(p)
  )

plot_dat <- plot_dat %>%
  mutate(plot_position = row_id)

## One stacked panel per trait
traits <- unique(plot_dat$trait)
n_traits <- length(traits)

## Significance threshold
threshold <- -log10(0.05 / nrow(dat$snp_info))

## One color per block
blocks <- unique(plot_dat$block)
#block_cols <- setNames(
#  hcl.colors(length(blocks), palette = "Viridis"),
#  as.character(blocks)
#)

## Save current graphics settings and restore them afterward
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))

par(
  mfrow = c(n_traits, 1),
  mar = c(2.3, 4.3, 1.5, 1),
  oma = c(3, 1, 3, 0),
  las = 1
)

for (i in seq_along(traits)) {
  tr <- traits[i]
  d <- plot_dat[plot_dat$trait == tr, ]

  plot(
    d$plot_position,
    d$logp,
    pch = 16,
    cex = 0.7,
 #   col = adjustcolor(
 #     block_cols[as.character(d$block)],
 #     alpha.f = 0.75
 #   ),
    xlab = "",
    ylab = if (i == ceiling(n_traits / 2)) expression(-log[10](p)) else "",
    main = tr,
    axes = FALSE
  )

  axis(2)
  box()

  ## Display the x-axis only on the bottom panel
  if (i == n_traits) {
    axis(1)
  }

  abline(h = threshold, lty = 2, col = "red")
}

mtext("Variant index", side = 1, outer = TRUE, line = 1)
mtext(
  "Simulated trait manhattans",
  side = 3,
  outer = TRUE,
  line = 1,
  font = 2
)



#trait_manhat <- ggplot(plot_dat, aes(plot_position, logp)) +
#  geom_point(aes(color = factor(block)), size = 0.7, alpha = 0.75) +
#  facet_wrap(~trait, ncol = 1, scales = "free_y") +
#  geom_hline(
#    yintercept = -log10(0.05 / nrow(dat$snp_info)),
#    linetype = "dashed",
#    color = "red"
#  ) +
#  scale_color_viridis_d(guide = "none") +
#  labs(
#    x = "Variant index",
#    y = expression(-log[10](p)),
#    title = "Association results for three traits"
#  ) +
#  theme_classic()
#trait_manhat

#fit regular gfa

Z_hat = with(dat, beta_hat/s_estimate)
ix <- sim_ld_prune(dat, R_LD = ref$R_LD, r2_thresh = 0.01)
Z_hat_indep <- Z_hat[ix,]

fit_toy <- gfa_fit(Z_hat = Z_hat_indep,
                   N = rep(50000, 3))

p_reg <- cbind(fit_toy$F_hat,fit_toy$F_hat_single) %>% 
  plot_factors(.) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2]) +
  xlab("Factor") +
  ylab("Trait") +
  ggtitle("Regular GFA Fit on Close-Peak Data") 
p_reg

# look for what is L_hat_marginal
j <- 1

f_norm <- sqrt(sum(fit_toy$fit$F_pm[, j]^2))  # happens so f_hat = f_pm / norm(f_pm)
scaling_factor <- f_norm / fit_toy$scale[1]  # since we have LF^T, and f_hat normalized by dividing by norm, have to multiply L column by norm of factor.  then divide by scale, not clear on why they did that

f_norm
scaling_factor

all.equal(
  fit_toy$F_hat[, j],
  fit_toy$fit$F_pm[, j] / f_norm
)

all.equal(
  fit_toy$L_hat[, j],
  fit_toy$fit$L_pm[, j] * scaling_factor
)
# so L_hat_marginal = fit_toy$fit$L_pm[, j] * scaling_factor
L_hat_marg <- fit_toy$fit$L_pm[, j] * scaling_factor
plot(L_hat_marg)
plot(fit_toy$L_hat[,1])

# susie part ---

# Step 1: run SuSiE-RSS on each trait separately, one LD block at a time.
#
# Inputs
#   beta_hat, se : J x T matrices (dat_fixed$beta_hat, dat_fixed$s_estimate)
#   ld_list      : list of LD (correlation) matrices, one per block, in the same
#                  variant order as the rows of beta_hat (ld_mat_list)
#   n            : GWAS sample size (50000 in the simulation)
#
# Output (list)
#   pip   : J x T matrix of per-variant PIPs (0 for blocks that were skipped)
#   cs    : data.frame with one row per credible set, per trait
#           (trait, block, cs, size, min_abs_corr, lead, lead_pip, members)
#           `lead` and `members` are row indices into beta_hat
#   fits  : per-trait list of per-block susie fits (only if keep_fits = TRUE)
#
# Notes
#   * Blocks whose strongest |z| is below `skip_z` are not fitted (PIPs set to 0).
#     This is a speed shortcut; set skip_z = 0 to fit every block. Because the
#     union over traits is taken later, a block skipped for one trait is still
#     picked up if it matters for another.
#   * estimate_residual_variance = FALSE is susie_rss's default and is kept
#     explicit here. With z-scores plus the population LD this is the safe choice.
#   * L is the maximum number of effects per block, not per genome.



# ---------------------------------------------------------------------------
# Usage with the GWASBrewer simulation
#
# should L be altered???
res1 <- GFA:::run_susie_per_trait(dat$beta_hat, dat$s_estimate, ref$R_LD,
                               n = 50000, L = 10, ncores = 1)
#   res1$cs                      # credible sets, per trait
#   colSums(res1$pip > 0.5)      # variants with PIP > 0.5 per trait
#
# Optional sanity check against the truth, if your simulation object stores the
# true joint effects (I believe the element is dat_fixed$beta_joint, J x T; adjust the
# name if yours differs):
   truth <- lapply(seq_len(ncol(dat$beta_joint)), function(t) which(dat$beta_joint[, t] != 0))
   for (t in seq_along(truth)) {
     cs_t <- subset(res1$cs, trait == colnames(res1$pip)[t])
     covered <- vapply(truth[[t]], function(v) any(vapply(cs_t$members, function(m) v %in% m, TRUE)), TRUE)
     cat(colnames(res1$pip)[t], ": ", sum(covered), "/", length(truth[[t]]),
         " true causal variants in a credible set; ", nrow(cs_t), " sets total\n", sep = "")
   }

# collect top snps from each credible set. For now ignore how big the pip is for the top snp
   # alternate: collect all snps with pip > C

#ix <- unique(sort(res1$cs$lead))
#R_ix <- GWASBrewer::sim_extract_ld(dat_fixed = dat_fixed, index = ix, R_LD = ld_mat_list)

a <- apply(res1$pip, 1, max)
ix <- which(a > 0.05) # threshold chosen based on number of variants. using 0.5 gives 44 variants. this gives 222 which seemed like a good number to test with.
blocks <- dat$snp_info$block[ix]
R_ix <- lapply(unique(blocks), function(b){
  ixb <- ix[blocks == b]
  if(length(ixb) == 1){
    R <- matrix(1, nrow = 1, ncol = 1)
    rownames(R) <- ixb
    colnames(R) <- ixb
    return(R)
  }
  sim_extract_ld(dat = dat, index = ixb, R_LD = ref$R_LD)
})

print(length(ix))

# whitening --

w <- GFA:::make_whitening(R_ix, keep_frac = 0)
fn  <- GFA:::susie_ebnm_fn(w$design, L = 10)
Z <- lapply(unique(blocks), function(b){
  ixb <- ix[blocks == b]
  Z_hat[ixb,]
})
Zwhitened <- Map(function(W, z){ W %*% z}, w$W, Z)
Zwhitened <- Reduce(rbind, Zwhitened)
fit1 <- gfa_fit(Z_hat = Zwhitened, N = rep(50000, 3), params = list(ebnm_fn_L = fn))

# construct scaled F_hat and L_hat
F_hat_marg <- matrix(data=NA,nrow=nrow(fit1$fit$F_pm),ncol=ncol(fit1$fit$F_pm))
L_hat_marg <- matrix(data=NA,nrow=nrow(fit1$fit$L_pm),ncol=ncol(fit1$fit$L_pm))

for (j in 1:ncol(fit1$fit$F_pm)){
  f_norm <- sqrt(sum(fit1$fit$F_pm[, j]^2))  # happens so f_hat = f_pm / norm(f_pm)
  scaling_factor <- f_norm / fit1$scale[1]  # since we have LF^T, and f_hat normalized by dividing by norm, have to multiply L column by norm of factor.  then divide by scale, not clear on why they did that

  # so L_hat_marginal = fit1$fit$L_pm[, j] * scaling_factor
  L_hat_marg[,j] <- fit1$fit$L_pm[, j] * scaling_factor

  F_hat_marg[,j] <- fit1$fit$F_pm[, j] / f_norm
}

plot(L_hat_marg[,1])

p_whiten <- F_hat_marg %>%
  plot_factors(.) +
  scale_fill_gradient2(low = viridis(3)[1], high = viridis(3)[2]) +
  xlab("Factor") +
  ylab("Trait") +
  ggtitle("Whitened GFA Fit on Close-Peak Data (Questionable SuSiE conv)")

p_whiten

dev.off()
