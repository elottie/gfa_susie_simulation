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
                          snp_effect_function_L = factor_inputs$snp_effect_function,
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
    logp <- -(log(2) + pnorm(
			     abs(beta / se),
			     lower.tail = FALSE,
			     log.p = TRUE
		 )) / log(10)
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



dev.off()
