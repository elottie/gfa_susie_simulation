library(GWASBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(GFA)
library(susieR)
library(Matrix)

set.seed(123)

# --- read gb input data obj ---
gb_inp_string <- "20Trt_FacProp0_NClose4_Poly0.5"
gb_inp_filename <- paste0('GBInp_',gb_inp_string,'.RDS')
gb_inp <- readRDS(gb_inp_filename)
pi_theta <- gb_inp$pi
str(pi_theta)
M <- ncol(pi_theta)

ref <- readRDS("/nfs/turbo/sph-jvmorr/GFA_susie/simulation/metab_susie/ld_ref/metsim_ld_ref_1M.RDS") 

# read in heritabilities so sim data can use them
met_dict <- readRDS("metDict.RDS")
herit <- met_dict$GREML_PC5_Herit
h2_used <- sample(herit, size = M)

pdf(paste0("Rplots_",gb_inp_string,".pdf"), width = 8, height = 10)

# --- options ---
# specify true facs
F_mat = matrix(0, nrow = M, ncol = 2)  # no facs

dat <- sim_lf(
  F_mat = F_mat,
  N = 10000,
  J = nrow(pi_theta),
  h2_trait = h2_used,
  omega = rep(0, M),        # No heritability mediated by factors
  pi_L = 0,                 # Placeholder; factor effects are inactive
  pi_theta = pi_theta,          # Direct-effect causal SNP proportions
  snp_effect_function_theta = gb_inp$snp_effect_function,
  R_LD = ref$R_LD,          # jean made from metsim
  af = ref$snp_info$AF,
  est_s = TRUE)
dat <- GWASBrewer:::calc_ld_scores(dat, ref$R_LD)
str(dat)

# --- trait manhat plots to see if looks like how we want ---
# make the plotting data
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
    logp = -(log(2) + pnorm(-abs(z),log.p = TRUE)) / log(10)
  )

plot_dat <- plot_dat %>%
  mutate(plot_position = row_id)

# added extra code, this plots the upper and lower bound of the pairs region so we can see close peaks
stopifnot(length(gb_inp$snp_ids) == nrow(dat$beta_hat))

# Extract chromosome and position from "1:777957:G:A"
snp_parts <- strsplit(gb_inp$snp_ids, ":", fixed = TRUE)

snp_chr <- vapply(snp_parts, `[`, character(1), 1L)
snp_pos <- as.numeric(vapply(snp_parts, `[`, character(1), 2L))

# Collect both sides of each pair, keeping their trait assignments
pairs <- as.data.frame(gb_inp$pairs)

cs <- unique(rbind(
  data.frame(trait_number = pairs$col_a, cs = pairs$cs_a),
  data.frame(trait_number = pairs$col_b, cs = pairs$cs_b)
))

cs_parts <- strsplit(cs$cs, "|", fixed = TRUE)

cs$chr <- vapply(cs_parts, `[`, character(1), 2L)
cs$lower <- as.numeric(vapply(cs_parts, `[`, character(1), 3L))
cs$upper <- as.numeric(vapply(cs_parts, `[`, character(1), 4L))

# Map each interval to the first and last matching SNP row
cs$x_bounds <- lapply(seq_len(nrow(cs)), function(j) {
  rows <- which(
    snp_chr == cs$chr[j] &
      snp_pos >= cs$lower[j] &
      snp_pos <= cs$upper[j]
  )

  if (length(rows)) range(rows) else numeric(0)
})

# now begin the plotting
## One stacked panel per trait
traits <- unique(plot_dat$trait)

## Significance threshold
threshold <- -log10(0.05 / nrow(dat$snp_info))

## Save current graphics settings and restore them afterward
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))

par(
  mfrow = c(5, 2),
  mar = c(2.3, 4.3, 1.5, 1),
  oma = c(3, 1, 3, 0),
  las = 1
)

plots_per_page <- 10L
n_plots <- length(traits)

panel_titles <- paste0(
  trait_names,
  ":  n_cs ", gb_inp$N,
  ", h2 ", signif(h2_used, 3)
)

for (i in seq_along(traits)) {
  tr <- traits[i]
  d <- plot_dat[plot_dat$trait == tr, ]

  # Trait column corresponding to this plot
  trait_number <- match(as.character(tr), trait_names)

  # Panel number within the current page: 1 through 10
  panel <- ((i - 1L) %% plots_per_page) + 1L

  plot(
    d$plot_position,
    d$logp,
    pch = 16,
    cex = 0.7,
    xlab = "",
    ylab = "",
    main = panel_titles[trait_number],
    axes = FALSE
  )

  bounds <- unlist(
    cs$x_bounds[cs$trait_number == trait_number],
    use.names = FALSE
  )

  if (length(bounds)) {
    abline(v = unique(bounds), col = "blue", lty = 3, lwd = 1)
  }

  axis(2)
  box()

  # Bottom row on a full page: panels 9 and 10
  page_start <- ((i - 1L) %/% plots_per_page) * plots_per_page
  panels_this_page <- min(plots_per_page, n_plots - page_start)

  # No further plot below this panel in its column
  if (panel + 2L > panels_this_page) {
    axis(1)
  }

  abline(h = threshold, lty = 2, col = "red")

  # Add outer labels before moving to the next page
  if (panel == plots_per_page || i == n_plots) {
    mtext("Variant index", side = 1, outer = TRUE, line = 1)
    mtext(
      expression(-log[10](p)),
      side = 2,
      outer = TRUE,
      line = 0
    )
    mtext(
      paste0("Simulated Trait Manhattans:  ",gb_inp_string),
      side = 3,
      outer = TRUE,
      line = 1,
      font = 2
    )
  }
}

dev.off()
