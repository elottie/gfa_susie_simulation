library(data.table)
library(Matrix)

## ---------- 1. Load per-chromosome objects into one pool ----------
# Stores only sparse info: which CS overlap, and lead-variant r2 for non-overlapping
# pairs above r2_min. Pass r2_close >= r2_min to make_gwasbrewer_inputs().
load_cs_pool <- function(dir, prefix, ref_info, chrs = 1:22, r2_min = 0.05) {
  mk <- function(pos, a, b) paste(pos, pmin(a, b), pmax(a, b), sep = ":")

  # reference variants keyed on chromosome + position + unordered alleles
  ref_info <- as.data.table(ref_info)
  ref_info[, chr_num := sub("^chr", "", as.character(chr))]
  ref_info[, key := mk(pos, A1, A2)]

  cs_all <- ov_all <- ed_all <- stats <- list()
  off <- 0L
  for (ch in chrs) {
    stem <- file.path(dir, paste0(prefix, ".chr", ch))
    f <- paste0(stem, c(".cs_list.RDS", ".cs_lead_ld.RDS"))
    if (!all(file.exists(f))) { warning("Missing files for chr", ch, "; skipping"); next }
    cs <- readRDS(f[1]); R <- readRDS(f[2])
    stopifnot(identical(names(cs), rownames(R)))

    ref_ch <- ref_info[chr_num == as.character(ch)]
    lookup <- setNames(ref_ch$SNP, ref_ch$key)

    ## ---- trim each credible set to variants in the reference ----
    n_var_before <- sum(lengths(lapply(cs, `[[`, "id")))
    cs_t <- lapply(cs, function(x) {
      sp  <- tstrsplit(x$id, ":", fixed = TRUE)
      rid <- unname(lookup[mk(as.integer(sp[[2]]), sp[[3]], sp[[4]])])
      k   <- !is.na(rid)
      if (!any(k)) return(NULL)
      x$lead_kept <- x$lead_id %in% x$id[k]        # was the original lead retained?
      p <- x$pip[k]
      if (sum(p) <= 0) p <- rep(1, length(p))
      x$id          <- rid[k]                       # now the reference SNP ids
      x$position    <- x$position[k]
      x$effect_size <- x$effect_size[k]
      x$pip         <- p / sum(p)                   # renormalize over kept variants
      j <- which.max(x$pip)
      x$lead_id  <- x$id[j]
      x$lead_pos <- x$position[j]
      x
    })
    keep <- !vapply(cs_t, is.null, TRUE)
    cs_t <- cs_t[keep]
    stats[[as.character(ch)]] <- data.table(
      chr = ch, n_cs = length(cs), n_cs_dropped = sum(!keep),
      n_cs_trimmed = sum(vapply(cs_t, function(x) length(x$id), 1L) <
                         lengths(lapply(cs[keep], `[[`, "id"))),
      n_var_before = n_var_before,
      n_var_after  = sum(lengths(lapply(cs_t, `[[`, "id"))))
    if (length(cs_t) == 0) next

    ## ---- overlap, recomputed from the trimmed sets ----
    ids  <- lapply(cs_t, `[[`, "id")
    vars <- unique(unlist(ids))
    M <- sparseMatrix(i = rep(seq_along(ids), lengths(ids)), j = match(unlist(ids), vars),
                      x = 1, dims = c(length(ids), length(vars)))
    O  <- as(tcrossprod(M), "generalMatrix")
    ov <- as.data.table(summary(O))[i != j, .(i, j)]

    ## ---- close pairs: lead-variant r2, only among credible sets whose original lead was retained ----
    R2 <- R[keep, keep, drop = FALSE]^2
    lk <- vapply(cs_t, `[[`, TRUE, "lead_kept")
    R2[!lk, ] <- 0; R2[, !lk] <- 0
    w  <- which(R2 > r2_min, arr.ind = TRUE)
    ed <- data.table(i = w[, 1], j = w[, 2], r2 = R2[w])[i != j]
    ed <- ed[!ov, on = .(i, j)]

    ov[, `:=`(i = i + off, j = j + off)]
    ed[, `:=`(i = i + off, j = j + off)]
    cs_all[[as.character(ch)]] <- cs_t
    ov_all[[as.character(ch)]] <- ov
    ed_all[[as.character(ch)]] <- ed
    off <- off + length(cs_t)
  }
  st <- rbindlist(stats)
  message(sum(st$n_cs_dropped), " of ", sum(st$n_cs), " credible sets dropped (no variants in reference); ",
          sum(st$n_cs_trimmed), " trimmed")
  list(cs = do.call(c, unname(cs_all)),
       overlaps = rbindlist(ov_all), edges = rbindlist(ed_all), stats = st)
}

## ---------- 2. Build pi and effect functions from selected variants ----------
# sel: data.table with columns col, id, chr, pos, effect (one row per causal variant)
finalize_inputs <- function(sel, n_col, snp_ids = NULL) {
  sel <- copy(sel)
  if (is.null(snp_ids))
    snp_ids <- unique(sel[order(suppressWarnings(as.integer(chr)), pos), id])
  sel[, row := match(id, snp_ids)]
  if (anyNA(sel$row)) stop("Some selected variants are not in snp_ids")

  pi <- matrix(0, length(snp_ids), n_col)
  pi[cbind(sel$row, sel$col)] <- 1

  fns <- lapply(seq_len(n_col), function(cc) {
    e <- sel[col == cc][order(row), effect]       # effects in row order
    function(n, ...) {
      if (n != length(e)) stop("Effect function expected ", length(e), " effects, asked for ", n)
      e
    }
  })
  list(pi = pi, snp_effect_function = fns, selected = sel, snp_ids = snp_ids)
}

## ---------- 3. Main function ----------
pick <- function(x) x[sample.int(length(x), 1)]   # safe even when length(x) == 1

## Resolve n_causal_fn into a vector of per-column counts.
##  - a single function: called once as f(n_col), returns n_col values (as before)
##  - a list of length n_col: element k is either a function, called as f(1)
##    to give the count for column k, or a fixed number
resolve_n_causal <- function(n_causal_fn, n_col) {
  if (is.function(n_causal_fn)) {
    N <- n_causal_fn(n_col)
  } else if (is.list(n_causal_fn)) {
    if (length(n_causal_fn) != n_col)
      stop("n_causal_fn is a list of length ", length(n_causal_fn), " but n_col = ", n_col)
    N <- vapply(n_causal_fn, function(f) {
      v <- if (is.function(f)) f(1) else f
      if (length(v) != 1) stop("Each element of n_causal_fn must give exactly one count")
      as.numeric(v)
    }, numeric(1))
  } else {
    stop("n_causal_fn must be a function or a list of functions/numbers")
  }
  if (length(N) != n_col || anyNA(N) || any(N < 0) || any(N != round(N)))
    stop("n_causal_fn must produce n_col non-negative integer counts")
  as.integer(N)
}


make_gwasbrewer_inputs <- function(pool, n_col, n_causal_fn,
                                   n_close = 0, r2_close = 0.5,
                                   snp_ids = NULL, exclude_cs = character(0),
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cs  <- pool$cs
  ncs <- length(cs)
  ed  <- pool$edges[r2 > r2_close]
  nbr <- split(ed$j, factor(ed$i, levels = seq_len(ncs)))                 # close neighbors
  ovl <- split(pool$overlaps$j, factor(pool$overlaps$i, levels = seq_len(ncs)))
  has_nbr <- lengths(nbr) > 0

  # blocked = already selected, overlapping a selected CS, excluded, or not usable
  blocked <- rep(FALSE, ncs)
  if (!is.null(snp_ids))
    blocked <- !vapply(cs, function(x) all(x$id %in% snp_ids), TRUE)
  block <- function(s) blocked[c(s, ovl[[s]])] <<- TRUE
  ex <- match(exclude_cs, names(cs))
  for (s in ex[!is.na(ex)]) block(s)

  ## Step 1: number of causal variants per column
  N <- resolve_n_causal(n_causal_fn, n_col)
  assign <- replicate(n_col, integer(0), simplify = FALSE)

  ## Step 2: close pairs across different columns
  cap <- N
  pairs <- list()
  for (b in seq_len(n_close)) {
    cols <- which(cap > 0)
    if (length(cols) < 2) { warning("Only placed ", b - 1, " of ", n_close, " close pairs (column capacity)"); break }
    ab <- sample(cols, 2)

    cand <- which(!blocked & has_nbr)
    s1 <- NA_integer_
    while (length(cand) > 0) {                       # reject CS with no free neighbor
      s <- pick(cand)
      free <- nbr[[s]][!blocked[nbr[[s]]]]
      if (length(free) > 0) { s1 <- s; break }
      cand <- setdiff(cand, s)
    }
    if (is.na(s1)) { warning("Only placed ", b - 1, " of ", n_close, " close pairs (no candidates left)"); break }
    s2 <- pick(free)

    block(s1); block(s2)
    assign[[ab[1]]] <- c(assign[[ab[1]]], s1)
    assign[[ab[2]]] <- c(assign[[ab[2]]], s2)
    cap[ab] <- cap[ab] - 1L
    pairs[[b]] <- data.table(col_a = ab[1], col_b = ab[2],
                             cs_a = names(cs)[s1], cs_b = names(cs)[s2],
                             r2 = ed[i == s1 & j == s2, r2])
  }

  ## Step 3: fill remaining causal variants with unconstrained, non-overlapping CS
  for (cc in sample.int(n_col)) {
    for (r in seq_len(N[cc] - length(assign[[cc]]))) {
      cand <- which(!blocked)
      if (length(cand) == 0) stop("Ran out of non-overlapping credible sets; reduce the number of causal variants")
      s <- pick(cand); block(s)
      assign[[cc]] <- c(assign[[cc]], s)
    }
  }

  ## Step 4: one variant per CS, drawn by PIP, effect = posterior mean given inclusion
  sel <- rbindlist(lapply(seq_len(n_col), function(cc) {
    rbindlist(lapply(assign[[cc]], function(s) {
      x <- cs[[s]]
      v <- sample.int(length(x$id), 1, prob = x$pip)
      data.table(col = cc, pool_cs_name = names(cs)[s], cs_key = x$cs_key, metabolite = x$metabolite,
                 chr = x$chromosome, id = x$id[v], pos = x$position[v],
                 effect = x$effect_size[v])
    }))
  }))

  c(finalize_inputs(sel, n_col, snp_ids),
    list(N = N, pairs = rbindlist(pairs)))
}
