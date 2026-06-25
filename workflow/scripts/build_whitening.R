log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Read a GRM written by gelex: <prefix>.bin holds the lower triangle (incl. diag)
# row-by-row as float32; <prefix>.id gives FID/IID sample order.
read_grm <- function(prefix) {
  ids_df <- read.table(paste0(prefix, ".id"), header = FALSE, sep = "\t",
                       colClasses = "character", stringsAsFactors = FALSE)
  ids <- paste(ids_df[[1]], ids_df[[2]], sep = ":")
  n <- length(ids)
  con <- file(paste0(prefix, ".bin"), "rb")
  on.exit(close(con))
  vals <- readBin(con, what = "numeric", n = n * (n + 1) / 2, size = 4)
  g <- matrix(0, n, n)
  pos <- 1L
  for (i in seq_len(n)) {
    g[i, seq_len(i)] <- vals[pos:(pos + i - 1L)]
    pos <- pos + i
  }
  full <- g + t(g)
  diag(full) <- diag(g)
  list(G = full, ids = ids)
}

# Reconstruct the LOCO GRM exactly as gelex does (loco_reader.cpp):
# (G_whole - G_chr) / (k_whole - k_chr), with k = trace / n.
loco_grm <- function(whole, chr_grm) {
  n <- nrow(whole)
  k_loco <- (sum(diag(whole)) - sum(diag(chr_grm))) / n
  if (k_loco <= 0) stop("LOCO denominator <= 0")
  (whole - chr_grm) / k_loco
}

prefix <- dirname(snakemake@input[["summary"]])
chr <- as.character(snakemake@params[["chr"]])

add_whole <- read_grm(sub("\\.bin$", "", snakemake@input[["grm_add"]][1]))
add_chr <- read_grm(sub("\\.bin$", "", snakemake@input[["grm_add_chr"]][1]))
dom_whole <- read_grm(sub("\\.bin$", "", snakemake@input[["grm_dom"]][1]))
dom_chr <- read_grm(sub("\\.bin$", "", snakemake@input[["grm_dom_chr"]][1]))

ids <- add_whole$ids
stopifnot(identical(ids, dom_whole$ids))
stopifnot(identical(ids, add_chr$ids))
stopifnot(identical(ids, dom_chr$ids))

g_add <- loco_grm(add_whole$G, add_chr$G)
g_dom <- loco_grm(dom_whole$G, dom_chr$G)

# Per-chromosome variance components from gelex reml --loco.
summ <- read.table(snakemake@input[["summary"]], header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
sub <- summ[as.character(summ$chr) == chr, ]
if (nrow(sub) == 0) stop(sprintf("chr %s not found in reml summary", chr))
if (any(sub$converged != 1)) stop(sprintf("reml did not converge for chr %s", chr))

genetic <- sub[sub$term != "Residual" & sub$type == "variance", ]
residual <- sub[sub$term == "Residual", ]
if (nrow(genetic) != 2L) stop("expected exactly 2 genetic variance components (add, dom)")
# gelex preserves --grm order: first = additive, second = dominance.
sa2 <- genetic$estimate[1]
sd2 <- genetic$estimate[2]
se2 <- residual$estimate[1]

n <- nrow(g_add)
v <- sa2 * g_add + sd2 * g_dom + se2 * diag(n)

upper <- tryCatch(chol(v), error = function(e) {
  jitter <- 1e-6 * mean(diag(v))
  warning(sprintf("V not positive-definite, adding jitter %g", jitter))
  chol(v + jitter * diag(n))
})
# V = L L', L = t(upper); whitening matrix W = L^-1.
linv <- solve(t(upper))

saveRDS(list(Linv = linv, ids = ids), snakemake@output[["linv"]])
