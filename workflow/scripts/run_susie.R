log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressMessages({
  library(data.table)
  library(susieR)
})

tol <- 1e-8
susie_cfg <- snakemake@params[["susie"]]
transform <- snakemake@params[["transform"]]

write_empty <- function() {
  fwrite(data.table(snp = character(), chr = character(), bp = integer(),
                    effect = character(), pip = numeric()),
         snakemake@output[["pip"]], sep = "\t")
  fwrite(data.table(cs = integer(), snp = character(), chr = character(),
                    bp = integer(), effect = character(), pip = numeric()),
         snakemake@output[["cs"]], sep = "\t")
  saveRDS(NULL, snakemake@output[["rds"]])
}

# whitening matrix and canonical sample order
wh <- readRDS(snakemake@input[["linv"]])
linv <- wh$Linv
ids <- wh$ids
n <- length(ids)

# transformed phenotype, aligned to ids
phen <- fread(snakemake@input[["phenotype"]], header = TRUE)
trait <- colnames(phen)[3]
phen_key <- paste(phen[[1]], phen[[2]], sep = ":")
pidx <- match(ids, phen_key)
if (anyNA(pidx)) stop("transformed phenotype missing samples present in GRM ids")
y <- as.numeric(phen[[trait]][pidx])

# fixed design Xf: intercept only for iint (PCs already absorbed); else intercept + PCs
if (transform == "iint") {
  xf <- matrix(1, n, 1)
} else {
  qcov <- fread(snakemake@input[["qcovar"]], header = TRUE)
  qkey <- paste(qcov[[1]], qcov[[2]], sep = ":")
  qidx <- match(ids, qkey)
  if (anyNA(qidx)) stop("qcovar missing samples present in GRM ids")
  pcs <- as.matrix(qcov[qidx, -(1:2), with = FALSE])
  xf <- cbind(1, pcs)
}

# region genotypes (plink --recode A); first 6 cols are FID IID PAT MAT SEX PHENOTYPE
geno <- fread(snakemake@input[["geno"]], header = TRUE, na.strings = "NA")
if (ncol(geno) <= 6L) {
  write_empty()
  quit(save = "no")
}
geno_key <- paste(geno[[1]], geno[[2]], sep = ":")
gidx <- match(ids, geno_key)
if (anyNA(gidx)) stop("region genotypes missing samples present in GRM ids")
raw <- as.matrix(geno[gidx, 7:ncol(geno), with = FALSE])
raw_cols <- colnames(geno)[7:ncol(geno)]

# map plink --recode A column "<snpid>_<countedallele>" back to bim snp/chr/bp
bim <- fread(snakemake@input[["bim"]], header = FALSE)
colnames(bim) <- c("chr", "snp", "cm", "bp", "a1", "a2")
bim_key <- paste(bim$snp, bim$a1, sep = "_")
bidx <- match(raw_cols, bim_key)
if (anyNA(bidx)) stop("could not map some genotype columns to .bim entries")
snp_id <- bim$snp[bidx]
snp_chr <- as.character(bim$chr[bidx])
snp_bp <- bim$bp[bidx]

# NC (NOIACenter) encoding: empirical freqs, center only, NOIA dominance, missing -> 0
m <- ncol(raw)
add_list <- vector("list", m)
dom_list <- vector("list", m)
add_ok <- logical(m)
dom_ok <- logical(m)
for (j in seq_len(m)) {
  g <- raw[, j]
  obs <- g[!is.na(g)]
  nn <- length(obs)
  if (nn == 0) next
  pAA <- sum(obs == 2) / nn  # A1A1
  pAa <- sum(obs == 1) / nn  # A1A2
  paa <- sum(obs == 0) / nn  # A2A2
  p <- pAA + 0.5 * pAa

  add <- g - 2 * p
  add[is.na(add)] <- 0
  if (stats::var(add) > tol) {
    add_list[[j]] <- add
    add_ok[j] <- TRUE
  }

  denom <- pAA + paa - (pAA - paa)^2
  if (denom >= tol) {
    cAA <- -2 * paa * pAa / denom
    cAa <- 4 * pAA * paa / denom
    caa <- -2 * pAA * pAa / denom
    code <- c(caa, cAa, cAA)[g + 1L]  # g in {0,1,2}
    mean_dom <- paa * caa + pAa * cAa + pAA * cAA
    dom <- code - mean_dom
    dom[is.na(dom)] <- 0
    if (stats::var(dom) > tol) {
      dom_list[[j]] <- dom
      dom_ok[j] <- TRUE
    }
  }
}

cols <- list()
meta_snp <- character()
meta_chr <- character()
meta_bp <- integer()
meta_eff <- character()
for (j in which(add_ok)) {
  cols[[length(cols) + 1L]] <- add_list[[j]]
  meta_snp <- c(meta_snp, snp_id[j]); meta_chr <- c(meta_chr, snp_chr[j])
  meta_bp <- c(meta_bp, snp_bp[j]); meta_eff <- c(meta_eff, "add")
}
for (j in which(dom_ok)) {
  cols[[length(cols) + 1L]] <- dom_list[[j]]
  meta_snp <- c(meta_snp, snp_id[j]); meta_chr <- c(meta_chr, snp_chr[j])
  meta_bp <- c(meta_bp, snp_bp[j]); meta_eff <- c(meta_eff, "dom")
}
if (length(cols) == 0L) {
  write_empty()
  quit(save = "no")
}
z <- do.call(cbind, cols)

# whiten, then project out the whitened fixed-effect columns (GLS equivalent of gelex's P)
xw <- linv %*% z
yw <- as.numeric(linv %*% y)
cw <- linv %*% xf
resid_out <- function(mat) mat - cw %*% solve(crossprod(cw), crossprod(cw, mat))
xw <- resid_out(xw)
yw <- as.numeric(resid_out(matrix(yw, ncol = 1)))

fit <- susie(xw, yw,
             L = susie_cfg[["L"]],
             coverage = susie_cfg[["coverage"]],
             min_abs_corr = susie_cfg[["min_abs_corr"]],
             standardize = susie_cfg[["standardize"]],
             intercept = FALSE)

pip <- data.table(snp = meta_snp, chr = meta_chr, bp = meta_bp,
                  effect = meta_eff, pip = as.numeric(fit$pip))
fwrite(pip, snakemake@output[["pip"]], sep = "\t")

cs_rows <- list()
if (!is.null(fit$sets$cs)) {
  for (k in seq_along(fit$sets$cs)) {
    idx <- fit$sets$cs[[k]]
    cs_name <- names(fit$sets$cs)[k]
    cs_rows[[length(cs_rows) + 1L]] <- data.table(
      cs = cs_name, snp = meta_snp[idx], chr = meta_chr[idx],
      bp = meta_bp[idx], effect = meta_eff[idx], pip = as.numeric(fit$pip[idx]))
  }
}
cs_dt <- if (length(cs_rows)) rbindlist(cs_rows) else data.table(
  cs = character(), snp = character(), chr = character(),
  bp = integer(), effect = character(), pip = numeric())
fwrite(cs_dt, snakemake@output[["cs"]], sep = "\t")

saveRDS(fit, snakemake@output[["rds"]])
