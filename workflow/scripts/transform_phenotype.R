log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

rint <- function(x, k = 0.375) {
  n <- length(x)
  qnorm((rank(x, ties.method = "average") - k) / (n - 2 * k + 1))
}

phen <- read.table(snakemake@input[["phenotype"]], header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
qcov <- read.table(snakemake@input[["qcovar"]], header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
trait <- colnames(phen)[3]
transform <- snakemake@params[["transform"]]

# qcovar samples are the realized GRM/PCA set (phenotype is a superset because plink
# --keep intersects with the bfile). Use qcovar as the canonical sample set and align
# the phenotype to it, mirroring how gelex intersects samples before transforming.
key_p <- paste(phen$FID, phen$IID, sep = ":")
key_q <- paste(qcov$FID, qcov$IID, sep = ":")
pidx <- match(key_q, key_p)
keep <- !is.na(pidx) & !is.na(phen[[trait]][pidx])
if (!any(keep)) stop("no overlapping non-missing samples between phenotype and qcovar")

fid <- qcov$FID[keep]
iid <- qcov$IID[keep]
y <- phen[[trait]][pidx[keep]]
pcs <- as.matrix(qcov[keep, -(1:2), drop = FALSE])

if (transform == "iint") {
  y <- rint(residuals(lm(y ~ pcs)))
} else if (transform == "dint") {
  y <- rint(y)
} else if (transform != "none") {
  stop(sprintf("unknown transform: %s", transform))
}

out <- data.frame(FID = fid, IID = iid, check.names = FALSE)
out[[trait]] <- y
write.table(out, snakemake@output[["phenotype"]], sep = "\t",
            quote = FALSE, row.names = FALSE)
