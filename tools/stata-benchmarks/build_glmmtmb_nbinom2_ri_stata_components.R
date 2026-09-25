# Run from tools/stata-benchmarks with the returned revision-1 log as an argument.
args <- commandArgs(trailingOnly = TRUE)
log_file <- if (length(args)) args[1] else "suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark.log"
lines <- readLines(log_file, warn = FALSE)
starts <- grep("^=== LAPLACE COMPONENT: ", lines)
outcomes <- c("y1", "y2", "y1_partial", "y2_partial", "y1_left")
trailing_numbers <- function(s) {
  tokens <- strsplit(trimws(s), "[[:space:]]+")[[1L]]
  if (length(tokens) < 2) return(numeric())
  values <- suppressWarnings(as.numeric(tokens[-1L]))
  if (anyNA(values)) numeric() else values
}
models <- list()
for (j in seq_along(starts)) {
  end <- if (j < length(starts)) starts[j + 1L] - 1L else length(lines)
  block <- lines[starts[j]:end]
  outcome <- sub("^=== LAPLACE COMPONENT: (.*) ===$", "\\1", block[1])
  if (!outcome %in% outcomes) next
  ib <- grep("^e\\(b\\)\\[1,5\\]", block)
  iv <- grep("^symmetric e\\(V\\)\\[5,5\\]", block)
  b <- unlist(lapply(grep("^y1 +", block[ib:iv], value = TRUE), trailing_numbers))
  rows <- lapply(block[(iv+1):length(block)], trailing_numbers)
  rows <- rows[lengths(rows)>0][1:5]
  stopifnot(length(b)==5, identical(lengths(rows), 1:5))
  V <- matrix(0,5,5)
  for (i in 1:5) V[i,1:i] <- rows[[i]]
  V <- V + t(V) - diag(diag(V))
  J <- diag(c(1,1,1,-1,1/(2*b[5])))
  V <- J %*% V %*% J
  b[4] <- -b[4]
  b[5] <- .5*log(b[5])
  order <- c(3,1,2,4,5)
  b <- b[order]; V <- V[order,order]
  names(b) <- c("(Intercept)", "x", "z", "log_phi", "log_sigma")
  dimnames(V) <- list(names(b),names(b))
  ll <- as.numeric(sub(".*LL= *", "", grep("^N=.*LL=", block, value=TRUE)))
  models[[outcome]] <- list(coefficients=b,vcov=V,logLik=ll)
  if (outcome %in% c("y1","y2")) {
    ibm <- grep("^symmetric r\\(b\\)",block)
    ivm <- grep("^symmetric r\\(V\\)",block)
    bm <- vapply(ibm,function(k) as.numeric(tail(strsplit(trimws(block[k+2])," +")[[1]],1)),numeric(1))
    vm <- vapply(ivm,function(k) as.numeric(tail(strsplit(trimws(block[k+2])," +")[[1]],1)),numeric(1))
    models[[outcome]]$native_margins <- list(estimates=bm,variances=vm)

  }
}
dat <- read.csv("suest_r_glmmtmb_nbinom2_ri_crosslang_benchmark.csv")
saveRDS(list(data=dat, models=models, metadata=list(Stata="19.5",revision=1L,
  source_log=basename(log_file),boundary_case="y2_right; Stata r(430)")),
  "../../tests/testthat/fixtures/glmmtmb-nbinom2-ri-stata-components.rds")
right <- dat[dat$id>40,c("id","time","x","z","y2")]
right$id <- factor(right$id)
saveRDS(right,"../../tests/testthat/fixtures/glmmtmb-nbinom2-ri-boundary.rds")
cat("NBINOM2_REVISION1_COMPONENT_AND_BOUNDARY_FIXTURES_COMPLETE=1\n")
