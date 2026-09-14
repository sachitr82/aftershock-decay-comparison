# Gutenberg-Richter b-value estimators and completeness (Mc) tools.
# `dm` is the analysis bin / cutoff spacing used by completeness methods.
# `dmag` is the catalogue magnitude discretisation (0.01 for the SCEDC data).

#' Frequency-magnitude (Gutenberg-Richter) counts
#'
#' Bins magnitudes and returns the incremental (per-bin) and cumulative
#' (\eqn{N \ge M}) counts used for GR / frequency-magnitude plots.
#'
#' @param mags Numeric vector of magnitudes.
#' @param dm Magnitude bin width (default 0.1).
#' @return A data frame with columns `mag` (bin centre), `incremental` (count
#'   in bin) and `cumulative` (cumulative count from that magnitude bin upward).

gr_data <- function(mags, dm = 0.1) {
  mags <- mags[is.finite(mags)]
  brks <- seq(floor(min(mags) / dm) * dm - dm / 2, max(mags) + dm, by = dm)
  h <- hist(mags, breaks = brks, plot = FALSE)
  data.frame(mag = h$mids,
             incremental = h$counts,
             cumulative  = rev(cumsum(rev(h$counts))))
}

#' Maximum-curvature estimate of the completeness magnitude (Mc)
#'
#' Maximum-curvature (MAXC) method: the modal bin of the non-cumulative
#' frequency-magnitude distribution.
#'
#' @param mags Numeric vector of magnitudes.
#' @param dm Magnitude bin width (default 0.1).
#' @return The estimated magnitude of completeness (numeric scalar).
#' @seealso [gr_data()]

estimate_mc_maxc <- function(mags, dm = 0.1) {
  g <- gr_data(mags, dm)
  g$mag[which.max(g$incremental)]
}

#' Aki-Utsu maximum-likelihood b-value (binning-corrected)
#'
#' Maximum-likelihood Gutenberg-Richter b-value above a completeness magnitude
#' (Aki, 1965; Utsu, 1965), with the  Bender (1983) half-bin correction applied
#' by default, and the Shi & Bolt (1982) standard error.
#'
#' @param mags Numeric vector of magnitudes.
#' @param mc Magnitude cut-off; events in the `mc` magnitude class and above
#'   are used.
#' @param dmag Catalogue magnitude discretisation (default 0.01).
#' @param correction Half-bin correction used in the denominator
#'   (default `dmag/2`).
#' @return Named numeric vector: `b`, `b_err` (Shi & Bolt error) and `n` (events
#'   used). `b`/`b_err` are `NA` if fewer than two events are retained.
#' @examples b_aki(c(3.2, 4.1, 0.5), mc = 3.0)

b_aki <- function(mags, mc, dmag = 0.01, correction = dmag / 2) {
  m <- mags[is.finite(mags) & mags > mc - dmag / 2]
  n <- length(m)
  if (n < 2) return(c(b = NA_real_, b_err = NA_real_, n = n))
  mbar <- mean(m)
  b <- log10(exp(1)) / (mbar - (mc - correction))
  b_err <- log(10) * b^2 * sqrt(sum((m - mbar)^2) / (n * (n - 1)))
  c(b = b, b_err = b_err, n = n)
}

#' b-value stability curve (MBS)
#'
#' b-value as a function of cut-off magnitude for the "Mc by b-value stability"
#' method (Cao & Gao 2002; Woessner & Wiemer 2005).
#'  Uses the binning-corrected Aki MLE ([b_aki()]).
#'
#' @param mags Numeric vector of magnitudes.
#' @param mc_range Optional cut-offs. Defaults to min(mags) .. max-5*dm by `dm`.
#' @param dm Candidate-Mc spacing (default 0.1).
#' @param dmag Catalogue magnitude discretisation (default 0.01).
#' @param correction Half-bin correction passed to [b_aki()] (default `dmag/2`).
#' @return Data frame, one row per cut-off: `mc`, `b` (Aki/Utsu MLE), `b_err`
#'   (Shi & Bolt) and `n`.
#' @seealso [b_aki()], [mc_mbs()]

b_stability <- function(mags, mc_range = NULL, dm = 0.1,
                        dmag = 0.01, correction = dmag / 2) {
  mags <- mags[is.finite(mags)]
  if (is.null(mc_range)) {
    mc_range <- seq(floor(min(mags) / dm) * dm, max(mags) - 5 * dm, by = dm)
  }
  do.call(rbind, lapply(mc_range, function(mc) {
    a <- b_aki(mags, mc, dmag = dmag, correction = correction)
    data.frame(mc = mc, b = a[["b"]], b_err = a[["b_err"]], n = a[["n"]])
  }))
}

#' Mc by b-value stability (Woessner & Wiemer 2005 criterion)
#'
#' Returns the lowest cut-off at which the local b-value agrees with the mean of
#' the b-values over the next half-magnitude window (`dM`), to within its Shi &
#' Bolt uncertainty: Mc = min{Mco : |bave - b| <= b_err}. `bave` is the mean of
#' the `dM/dm + 1` cut-offs from Mco through Mco+dM (6 for dM = 0.5, dm = 0.1).
#'
#' @param mags Numeric vector of magnitudes.
#' @param dm Candidate-Mc spacing (default 0.1).
#' @param dmag Catalogue magnitude discretisation (default 0.01).
#' @param dM Averaging window for `bave` (default 0.5).
#' @return List with `mc` (NA if the criterion is never met) and `curve`
#'   (the [b_stability()] data frame augmented with `bave` and `dbdiff`).
#' @seealso [b_stability()]

mc_mbs <- function(mags, dm = 0.1, dmag = 0.01, dM = 0.5) {
  cur <- b_stability(mags, dm = dm, dmag = dmag)
  w <- round(dM / dm) + 1                                 
  n <- nrow(cur)
  cur$bave <- vapply(seq_len(n), function(i) {
    if (i + w - 1 > n) return(NA_real_)               
    mean(cur$b[i:(i + w - 1)])            
  }, numeric(1))
  cur$dbdiff <- abs(cur$bave - cur$b)
  ok <- which(cur$dbdiff <= cur$b_err)
  list(mc = if (length(ok)) cur$mc[min(ok)] else NA_real_, curve = cur)
}

#' Goodness-of-fit completeness estimate (GFT)
#'
#' Wiemer & Wyss (2000): Mc is the lowest cut-off whose synthetic-vs-observed
#' fit R first reaches `level`% (not where R is maximal).
#'
#' @param mags Numeric vector of magnitudes.
#' @param dm Candidate-Mc / FMD evaluation spacing (default 0.1).
#' @param dmag Catalogue magnitude discretisation (default 0.01).
#' @param level Fit threshold in percent (typically 90 or 95).
#' @param min_n Minimum events above a cut-off required to attempt a fit.
#' @return List: `mc` (NA if `level` never reached), `cutoffs`, `R`, `level`.
#' @seealso [b_aki()]

mc_gft <- function(mags, dm = 0.1, dmag = 0.01, level = 90, min_n = 50)  {
  mags <- mags[is.finite(mags)]
  cutoffs <- seq(floor(min(mags) / dm) * dm, max(mags) - 10 * dm, by = dm)
  R <- rep(NA_real_, length(cutoffs))
  for (i in seq_along(cutoffs)) {
    mco <- cutoffs[i]
    ms <- mags[mags > mco - dmag / 2]
    if (length(ms) < min_n) next
    b <- b_aki(ms, mco, dmag = dmag)[["b"]]           
    a <- log10(length(ms)) + b * mco
    bins <- seq(mco, max(ms), by = dm)
    B <- vapply(bins, function(mm) sum(ms > mm - dmag / 2), numeric(1))
    S <- 10^(a - b * bins)
    R[i] <- 100 * (1 - sum(abs(B - S)) / sum(B))
  }
  ok <- which(R >= level)
  list(mc = if (length(ok)) cutoffs[min(ok)] else NA_real_,
       cutoffs = cutoffs, R = R, level = level)
}