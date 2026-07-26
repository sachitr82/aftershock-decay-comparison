# Utility functions reused across the EDA notebook: frequency-magnitude
# tables, an initial estimate of the completeness magnitude (Mc), sequence
# extraction, time-since-mainshock, and tracking Mc over time.

#' Frequency-magnitude (Gutenberg-Richter) counts
#'
#' Bins magnitudes and returns the incremental (per-bin) and cumulative
#' (\eqn{N \ge M}) counts used for GR / frequency-magnitude plots.
#'
#' @param mags Numeric vector of magnitudes.
#' @param dm Magnitude bin width (default 0.1).
#' @return A data frame with columns `mag` (bin centre), `incremental` (count
#'   in bin) and `cumulative` (count of events with magnitude >= `mag`).
#' @export
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
#' @export
estimate_mc_maxc <- function(mags, dm = 0.1) {
  g <- gr_data(mags, dm)
  g$mag[which.max(g$incremental)]
}

#' Extract a space-time-magnitude sub-catalogue (the sequence)
#'
#' Subsets a catalogue to a spatial box, an optional time window and an
#' optional minimum magnitude, then re-numbers the retained events.
#'
#' @param cat A catalogue data frame with `lon`, `lat`, `datetime` and `mag`.
#' @param lon_range,lat_range Length-2 numeric vectors giving the box bounds.
#' @param t_range Optional length-2 POSIXct vector `c(start, end)`.
#' @param mmin Optional minimum magnitude.
#' @return The subset data frame, with `event_num` re-numbered from 1 and row
#'   names reset.
#' @export
extract_sequence <- function(cat, lon_range, lat_range, t_range = NULL,
                             mmin = NULL) {
  keep <- cat$lon >= lon_range[1] & cat$lon <= lon_range[2] &
    cat$lat >= lat_range[1] & cat$lat <= lat_range[2]
  if (!is.null(t_range)) {
    keep <- keep & cat$datetime >= t_range[1] & cat$datetime <= t_range[2]
  }
  if (!is.null(mmin)) {
    keep <- keep & cat$mag >= mmin
    }
  out <- cat[keep, , drop = FALSE]
  out$event_num <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

#' Days elapsed since a reference time
#'
#' @param datetime A POSIXct vector.
#' @param t0 A reference POSIXct (e.g. the mainshock origin time).
#' @return Numeric vector of signed differences in days.
#' @export
days_since <- function(datetime, t0) {
  as.numeric(difftime(datetime, t0, units = "days"))
}

#' Windowed completeness magnitude to track short-term incompleteness (STAI)
#'
#' Estimates Mc (maximum curvature) in successive blocks of `block` events,
#' revealing the elevated completeness magnitude immediately after a mainshock.
#'
#' @param cat A catalogue data frame with `datetime` and `mag`.
#' @param block Number of events per block (default 250).
#' @param dm Magnitude bin width (default 0.1).
#' @return A data frame with `t_mid` (block mid-time), `mc`, and `n`
#'   (events per block); empty if fewer than `block` events.
#' @seealso [estimate_mc_maxc()]
#' @export
mc_in_blocks <- function(cat, block = 250, dm = 0.1) {
  cat <- cat[order(cat$datetime), ]
  n <- nrow(cat)
  if (n < block) return(data.frame())
  starts <- seq(1, n - block + 1, by = block)
  do.call(rbind, lapply(starts, function(s) {
    idx <- s:(s + block - 1)
    data.frame(
      t_mid = cat$datetime[s + block %/% 2],
      mc    = estimate_mc_maxc(cat$mag[idx], dm),
      n     = block
    )
  }))
}