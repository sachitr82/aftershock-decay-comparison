#' Read and clean a raw SCEDC/SCSN catalogue 
#'
#' Reads a raw catalogue file from the SCEDC catalgoue search ("SCEDC" output
#' format), discards non-data lilnes, names columns, and returns a tidy data
#' frame with parsed times and convenience columns for exploratory analysis.
#' 
#' @details
#' The raw file contains multiple non-data lines: a title, an HTML <PRE> tag, 
#' a DOI line and a closing </PRE> tag. Genuine data rows begin with a date 
#' (YYYY/MM/DD). Columns follow the SCEC_DC format specification:
#' (<https://scedc.caltech.edu/eq-catalogs/docs/scec_dc.html>):
#'   date time et gt mag magtype lat lon depth q evid nph ngrm
#' Event times are UTC.
#' 
#' `datetime` is kept in UTC, while local-time convenience columns are derived
#' in the timezone `tz`.
#' 
#' @param path Path to the raw SCEDC catalogue file
#' @param tz Timezone used to derive the local-time convenience columns;
#'  defaults to '"America/Los_Angeles"'
#' @param col_names Character vector of column names in file order; defaults to
#'   the standard SCEC_DC columns.
#' @return A data frame with one row per event: the raw columns plus
#'   `datetime` (POSIXct, UTC), `event_num` (time-ordered index), `local`
#'   (`datetime` in `tz`), and the local calendar features `year`, `month`,
#'   `month_start`, `hour`, and `wday` (`month` and `wday` are ordered factors).
#'
#' @importFrom utils read.table
#' @importFrom lubridate with_tz year month floor_date hour wday
#' @export

read_catalogue <- function(path, tz = "America/Los_Angeles",
                           col_names = c("date", "time", "et", "gt", "mag",
                                         "magtype", "lat", "lon", "depth", "q",
                                         "evid", "nph", "ngrm")) {
  stopifnot(file.exists(path))
  
  lines <- readLines(path, warn = FALSE)
  
  # Retain rows with genuine data (those starting with a date), dropping the 
  # title / <PRE> / DOI / <\PRE> lines
  data_lines <- grep("^\\s*\\d{4}/\\d{2}/\\d{2}", lines, value = TRUE)
  if(length(data_lines)==0L){
    stop("No data rows found - check the file path/format.", call. = FALSE)
  }
  
  # Robustly convert filtered text into named data frame (tolerates minor
  # format deviations, such as unexpected column counts)
  cat <- read.table(text = data_lines, header = FALSE, fill = TRUE,
                           stringsAsFactors = FALSE)
  if (ncol(cat) != length(col_names)){
    warning(sprintf("Expected %d columns but read %d - check file format.", 
                      length(col_names), ncol(cat)), call.=FALSE)
  }
  colnames(cat) <- col_names[seq_len(ncol(cat))]
  
  # Parse the UTC origin time (date time field with fractional seconds via %oS)
  cat$datetime <- as.POSIXct(paste(cat$date, cat$time),
                             format = "%Y/%m/%d %H:%M:%OS", tz = "UTC")
  
  # Convert relevant fields to numeric
  num <- intersect(c("mag", "lat", "lon", "depth", "nph", "ngrm"), names(cat))
  cat[num] <- lapply(cat[num], as.numeric)

  # Drop rows that failed to parse (no time or no magnitude)
  ok <- !is.na(cat$datetime) & !is.na(cat$mag)
  if (sum(!ok) > 0L) {
    message(sprintf("Dropped %d rows with unparseable time/magnitude.",
                    sum(!ok)))
  }
  cat <- cat[ok, , drop = FALSE]
  
  # Order by natural time and add convenience columns for the EDA
  # datetime, month, year stays UTC; local-time features are derived in `tz`.
  cat <- cat[order(cat$datetime), , drop = FALSE] 
  cat$event_num <- seq_len(nrow(cat))   
  cat$local <- with_tz(cat$datetime,tz) # UTC
  cat$year      <- year(cat$datetime)   # UTC
  cat$month     <- month(cat$datetime, label = TRUE, abbr = TRUE)  # UTC
  cat$month_start <- floor_date(cat$datetime, "month")  # UTC
  cat$hour      <- hour(cat$local) # local
  cat$wday      <- wday(cat$local, label = TRUE, abbr = FALSE, week_start = 1) #
  rownames(cat) <- NULL
  cat
}