#===============================================================================
# Session information
#===============================================================================

library(here)

out_file <- here("outputs", "session_info.txt")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

writeLines(capture.output(sessionInfo()), out_file)