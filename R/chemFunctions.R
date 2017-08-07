
library(stringr)

## Return a table of atomcounts from a chemical formula
atomCounts <- function(form) {
  ## Split the formula string
  f.split <- str_match_all(form, '([A-Z][a-z]*)(\\d*)')
  df.split <- as.data.frame(f.split[[1]], stringsAsFactors = FALSE)
  names(df.split) <- c('Symbol', 'Element', 'Count')
  df.split[[3]] <- as.numeric(df.split[[3]])
  df.split[is.na(df.split)] <- 1
  return(df.split[2:3])
}
