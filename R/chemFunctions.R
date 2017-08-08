
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

## Split a formula into chunks and return atom counts
formulaSplit <- function(form) {
  out <- list()
  formulaSplitInternal <- function(form, mult) {
    form.section <- str_extract_all(form, '[A-Z][a-z]?\\d*|\\((?:[^()]*(?:\\(.*\\))?[^()]*)+\\)\\d+')[[1]]
    send.out <- form.section[!grepl('\\(', form.section)]
    send.recurse <- form.section[grepl('\\(', form.section)]
    for (i in send.out) out[[length(out) + 1]] <<- c(mult, i)
    if (length(send.recurse) == 0) {
    } else {
      for (i in send.recurse) {
        newMultAll <- str_extract_all(i, '\\)\\d+')[[1]]
        newMultLast <- newMultAll[length(newMultAll)]
        if (nchar(newMultLast) == 1) {
          m <- mult
        } else {
          m <- mult * as.integer(substring(newMultLast, 2))
        }
        newForm <- substring(i, 2, nchar(i)-nchar(newMultLast))
        Recall(newForm, m)
      }
    }
  }
  formulaSplitInternal(form, 1)
  return(out)
}


library(dplyr)
## Return atom counts
getCount <- function(form, type = 'table') {
  split.formula <- formulaSplit(form)
  fAll <- lapply(split.formula, function(x) {
    df <- atomCounts(x[2])
    df[[2]] <- df[[2]] * as.numeric(x[1])
    df
  })
  df.form <- bind_rows(fAll)
  df.form <- df.form %>%
    group_by(Element) %>%
    summarise(Count = sum(Count)) %>%
    arrange(Element)
  if (type == 'table') {
    return(df.form)
  } else {
    return(paste0(paste0(df.form[[1]], df.form[[2]], collapse = '')))
  }
}
