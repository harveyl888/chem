#' Count element occurrence
#'
#' Return a table of atomcounts from a chemical formula
#'
#' @param form Input formula
#'
#' @return data frame consisting of molecular symbol, element and count
#'
#' @import stringr
atomCounts <- function(form) {
  ## Split the formula string
  f.split <- str_match_all(form, '([A-Z][a-z]*)(\\d*)')
  df.split <- as.data.frame(f.split[[1]], stringsAsFactors = FALSE)
  names(df.split) <- c('Symbol', 'Element', 'Count')
  df.split[[3]] <- as.numeric(df.split[[3]])
  df.split[is.na(df.split)] <- 1
  return(df.split[2:3])
}

#' Split chemical formula
#'
#' Split a formula into chunks and return atom counts
#'
#' @param form Input formula
#'
#' @import stringr
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


#' Return a molecular formula
#'
#' Return a standardized molecular formula
#'
#' @param form Input formula
#' @param output output type - either table (columns = Element and Count) or formula (string)
#' @param order atom ordering - either alpha or hill notation
#'
#' @return a tibble or string corresponding to the molecular formula
#'
#' @import dplyr
#' @export
getMolecularFormula <- function(form, output = 'table', order = 'hill') {
  split.formula <- formulaSplit(form)
  fAll <- lapply(split.formula, function(x) {
    df <- atomCounts(x[2])
    df[[2]] <- df[[2]] * as.numeric(x[1])
    df
  })
  df.form <- bind_rows(fAll)
  df.form <- df.form %>%
    group_by(Element) %>%
    summarise(Count = sum(Count))
  hillOrder <- c('C', 'H', df.form[['Element']][!df.form[['Element']] %in% c('C', 'H')])
  if (order == 'hill') {
    df.form <- df.form %>% slice(match(hillOrder, Element))
  } else {
    df.form <- df.form %>% arrange(Element)
  }
  if (output == 'table') {
    return(df.form)
  } else {
    return(paste0(paste0(df.form[[1]], ifelse(df.form[[2]] > 1, df.form[[2]], ''), collapse = '')))
  }
}


#' Calculate exact mass
#'
#' Calculate exact mass from formula
#'
#' @param form Input formula
#'
#' @return exact mass (numeric)
#'
#' @export
exactMass <- function(form) {
  df.form <- getMolecularFormula(form, output = 'table', order = 'alpha')
  atomMissing <- which(!df.form[[1]] %in% df.mass$atom)
  if (length(atomMissing) > 0) return (paste0('atoms not recognized: ', paste0(df.form[atomMissing, 1], collapse=', ')))
  exact_mass <- sum(apply(df.form, 1, function(x) as.numeric(x[2]) * df.mass[df.mass$atom == x[1], ]$mass))
  return(exact_mass)
}
