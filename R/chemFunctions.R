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
getMolecularFormula <- function(form = '', output = 'table', order = 'hill') {

  ## initial check - NA, NULL or empty string
  if (is.null(form)) return ('no formula')
  if (is.na(form) | nchar(form) == 0 | form == '.' | form == 'no formula') return ('no formula')

  ## check for polymer
  if (grepl('\\)n', form) | form == 'polymer') return ('polymer')

  ## check for R group
  if (grepl('R(?![anehgbuf])', form, perl = TRUE) | form == 'R group') return ('R group')

  ## check for X.Y formula
  v.formulae <- str_trim(unlist(str_split(form, '\\.')))
  v.multiplier <- str_extract(v.formulae, '^\\d+')

  ## split formula for each section
  fAll <- lapply(seq(v.formulae), function(x) {
    split.formula <- formulaSplit(v.formulae[x])
    df.out <- lapply(split.formula, function(y) {
      df <- atomCounts(y[2])
      df[[2]] <- df[[2]] * as.numeric(y[1]) * ifelse(is.na(v.multiplier[x]), 1, as.numeric(v.multiplier[x]))
      df
    })
    bind_rows(df.out)
  })

  df.form <- bind_rows(fAll)
  df.form <- df.form %>%
    group_by(Element) %>%
    summarise(Count = sum(Count)) %>%
    arrange(Element)

  if ('C' %in% df.form[['Element']]) {
    hillOrder <- c('C', 'H', df.form[['Element']][!df.form[['Element']] %in% c('C', 'H')])
  } else {
    hillOrder <- df.form[['Element']]
  }
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


#' Check cas number is correct
#'
#' Check cas number is correct and return true if valid or false if not
#'
#' @param cas Input cas
#'
#' @return boolean
#'
#' @export
cas_check <- function(cas) {
  cas_split <- str_split(cas, '-')[[1]]
  cas_validation <- as.numeric(cas_split[3])
  cas_to_check <- as.numeric(rev(strsplit(paste0(cas_split[1], cas_split[2]), '')[[1]]))
  cas_sum <- sum(cas_to_check * seq(length(cas_to_check)))
  return(cas_sum %% 10 == cas_validation)
}

#' return smiles string from cas
#'
#' return a smiles string from a cas identifier using the NCI/CADD cactus service.
#'     If no smiles string is availble then NA is returned.
#'
#' @param cas Input cas
#'
#' @return string
#'
#' @export
cas_to_smiles <- function(cas) {
  url <- sprintf('https://cactus.nci.nih.gov/chemical/structure/%s/smiles', cas)
  tryCatch({
    smiles <- readLines(url, warn = FALSE)
    smiles[1]
  },
  error = function(e) {
    NA
  })
}


#' Extact data from PubChem
#'
#' Extact data from PubChem using the cas number as a lookup.
#'     If more than one CID are matched, data for the highest ranked one are returned
#'
#' @param cas Input cas
#'
#' @return named list of properties
#'
#' @import jsonlite
#' @import dplyr
#'
#' @export
#'
data_from_pubchem <- function(cas) {
  url_txt <- sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/xref/RN/%s/cids/json', cas)
  top_cid <- tryCatch({
    grab_cids <- jsonlite::fromJSON(url(url_txt), flatten = T)
    df_cids <- grab_cids$InformationList$Information
    df_cids$CID <- sapply(df_cids$CID, function(x) if_else(is.null(x), NA_integer_, x))
    df_cids %>%
      dplyr::filter(!is.na(CID)) %>%
      dplyr::group_by(CID) %>%
      dplyr::summarise(count=n()) %>%
      dplyr::arrange(desc(count)) %>%
      dplyr::slice(1) %>%
      dplyr::select(CID) %>%
      dplyr::pull()
  },
  error = function(e) {
    NA
  })
  if (!is.na(top_cid)) {
    url_txt <- sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/CanonicalSMILES,IsomericSMILES,InChI,IUPACName,MolecularFormula/json', top_cid)
    pubchem_out <- jsonlite::fromJSON(url(url_txt))
    l.pubchem <- pubchem_out$PropertyTable$Properties
    l.pubchem[['cas']] <- cas
    l.pubchem
  } else {
    list(cas = cas)
  }
}
