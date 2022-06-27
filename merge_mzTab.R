library(MSnbase)

make.pattern <- function(prefix) {
  ## escape characters that would have special meaning in reg. exp.:
  prefix.esc <- gsub("]", "\\]", gsub("[", "\\[", prefix, fixed=TRUE), fixed=TRUE)
  paste0("^", prefix.esc, "\\[(\\d+)\\](.*)$")
}

get.max.index <- function(values, prefix) {
  pattern <- make.pattern(prefix)
  matches <- grep(pattern, values, value=TRUE)
  if (length(matches) == 0)
    return(NULL)
  index <- max(as.integer(sub(pattern, "\\1", matches)))
  attr(index, "matches") <- matches
  index
}

adjust.indexes <- function(values, prefix, offset) {
  pattern <- make.pattern(prefix)
  ind <- grep(pattern, values)
  if (length(ind) == 0) # no such fields
    return(values)
  ## find numbers we want to replace:
  indexes <- as.integer(sub(pattern, "\\1", values[ind]))
  indexes <- indexes + offset
  ## replace:
  values[ind] <- mapply(function(value, num) {
    sub(pattern, paste0(prefix, "[", num, "]\\2"), value)
  }, values[ind], indexes)
  values
}

select.non.missing <- function(values) {
  if (is.na(values[1]))
    values[2]
  else
    values[1]
}

paste.non.missing <- function(values, sep="|") {
  paste(na.omit(values), collapse=sep)
}

min.or.na <- function(values) {
  min <- suppressWarnings(min(values, na.rm=TRUE))
  if (min == Inf) # all values were NA
    return(NA)
  min
}

combine.columns <- function(merged, field, fun=select.non.missing, ...) {
  fields <- paste0(field, c(".x", ".y"))
  if (!all(fields %in% names(merged)))
    return(merged)
  merged[[fields[1]]] <- apply(merged[, fields], 1, fun, ...)
  names(merged)[names(merged) == fields[1]] <- field
  merged[, setdiff(names(merged), fields[2])]
}


collapse.protein.rows <- function(mztab) {
  prot <- mztab@Proteins
  if (!any(duplicated(prot$accession)))
    return(mztab)
  type.col <- "opt_global_result_type"
  npep.col <- "opt_global_nr_found_peptides"
  decoy.col <- "opt_global_cv_PRIDE:0000303_decoy_hit"
  if ((type.col %in% names(prot)) && (("protein_details" %in% prot[[type.col]]) && (("single_protein" %in% prot[[type.col]]) || ("indistinguishable_protein_group" %in% prot[[type.col]])))) {
    ## collapse rows belonging to the same protein:
    parts <- split(prot, prot$accession)
    parts <- lapply(parts, function(part) {
      if (nrow(part) == 1)
        return(part)
      row.ind <- part[[type.col]] == "protein_details"
      if (npep.col %in% names(part))
        part[!row.ind, npep.col] <- part[row.ind, npep.col]
      if (decoy.col %in% names(part))
        part[!row.ind, decoy.col] <- part[row.ind, decoy.col]
      part[!row.ind, ]
    })
    prot <- do.call(rbind, parts)
    rownames(prot) <- NULL
    mztab@Proteins <- prot
  }
  else
    stop("Unexpected duplicate protein accessions found")
  mztab
}


collapse.peptide.rows <- function(mztab) {
  pep <- mztab@Peptides
  cols <- c("sequence", "accession", "modifications", "charge")
  if (!any(duplicated(pep[, cols])))
    return(mztab)
  id <- apply(pep[, cols], 1, paste, collapse="_")
  parts <- split(pep, id)
  ## fields to collapse via "min":
  min.fields <- intersect(c("opt_global_Posterior_Error_Probability_score", "opt_global_q-value"),
                          names(pep))
  min.fields <- c(min.fields, grep("search_engine_score", names(pep), value=TRUE))
  ## fields to combine via "paste":
  ## TODO: avoid conversion of "opt_global_feature_id" values to floating point numbers
  paste.fields <- intersect(c("spectra_ref", "retention_time", "opt_global_feature_id"), names(pep))
  ## number of study variables:
  n.vars <- get.max.index(names(pep), "peptide_abundance_study_variable")
  parts <- lapply(parts, function(part) {
    if (nrow(part) == 1)
      return(part)
    part[1, "mass_to_charge"] <- median(part$mass_to_charge, na.rm=TRUE)
    for (field in min.fields)
      part[1, field] <- min.or.na(part[[field]])
    for (field in paste.fields)
      part[1, field] <- paste(part[[field]], collapse="|")
    ## set other fields depending on row with highest peptide abundance:
    for (i in 1:n.vars) {
      ab.var <- paste0("peptide_abundance_study_variable[", i, "]")
      ind <- which.max(pep[[ab.var]])
      if ((length(ind) == 0) || (ind == 1)) # all NA, or first row is already correct
        next
      pep[1, ab.var] <- pep[ind, ab.var]
      for (var in c("peptide_abundance_stdev", "peptide_abundance_std_error", "opt_global_mass_to_charge", "opt_global_retention_time")) {
        field <- paste0(var, "_study_variable[", i, "]")
        if (field %in% names(pep))
          pep[1, field] <- pep[ind, field]
      }
    }
    part[1, ]
  })
  pep <- do.call(rbind, parts)
  rownames(pep) <- NULL
  mztab@Peptides <- pep
  mztab
}


move.psm.columns <- function(mztab) {
  if (!.hasSlot(mztab, "PSMs"))
    return(mztab)
  psms <- mztab@PSMs
  parts <- split(psms, list(psms$sequence, psms$modifications, psms$charge), drop=TRUE)
  ## columns to remove from PSM section:
  rm.cols <- c("accession", "unique", "database", "database_version", "calc_mass_to_charge", "pre", "post", "start", "end", "opt_global_cv_MS:1002217_decoy_peptide")
  sel.cols <- c("sequence", "modifications", "charge", rm.cols)
  if (.hasSlot(mztab, "Peptides")) { # don't care about the scores on PSM level
    parts <- lapply(parts, function(part) part[1, sel.cols])
  } else { # select row with best score per modified sequence:
    sel.cols <- c(sel.cols, "search_engine_score[1]")
    parts <- lapply(parts, function(part) part[which.min(part[["search_engine_score[1]"]]), sel.cols])
  }
  part <- do.call(rbind, parts)
  if (.hasSlot(mztab, "Peptides")) { # merge columns into existing PEP section
    pep <- mztab@Peptides
    names(pep)[names(pep) == "accession"] <- "opt_global_primary_accession"
    pep <- merge(pep, part, by=c("sequence", "modifications", "charge"), all.x=TRUE)
    for (col in c("unique", "database", "database_version", "opt_global_cv_MS:1002217_decoy_peptide"))
      pep <- combine.columns(pep, col)
    ## move "opt_..." columns to the end:
    opt.cols <- grepl("^opt_", names(pep))
    pep <- pep[, c(which(!opt.cols), which(opt.cols))]
    mztab@Peptides <- pep
  } else { # create new PEP section
    names(part)[length(names(part))] <- "best_search_engine_score[1]"
    mztab@Peptides <- part
  }
  mztab@PSMs <- psms[, setdiff(names(psms), rm.cols)]
  mztab
}


annotate.ms.runs <- function(mztab, n.ms2, paths=NULL) {
  if (!is.null(paths))
    stopifnot(length(n.ms2) == length(paths))
  n.runs <- get.max.index(names(mztab@Metadata), "ms_run")
  if (length(n.ms2) != n.runs)
    stop("Number of 'ms_run[...]' entries doesn't match number of MS2 counts")
  for (i in 1:length(n.ms2)) {
    run <- paste0("ms_run[", i, "]")
    if (!is.null(paths)) {
      mztab@Metadata[[paste0(run, "-location")]] <- paste0("file://", paths[i])
    }
    new.entry <- list(n.ms2[i])
    names(new.entry) <- paste0(run, "-num_ms2")
    mztab@Metadata <- c(mztab@Metadata, new.entry)
  }
  mztab
}


merge.mztab.metadata <- function(meta1, meta2, variable.fields=c("ms_run", "assay", "study_variable")) {
  ## expect that data was analysed in the same way -> same meta data:
  meta.fields <- union(names(meta1), names(meta2))
  ## these entries should differ:
  pattern <- paste0("^", paste0("(", variable.fields, ")", collapse="|"), "\\[\\d+\\]-")
  var.fields <- grep(pattern, meta.fields, value=TRUE)
  ## sanity check:
  for (field in setdiff(meta.fields, var.fields)) {
    if (is.null(meta1[[field]]) || is.null(meta2[[field]])) {
      warning("Metadata entry for field '", field, "' is missing")
    } else if (meta1[[field]] != meta2[[field]]) {
      warning("Metadata entries for field '", field, "' differ")
    }
  }
  merged <- meta1[!(names(meta1) %in% var.fields)]
  ## for "variable" fields, we need to update index numbers in `meta2`:
  ## (e.g. if `meta1` has `ms_run[1]` to `ms_run[5]`, `meta2` needs to start at `ms_run[6]`)
  ## 1. update references inside metadata entries:
  for (field in c("ms_run", "assay")) {
    offset <- get.max.index(unlist(meta1), field)
    adjusted <- adjust.indexes(unlist(meta2), field, offset)
    meta2 <- relist(adjusted, meta2)
  }
  ## 2. update metadata names:
  for (field in variable.fields) {
    ## find max. index for this set of fields in `meta1`:
    ## TODO: use functions `get.max.index` and `adjust.indexes` here?
    pattern <- paste0("^", field, "\\[(\\d+)\\]-.*$")
    fields <- grep(pattern, names(meta1), value=TRUE)
    merged <- c(merged, meta1[fields])
    ind <- grep(pattern, names(meta2))
    if (length(ind) == 0) # no corresponding fields in `meta2`
      next
    offset <- max(as.integer(sub(pattern, "\\1", fields)))
    ## find numbers we want to replace:
    indexes <- as.integer(sub(pattern, "\\1", names(meta2)[ind]))
    indexes <- indexes + offset
    ## replace:
    names(meta2)[ind] <- mapply(function(field, num) {
      sub("\\d+", num, field)
    }, names(meta2)[ind], indexes)
    merged <- c(merged, meta2[ind])
  }
  return(merged)
}


merge.mztab.proteins <- function(prot1, prot2,
                                 variable.fields=c("protein_abundance_assay", "protein_abundance_study_variable", "protein_abundance_stdev_study_variable", "protein_abundance_std_error_study_variable")) {
  for (field in variable.fields) {
    offset <- get.max.index(names(prot1), field)
    if (!is.null(offset)) {
      names(prot2) <- adjust.indexes(names(prot2), field, offset)
    }
  }
  merged <- merge(prot1, prot2, by="accession", all=TRUE)
  ## TODO: support more than one protein score?
  merged <- combine.columns(merged, "best_search_engine_score[1]", min, na.rm=TRUE)
  ## TODO: recalculate values in the columns below? or keep multiple columns ("...[1]" etc.)?
  merged <- combine.columns(merged, "protein_coverage", max, na.rm=TRUE)
  merged <- combine.columns(merged, "opt_global_nr_found_peptides", max, na.rm=TRUE)
  merged <- combine.columns(merged, "opt_global_Posterior_Probability_score", min.or.na)
  merged <- combine.columns(merged, "opt_global_result_type", function(row) {
    if ("indistinguishable_protein_group" %in% row)
      return("indistinguishable_protein_group")
    if ("single_protein" %in% row)
      return("single_protein")
    "protein_details"
  })
  ## other columns to combine:
  cols1 <- grep("\\.x$", names(merged), value=TRUE)
  cols1 <- sub("\\.x$", "", cols1)
  cols2 <- grep("\\.y$", names(merged), value=TRUE)
  cols2 <- sub("\\.y$", "", cols2)
  for (col in intersect(cols1, cols2)) {
    merged <- combine.columns(merged, col)
  }
  merged
}


merge.mztab.peptides <- function(pep1, pep2, variable.fields=c("search_engine_score[1]_ms_run", "peptide_abundance_study_variable", "peptide_abundance_stdev_study_variable", "peptide_abundance_std_error_study_variable", "opt_global_mass_to_charge_study_variable", "opt_global_retention_time_study_variable")) {
  for (field in variable.fields) {
    offset <- get.max.index(names(pep1), field)
    if (!is.null(offset)) {
      names(pep2) <- adjust.indexes(names(pep2), field, offset)
    }
  }
  ## adjust "ms_run" indexes in "spectra_ref" column:
  refs <- strsplit(pep1$spectra_ref, "|", fixed=TRUE)
  offset <- get.max.index(unlist(refs), "ms_run")
  refs <- strsplit(pep2$spectra_ref, "|", fixed=TRUE)
  adjusted <- adjust.indexes(unlist(refs), "ms_run", offset)
  pep2$spectra_ref <- sapply(relist(adjusted, refs), paste, collapse="|")

  merged <- merge(pep1, pep2, by=c("sequence", "modifications", "charge"), all=TRUE)
  ## TODO: support more than one protein score?
  for (field in c("best_search_engine_score[1]", "opt_global_Posterior_Error_Probability_score", "opt_global_q-value")) {
  merged <- combine.columns(merged, field, min, na.rm=TRUE)
  }
  for (field in c("spectra_ref", "retention_time", "opt_global_feature_id")) {
    merged <- combine.columns(merged, field, paste.non.missing)
  }
  fields <- grep("^opt_global_mass_to_charge_study_variable", names(merged))
  if (length(fields) > 0) { # recalculate overall mass-to-charge
    merged$mass_to_charge.x <- rowMeans(merged[, fields], na.rm=TRUE)
  }
  ## other columns to combine:
  cols1 <- grep("\\.x$", names(merged), value=TRUE)
  cols1 <- sub("\\.x$", "", cols1)
  cols2 <- grep("\\.y$", names(merged), value=TRUE)
  cols2 <- sub("\\.y$", "", cols2)
  for (col in intersect(cols1, cols2)) {
    merged <- combine.columns(merged, col)
  }
  ## move "opt_..." columns to the end:
  opt.cols <- grep("^opt_", names(merged), value=TRUE)
  merged <- merged[, c(setdiff(names(merged), opt.cols), opt.cols)]
  merged
}


merge.mztab.psms <- function(psms1, psms2) {
  ## adjust "PSM_ID":
  psms2$PSM_ID <- psms2$PSM_ID + max(psms1$PSM_ID) + 1 # index starts at 0
  ## adjust "ms_run" index in "spectra_ref" column:
  offset <- get.max.index(psms1$spectra_ref, "ms_run")
  psms2$spectra_ref <- adjust.indexes(psms2$spectra_ref, "ms_run", offset)
  ## adjust index in "opt_global_map_index" column:
  ## ("..._map_index" is NA if "..._feature_id" is "not mapped")
  col <- "opt_global_map_index"
  if ((col %in% names(psms1)) && (col %in% names(psms2)))
    psms2[[col]] <- psms2[[col]] + max(psms1[[col]], na.rm=TRUE) + 1 # index starts at 0
  rbind(psms1, psms2)
}



## Merge data from `mztab2` into `mztab1`
merge.mztab.data <- function(mztab1, mztab2) {
  mztab1@Metadata <- merge.mztab.metadata(mztab1@Metadata, mztab2@Metadata)
  mztab1@Proteins <- merge.mztab.proteins(mztab1@Proteins, mztab2@Proteins)
  mztab1@Peptides <- merge.mztab.peptides(mztab1@Peptides, mztab2@Peptides)
  mztab1@PSMs <- merge.mztab.psms(mztab1@PSMs, mztab2@PSMs)
  mztab1@Filename <- paste(mztab1@Filename, mztab2@Filename, sep="|")
  mztab1
}


load.mztab.file <- function(path) {
  mztab <- MzTab(path)
  mztab <- collapse.protein.rows(mztab)
  if (.hasSlot(mztab, "Peptides"))
    mztab <- collapse.peptide.rows(mztab)
  mztab <- move.psm.columns(mztab)
  mztab
}


merge.mztab.files <- function(paths) {
  if (length(paths) == 0)
    stop("list of input paths is empty")
  cat("Working on file", paths[1], "\n")
  merged <- load.mztab.file(paths[1])
  for (path in paths[-1]) {
    cat("Working on file", path, "\n")
    mztab <- load.mztab.file(path)
    merged <- merge.mztab.data(merged, mztab)
  }
  cat("Done.\n")
  merged
}


test <- function() {
  mztab.dir <- "R:/Bioinformatics/Proteomics/SST/QuantMS/single_runs"
  paths <- list.files(mztab.dir, "\\.mzTab$", full.names=TRUE)
  merged <- merge.mztab.files(paths[1:3])
}
