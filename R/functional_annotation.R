#' Annotate an input table against the four bundled functional databases
#'
#' Matches each input row against the database Species_Name column using the
#' available rank1...rankN values. No random assignment is performed. The input
#' rows and sample-abundance columns are retained and one XLSX and TSV are
#' written for each requested database.
#'
#' @param input_file Path to an XLSX or TSV/TXT file containing ASV/component,
#'   an identifier column and taxonomy/species values in any subsequent columns.
#' @param output_dir Output directory.
#' @param db_names Database names to process.
#' @param db_dir Directory containing the bundled RDS files. If NULL, uses the
#'   installed package's inst/db directory.
#' @param write_xlsx Whether to write Excel files in addition to TSV.
#' @return Invisibly, a named list of output file paths and matching summaries.
#' @export
annotate_functional_databases <- function(
    input_file,
    output_dir = "functional_analysis/annotations",
    db_names = c("copepoda", "mixoplankton", "phytoplankton", "habs"),
    db_dir = NULL,
    write_xlsx = TRUE
) {
    if (!file.exists(input_file)) stop("Input file not found: ", input_file, call. = FALSE)
    ext <- tolower(tools::file_ext(input_file))
    if (ext == "xlsx") {
        input <- readxl::read_excel(input_file)
    } else if (ext %in% c("tsv", "txt", "csv")) {
        delim <- if (ext == "csv") "," else "\t"
        input <- readr::read_delim(input_file, delim = delim, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character()))
        # Restore numeric sample columns after reading all fields as text.
        sample_guess <- plankannot_detect_sample_columns(names(input))
        for (s in sample_guess) input[[s]] <- suppressWarnings(as.numeric(input[[s]]))
    } else {
        stop("Unsupported input format. Use .xlsx, .tsv, .txt or .csv.", call. = FALSE)
    }
    names(input) <- plankannot_safe_utf8(names(input))
    input[] <- lapply(input, function(x) if (is.character(x)) plankannot_safe_utf8(x) else x)

    # Input columns are interpreted by position/content, not by their names.
    # The first column is retained as the ASV/feature identifier. Every cell
    # from the remaining columns is eligible for taxonomy-name matching.
    # Semicolon-delimited taxonomy strings are split into individual taxa.
    sample_cols <- if (ncol(input) >= 3L) names(input)[2:(ncol(input) - 1L)] else character()
    taxonomy_cols <- if (ncol(input) >= 2L) names(input)[-1L] else names(input)
    if (ncol(input) < 2L) {
        stop("Input must contain at least an identifier column and taxonomy/data column.", call. = FALSE)
    }

    if (is.null(db_dir)) db_dir <- plankannot_find_db_dir()
    if (!dir.exists(db_dir)) stop("Database directory not found: ", db_dir, call. = FALSE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stem <- tools::file_path_sans_ext(basename(input_file))
    output_paths <- list()
    summaries <- list()

    for (db_name in db_names) {
        db_path <- file.path(db_dir, paste0(db_name, ".rds"))
        if (!file.exists(db_path)) stop("Database not found: ", db_path, call. = FALSE)
        db <- readRDS(db_path)
        if (!is.data.frame(db)) stop("Database '", db_name, "' is not a data.frame.", call. = FALSE)
        names(db) <- plankannot_safe_utf8(names(db))
        db[] <- lapply(db, function(x) if (is.character(x) || is.factor(x)) plankannot_safe_utf8(as.character(x)) else x)
        if (!"Species_Name" %in% names(db)) stop("Database '", db_name, "' has no Species_Name column.", call. = FALSE)

        norm <- function(x) tolower(trimws(gsub("[[:space:]]+", " ", plankannot_safe_utf8(as.character(x)))))
        db_key <- norm(db$Species_Name)
        keep <- !is.na(db_key) & nzchar(db_key)
        db <- db[keep, , drop = FALSE]
        db_key <- db_key[keep]
        # First occurrence wins if the database contains duplicate species names.
        unique_idx <- !duplicated(db_key)
        db <- db[unique_idx, , drop = FALSE]
        db_key <- db_key[unique_idx]

        # Match species names against values in each input row, irrespective
        # of input column names. Taxonomy strings containing semicolons are
        # split; candidate taxa are checked from right to left so the most
        # specific (usually species-level) name is preferred.
        norm_match <- function(x) {
            x <- tolower(trimws(plankannot_safe_utf8(as.character(x))))
            x <- gsub("[ _]+", " ", x)
            x
        }
        db_key_match <- norm_match(db$Species_Name)

        matched_idx <- rep(NA_integer_, nrow(input))
        matched_rank <- rep(NA_character_, nrow(input))
        matched_taxon <- rep(NA_character_, nrow(input))

        # Scan columns from right to left. This naturally prioritizes a final
        # taxonomy column while also supporting files with taxonomy split
        # across several unnamed or differently named columns.
        for (col_i in rev(seq_along(taxonomy_cols))) {
            col_name <- taxonomy_cols[[col_i]]
            vals <- as.character(input[[col_name]])
            for (row_i in seq_len(nrow(input))) {
                if (!is.na(matched_idx[[row_i]])) next
                cell <- vals[[row_i]]
                if (is.na(cell) || !nzchar(trimws(cell))) next

                candidates <- trimws(unlist(strsplit(cell, ";", fixed = TRUE), use.names = FALSE))
                candidates <- candidates[nzchar(candidates)]
                if (!length(candidates)) next

                # Prefer the last (most specific) taxonomy element.
                candidates <- rev(candidates)
                candidate_keys <- norm_match(candidates)
                hit <- match(candidate_keys, db_key_match)
                hit_pos <- which(!is.na(hit))[1]
                if (length(hit_pos) && !is.na(hit_pos)) {
                    matched_idx[[row_i]] <- hit[[hit_pos]]
                    matched_rank[[row_i]] <- col_name
                    matched_taxon[[row_i]] <- candidates[[hit_pos]]
                }
            }
        }

        extra <- db[matched_idx, setdiff(names(db), "Species_Name"), drop = FALSE]
        names(extra) <- make.unique(paste0(db_name, "_", names(extra)))
        result <- dplyr::bind_cols(input, extra)
        result$Matched <- !is.na(matched_idx)
        result$Matched_taxon <- matched_taxon
        result$Matched_rank <- matched_rank

        db_out <- file.path(output_dir, db_name)
        dir.create(db_out, recursive = TRUE, showWarnings = FALSE)
        tsv_path <- file.path(db_out, paste0(stem, "_", db_name, "_annotated.tsv"))
        readr::write_tsv(plankannot_clean_for_export(result), tsv_path, na = "")
        xlsx_path <- NULL
        if (isTRUE(write_xlsx)) {
            xlsx_path <- file.path(db_out, paste0(stem, "_", db_name, "_annotated.xlsx"))
            writexl::write_xlsx(list(annotated_dataset = plankannot_clean_for_export(result)), xlsx_path)
        }
        summaries[[db_name]] <- data.frame(database = db_name, rows = nrow(result), matched = sum(result$Matched), unmatched = sum(!result$Matched), match_percent = round(100 * mean(result$Matched), 2))
        output_paths[[db_name]] <- c(tsv = tsv_path, xlsx = if (is.null(xlsx_path)) NA_character_ else xlsx_path)
        message("[", db_name, "] matched ", sum(result$Matched), "/", nrow(result), " rows -> ", db_out)
    }
    summary <- do.call(rbind, summaries)
    utils::write.csv(summary, file.path(output_dir, paste0(stem, "_matching_summary.csv")), row.names = FALSE)
    invisible(list(files = output_paths, summary = summary, sample_columns = sample_cols, taxonomy_columns_scanned = taxonomy_cols))
}

#' @keywords internal
plankannot_detect_sample_columns <- function(nms) {
    nms[grepl("^Sample[ _]?[0-9]+$", nms, ignore.case = FALSE)]
}

#' @keywords internal
plankannot_find_db_dir <- function() {
    installed <- system.file("db", package = "plankannot")
    if (nzchar(installed) && dir.exists(installed)) return(installed)
    candidates <- c("inst/db", file.path(getwd(), "inst/db"))
    hit <- candidates[dir.exists(candidates)]
    if (length(hit) > 0) return(hit[[1]])
    stop("Could not locate plankannot functional databases.", call. = FALSE)
}

#' @keywords internal
plankannot_clean_for_export <- function(df) {
    df[] <- lapply(df, function(x) {
        if (!is.character(x)) return(x)
        x <- iconv(x, from = "", to = "UTF-8", sub = "")
        gsub("[\\x00-\\x08\\x0B\\x0C\\x0E-\\x1F]", "", x, perl = TRUE)
    })
    df
}
