library(readr)
library(dplyr)
library(stringr)
library(tools)


# standardisation of taxa name (e.g. remove p__, s__, ...)
remove_prefix <- function(x) {
    str_remove(x, "^[a-z]__")
}

# clean data, remove typos like commas and space
clean_taxon <- function(x) {
    x %>%
        as.character() %>%
        remove_prefix() %>%
        str_replace_all(",", "") %>%
        str_replace_all(" ", "") %>%
        tolower()
}

# standardisation of rank names
rename_ordo_to_order <- function(df) {
    if ("Ordo" %in% names(df)) {
        names(df)[names(df) == "Ordo"] <- "Order"
    }
    df
}

# functional annotation
annotate_user_dataset <- function(dataset_path, db_name, db_dir = "data/", base_output_dir = "annotated_results") {
    input_name <- file_path_sans_ext(basename(dataset_path))
    ds <- read_tsv(dataset_path, show_col_types = FALSE)

    # === NEW: load internal DB from inst/db/
    db_path <- system.file("db", paste0(db_name, ".rds"), package = "plankannot")
    if (db_path == "") stop("Database not found in package: ", db_name)
    db <- readRDS(db_path)

    ds <- rename_ordo_to_order(ds)
    db <- rename_ordo_to_order(db)

    match_col <- if (db_name == "major") "Group" else "Species_Name"
    if (!match_col %in% names(db)) {
        stop(paste("Missing expected match column in DB:", match_col))
    }

    db$.match_value <- clean_taxon(db[[match_col]])

    output_rows <- vector("list", nrow(ds))
    ds_cols <- names(ds)

    for (i in seq_len(nrow(ds))) {
        ds_row <- ds[i, , drop = FALSE]
        matched <- FALSE

        for (col in ds_cols) {
            val <- ds_row[[col]]
            if (is.na(val) || val == "") next
            val_clean <- clean_taxon(val)

            matched_row <- db %>% filter(.match_value == val_clean)

            if (nrow(matched_row) > 0) {
                output_rows[[i]] <- bind_cols(
                    ds_row,
                    matched_row[1, setdiff(names(matched_row), ".match_value"), drop = FALSE],
                    tibble(Matched_Rank = col)
                )
                matched <- TRUE
                break
            }
        }

        if (!matched) {
            output_rows[[i]] <- bind_cols(ds_row, tibble(Matched_Rank = NA))
        }
    }

    final <- bind_rows(output_rows)

    output_dir <- file.path(base_output_dir, db_name)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    output_file <- file.path(output_dir, paste0(input_name, "_", db_name, "_annotated.tsv"))
    write_tsv(final, output_file)

    message(paste("Saved annotated file:", output_file))
}

