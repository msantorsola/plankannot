library(readr)
library(dplyr)
library(stringr)
library(tools)

# ============================================================
# 1. FUNZIONI DI PULIZIA (da annotate_dataset.R)
# ============================================================

# Rimuove prefissi tassonomici tipo p__, s__, c__, ...
remove_prefix <- function(x) {
    str_remove(x, "^[a-z]__")
}

# Normalizza nomi tassonomici: rimuove prefissi, virgole, spazi, lowercase
clean_taxon <- function(x) {
    x %>%
        as.character() %>%
        remove_prefix() %>%
        str_replace_all(",", "") %>%
        str_replace_all(" ", "") %>%
        tolower()
}

# Rinomina colonna "Ordo" in "Order" se presente
rename_ordo_to_order <- function(df) {
    if ("Ordo" %in% names(df)) {
        names(df)[names(df) == "Ordo"] <- "Order"
    }
    df
}


# ============================================================
# 2. ANNOTAZIONE SINGOLO DATASET + SINGOLO DB (da annotate_dataset.R)
# ============================================================

annotate_user_dataset <- function(dataset_path, db_name,
                                  db_dir = "data/",
                                  base_output_dir = "annotated_results") {

    input_name <- file_path_sans_ext(basename(dataset_path))
    ds <- read_tsv(dataset_path, show_col_types = FALSE)

    # Carica il database RDS dal pacchetto installato
    db_path <- file.path("inst/db", paste0(db_name, ".rds"))
    if (!file.exists(db_path)) stop("Database not found: ", db_path)
    db <- readRDS(db_path)

    ds <- rename_ordo_to_order(ds)
    db <- rename_ordo_to_order(db)

    # Colonna di match: "Group" per il db "major", altrimenti "Species_Name"
    match_col <- if (db_name == "major") "Group" else "Species_Name"
    if (!match_col %in% names(db)) {
        stop(paste("Missing expected match column in DB:", match_col))
    }

    # Aggiunge colonna normalizzata al db per il confronto
    db$.match_value <- clean_taxon(db[[match_col]])

    output_rows <- vector("list", nrow(ds))
    ds_cols <- names(ds)

    # Loop riga per riga: cerca match su ogni colonna del dataset
    for (i in seq_len(nrow(ds))) {
        ds_row  <- ds[i, , drop = FALSE]
        matched <- FALSE

        for (col in ds_cols) {
            val <- ds_row[[col]]
            if (is.na(val) || val == "") next

            val_clean   <- clean_taxon(val)
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

    final      <- bind_rows(output_rows)
    output_dir <- file.path(base_output_dir, db_name)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    output_file <- file.path(output_dir, paste0(input_name, "_", db_name, "_annotated.tsv"))
    write_tsv(final, output_file)
    message("Saved: ", output_file)
}


# ============================================================
# 3. PARSE FILE PR2 .tax
#
# Input:  file tab-separato con 2 colonne:
#           seq_id   (es. AB353770.1.1740_U)
#           taxonomy (es. Eukaryota;TSAR;Alveolata;...;Species;)
# Output: data frame con colonne esplicite per i 8 livelli tassonomici
# ============================================================

parse_pr2_tax <- function(tax_path) {

    stopifnot(file.exists(tax_path))

    raw <- read_tsv(
        tax_path,
        col_names      = c("seq_id", "taxonomy"),
        col_types      = cols(.default = col_character()),
        show_col_types = FALSE
    )

    tax_levels <- c("Kingdom", "Supergroup", "Division", "Class",
                    "Order", "Family", "Genus", "Species")

    # Separa la stringa tassonomica nei livelli (l'ultimo campo dopo ; e' vuoto)
    tax_split           <- str_split_fixed(raw$taxonomy, ";", n = 9)[, 1:8]
    colnames(tax_split) <- tax_levels
    tax_split[tax_split == ""] <- NA

    result <- bind_cols(
        raw["seq_id"],
        as.data.frame(tax_split, stringsAsFactors = FALSE)
    )

    message("PR2 parsed: ", nrow(result), " sequenze")
    return(result)
}


# ============================================================
# 4. FUNZIONE PRINCIPALE
#
# Flusso:
#   1. Parsa il file .tax PR2 -> data frame con colonne tassonomiche
#   2. Salva TSV temporaneo (input atteso da annotate_user_dataset)
#   3. Per ogni db chiama annotate_user_dataset()
#   4. Restituisce tabella di status (ok / error) per ogni db
# ============================================================

annotate_pr2_with_dbs <- function(tax_path,
                                  db_names,
                                  output_dir = "annotated_results/",
                                  tmp_dir    = tempdir()) {

    stopifnot(file.exists(tax_path))
    stopifnot(length(db_names) > 0)

    # Step 1: parsing PR2
    message("\n=== Step 1: Parsing PR2 ===")
    pr2 <- parse_pr2_tax(tax_path)

    # Step 2: salva TSV temporaneo
    tmp_tsv <- file.path(
        tmp_dir,
        paste0(file_path_sans_ext(basename(tax_path)), ".tsv")
    )
    write_tsv(pr2, tmp_tsv)
    message("TSV temporaneo: ", tmp_tsv)

    # Step 3: annotazione con ogni db
    message("\n=== Step 2: Annotazione con i database ===")

    results <- lapply(db_names, function(db_name) {
        status <- tryCatch({
            annotate_user_dataset(
                dataset_path    = tmp_tsv,
                db_name         = db_name,
                base_output_dir = output_dir
            )
            "ok"
        }, error = function(e) {
            paste("error:", e$message)
        })

        message(sprintf("  [%s] -> %s", db_name, status))

        tibble::tibble(
            dataset  = basename(tax_path),
            database = db_name,
            status   = status
        )
    })

    dplyr::bind_rows(results)
}


# ============================================================
# ESECUZIONE
# ============================================================

result <- annotate_pr2_with_dbs(
    tax_path   = "pr2_version_5.0.0_SSU_mothur.tax",
    db_names   = c("copepoda", "habs", "mixoplankton", "phytoplankton"),
    output_dir = "annotated_results/"
)

print(result)
