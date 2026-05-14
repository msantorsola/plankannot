library(readr)
library(dplyr)
library(stringr)
library(tools)

# ============================================================
# 1. FUNZIONI DI PULIZIA GENERALI
# ============================================================

remove_prefix <- function(x) {
    str_remove(x, "^[a-z]__")
}

# Pulizia standard per il matching nei db
clean_taxon <- function(x) {
    x %>%
        as.character() %>%
        remove_prefix() %>%
        str_replace_all(",", "") %>%
        str_replace_all(" ", "") %>%
        tolower()
}

rename_ordo_to_order <- function(df) {
    if ("Ordo" %in% names(df)) {
        names(df)[names(df) == "Ordo"] <- "Order"
    }
    df
}

# Pulizia nomi MIDORI/taxon:
#   - rimuove prefisso "order_"
#   - rimuove taxid NCBI numerico alla fine (es. _2759)
clean_midori_name <- function(x) {
    x %>%
        str_remove("^order_") %>%      # order_Vannellidae_95227 -> Vannellidae_95227
        str_remove("_[0-9]+$") %>%     # Vannellidae_95227       -> Vannellidae
        str_squish()
}


# ============================================================
# 2. ANNOTAZIONE VETTORIZZATA (veloce)
# ============================================================

annotate_user_dataset <- function(dataset_path, db_name,
                                  db_dir = "inst/db/",
                                  base_output_dir = "annotated_results") {

    input_name <- file_path_sans_ext(basename(dataset_path))
    ds <- read_tsv(dataset_path, show_col_types = FALSE)

    db_path <- file.path(db_dir, paste0(db_name, ".rds"))
    if (!file.exists(db_path)) stop("Database not found: ", db_path)
    db <- readRDS(db_path)

    ds <- rename_ordo_to_order(ds)
    db <- rename_ordo_to_order(db)

    match_col <- if (db_name == "major") "Group" else "Species_Name"
    if (!match_col %in% names(db)) {
        stop("Missing expected match column in DB: ", match_col,
             "\nColonne disponibili: ", paste(names(db), collapse = ", "))
    }

    db$.match_value <- clean_taxon(db[[match_col]])
    db <- db[!duplicated(db$.match_value), ]

    # Colonne su cui cercare match (dal piu' specifico al meno specifico)
    search_cols <- c("Species", "Genus", "Order", "Family", "Class",
                     "Division", "Supergroup", "Kingdom")
    search_cols <- intersect(search_cols, names(ds))

    if (length(search_cols) == 0) {
        search_cols <- names(ds)[sapply(ds, is.character)]
    }

    message(sprintf("  [%s] Matching su: %s", db_name,
                    paste(search_cols, collapse = ", ")))

    matched_idx  <- rep(NA_integer_,   nrow(ds))
    matched_rank <- rep(NA_character_, nrow(ds))

    # Loop sulle COLONNE (vettorizzato) - non sulle righe
    for (col in search_cols) {
        unmatched <- which(is.na(matched_idx))
        if (length(unmatched) == 0) break

        vals_clean <- clean_taxon(ds[[col]][unmatched])
        idx        <- match(vals_clean, db$.match_value)
        found      <- !is.na(idx)

        matched_idx[unmatched[found]]  <- idx[found]
        matched_rank[unmatched[found]] <- col

        if (sum(found) > 0) {
            message(sprintf("  [%s] '%s': %d match (totale: %d/%d)",
                            db_name, col, sum(found),
                            sum(!is.na(matched_idx)), nrow(ds)))
        }
    }

    db_extra_cols        <- setdiff(names(db), c(match_col, ".match_value"))
    db_matched           <- db[matched_idx, db_extra_cols, drop = FALSE]
    rownames(db_matched) <- NULL

    final <- bind_cols(ds, db_matched, tibble(Matched_Rank = matched_rank))

    output_dir <- file.path(base_output_dir, db_name)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    output_file <- file.path(output_dir,
                             paste0(input_name, "_", db_name, "_annotated.tsv"))
    write_tsv(final, output_file)

    n_matched <- sum(!is.na(matched_rank))
    message(sprintf("  [%s] %d matchati (%.1f%%) | %d non matchati",
                    db_name, n_matched,
                    100 * n_matched / nrow(ds), nrow(ds) - n_matched))
    message("  Salvato: ", output_file)
}


# ============================================================
# 3. PARSER PR2 .tax
#
# Formato: seq_id \t Kingdom;Supergroup;Division;Class;Order;Family;Genus;Species;
# ============================================================

parse_pr2_tax <- function(tax_path) {

    raw <- read_tsv(tax_path,
                    col_names      = c("seq_id", "taxonomy"),
                    col_types      = cols(.default = col_character()),
                    show_col_types = FALSE)

    tax_levels          <- c("Kingdom", "Supergroup", "Division", "Class",
                             "Order", "Family", "Genus", "Species")
    tax_split           <- str_split_fixed(raw$taxonomy, ";", n = 9)[, 1:8]
    colnames(tax_split) <- tax_levels
    tax_split[tax_split == ""] <- NA

    result <- bind_cols(raw["seq_id"],
                        as.data.frame(tax_split, stringsAsFactors = FALSE))

    message("  Sequenze: ", nrow(result))
    return(result)
}


# ============================================================
# 4. PARSER MIDORI .tax e .taxon
#
# Formato: seq_id \t Kingdom_ID;Class_ID;Order_ID;Family_ID;Genus_ID;Species_ID;
# Pulizia: rimuove taxid numerici (Nome_12345 -> Nome) e prefisso "order_"
# ============================================================

parse_midori_tax <- function(tax_path) {

    raw <- read_tsv(tax_path,
                    col_names      = c("seq_id", "taxonomy"),
                    col_types      = cols(.default = col_character()),
                    show_col_types = FALSE)

    tax_levels          <- c("Kingdom", "Class", "Order", "Family", "Genus", "Species")
    tax_split           <- str_split_fixed(raw$taxonomy, ";", n = 7)[, 1:6]
    colnames(tax_split) <- tax_levels

    # Rimuove taxid da ogni cella
    tax_split           <- apply(tax_split, 2, clean_midori_name)
    tax_split[tax_split == "" | tax_split == "NA"] <- NA

    result <- bind_cols(raw["seq_id"],
                        as.data.frame(tax_split, stringsAsFactors = FALSE))

    message("  Sequenze: ", nrow(result))
    return(result)
}


# ============================================================
# 5. AUTO-DETECT FORMATO E PARSING
#
# Regole di rilevamento automatico dal nome file:
#   contiene "midori"  -> midori
#   estensione .taxon  -> midori (stesso formato)
#   tutto il resto     -> pr2
# ============================================================

detect_format <- function(tax_path) {
    fname <- tolower(basename(tax_path))
    if (str_detect(fname, "midori")) return("midori")
    if (str_ends(fname, "\\.taxon"))  return("midori")
    return("pr2")
}

parse_tax_file <- function(tax_path, format = NULL) {

    if (is.null(format)) format <- detect_format(tax_path)

    message("Formato: ", format, " | File: ", basename(tax_path))

    switch(format,
           pr2    = parse_pr2_tax(tax_path),
           midori = parse_midori_tax(tax_path),
           stop("Formato non supportato: '", format, "'. Usa 'pr2' o 'midori'.")
    )
}


# ============================================================
# 6. FUNZIONE PRINCIPALE
#
# Processa un singolo file .tax / .taxon contro tutti i db
# ============================================================

annotate_tax_with_dbs <- function(tax_path,
                                  db_names,
                                  format     = NULL,
                                  db_dir     = "inst/db/",
                                  output_dir = "annotated_results/",
                                  tmp_dir    = tempdir()) {

    stopifnot(file.exists(tax_path))
    stopifnot(length(db_names) > 0)

    # Step 1: parsing
    message("\n=== Step 1: Parsing ===")
    parsed  <- parse_tax_file(tax_path, format)

    # Step 2: salva TSV temporaneo
    tmp_tsv <- file.path(tmp_dir,
                         paste0(file_path_sans_ext(basename(tax_path)), ".tsv"))
    write_tsv(parsed, tmp_tsv)
    message("TSV temporaneo: ", tmp_tsv)

    # Step 3: annotazione con ogni db
    message("\n=== Step 2: Annotazione ===")

    results <- lapply(db_names, function(db_name) {
        status <- tryCatch({
            annotate_user_dataset(
                dataset_path    = tmp_tsv,
                db_name         = db_name,
                db_dir          = db_dir,
                base_output_dir = output_dir
            )
            "ok"
        }, error = function(e) paste("error:", e$message))

        tibble::tibble(dataset  = basename(tax_path),
                       database = db_name,
                       status   = status)
    })

    dplyr::bind_rows(results)
}


# ============================================================
# ESECUZIONE
# Trova automaticamente tutti i file .tax e .taxon nella root
# ============================================================

db_names <- c("copepoda", "habs", "mixoplankton", "phytoplankton")

# Cerca tutti i file .tax e .taxon nella cartella corrente
tax_files <- list.files(".", pattern = "\\.taxon$", full.names = TRUE)

if (length(tax_files) == 0) stop("Nessun file .tax o .taxon trovato nella cartella corrente")

message("File trovati: ", length(tax_files))
for (f in tax_files) message("  ", basename(f))

# Processa ogni file
all_results <- lapply(tax_files, function(f) {
    message("\n========================================")
    message("Elaborazione: ", basename(f))
    message("========================================")
    annotate_tax_with_dbs(
        tax_path   = f,
        db_names   = db_names,
        db_dir     = "inst/db/",
        output_dir = "annotated_results/"
    )
})

# Tabella riassuntiva finale
message("\n=== RIEPILOGO ===")
print(dplyr::bind_rows(all_results))
