library(readr)
library(dplyr)
library(stringr)
library(tools)

# ============================================================
# 1. FUNZIONI DI PULIZIA
# ============================================================

remove_prefix <- function(x) {
    str_remove(x, "^[a-z]__")
}

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


# ============================================================
# 2. ANNOTAZIONE VETTORIZZATA
#
# Invece di un loop riga per riga, usa match() su tutto
# il vettore in una sola operazione
# ============================================================

annotate_user_dataset <- function(dataset_path, db_name,
                                  db_dir = "inst/db/",
                                  base_output_dir = "annotated_results") {

    input_name <- file_path_sans_ext(basename(dataset_path))
    ds <- read_tsv(dataset_path, show_col_types = FALSE)

    # Carica db dalla cartella inst/db/ locale
    db_path <- file.path(db_dir, paste0(db_name, ".rds"))
    if (!file.exists(db_path)) stop("Database not found: ", db_path)
    db <- readRDS(db_path)

    ds <- rename_ordo_to_order(ds)
    db <- rename_ordo_to_order(db)

    # Colonna di match nel db
    match_col <- if (db_name == "major") "Group" else "Species_Name"
    if (!match_col %in% names(db)) {
        stop("Missing expected match column in DB: ", match_col,
             "\nColonne disponibili: ", paste(names(db), collapse = ", "))
    }

    # Colonna normalizzata nel db per il confronto
    db$.match_value <- clean_taxon(db[[match_col]])
    # Rimuove duplicati nel db sulla colonna di match (tiene prima occorrenza)
    db <- db[!duplicated(db$.match_value), ]

    # Colonne del dataset su cui cercare (solo quelle che esistono)
    search_cols <- c("Species", "Genus", "Order", "Family", "Class",
                     "Division", "Supergroup", "Kingdom")
    search_cols <- intersect(search_cols, names(ds))

    if (length(search_cols) == 0) {
        # Fallback: usa tutte le colonne carattere del dataset
        search_cols <- names(ds)[sapply(ds, is.character)]
    }

    message(sprintf("  [%s] Matching su colonne: %s", db_name,
                    paste(search_cols, collapse = ", ")))

    # Inizializza vettori di risultato
    matched_idx  <- rep(NA_integer_,   nrow(ds))  # indice riga nel db
    matched_rank <- rep(NA_character_, nrow(ds))  # colonna su cui e' avvenuto il match

    # Loop sulle COLONNE (non sulle righe!) - molto più veloce
    # Per ogni colonna, fa match vettorizzato su tutte le righe ancora non matchate
    for (col in search_cols) {

        # Righe ancora senza match
        unmatched <- which(is.na(matched_idx))
        if (length(unmatched) == 0) break  # tutte matchate, esci

        # Vettore pulito dei valori in questa colonna per le righe non ancora matchate
        vals_clean <- clean_taxon(ds[[col]][unmatched])

        # match() vettorizzato: trova posizione nel db per ogni valore
        idx <- match(vals_clean, db$.match_value)

        # Aggiorna solo le righe dove e' stato trovato un match
        found <- !is.na(idx)
        matched_idx[unmatched[found]]  <- idx[found]
        matched_rank[unmatched[found]] <- col

        n_found <- sum(found)
        if (n_found > 0) {
            message(sprintf("  [%s] '%s': %d match trovati (totale matchati: %d/%d)",
                            db_name, col, n_found,
                            sum(!is.na(matched_idx)), nrow(ds)))
        }
    }

    # Costruisce output finale
    # Colonne del db da aggiungere (esclude .match_value e match_col)
    db_extra_cols <- setdiff(names(db), c(match_col, ".match_value"))

    # Prende le righe del db corrispondenti (NA per le righe senza match)
    db_matched <- db[matched_idx, db_extra_cols, drop = FALSE]
    rownames(db_matched) <- NULL

    final <- bind_cols(
        ds,
        db_matched,
        tibble(Matched_Rank = matched_rank)
    )

    # Salva output
    output_dir <- file.path(base_output_dir, db_name)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    output_file <- file.path(output_dir,
                             paste0(input_name, "_", db_name, "_annotated.tsv"))
    write_tsv(final, output_file)

    n_matched   <- sum(!is.na(matched_rank))
    n_unmatched <- nrow(ds) - n_matched
    message(sprintf("  [%s] Totale: %d matchati (%.1f%%), %d non matchati",
                    db_name, n_matched,
                    100 * n_matched / nrow(ds), n_unmatched))
    message("  Salvato: ", output_file)
}


# ============================================================
# 3. PARSE FILE PR2 .tax
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
# ============================================================

annotate_pr2_with_dbs <- function(tax_path,
                                  db_names,
                                  db_dir     = "inst/db/",
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
                db_dir          = db_dir,
                base_output_dir = output_dir
            )
            "ok"
        }, error = function(e) {
            paste("error:", e$message)
        })

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
    db_dir     = "inst/db/",
    output_dir = "annotated_results/"
)

print(result)
