library(readr)
library(dplyr)
library(stringr)
library(tools)

# ============================================================
# 1. FUNZIONI DI PULIZIA
# ============================================================

clean_taxon <- function(x) {
    x %>%
        as.character() %>%
        str_squish() %>%
        tolower()
}

# Rimuove il suffisso "_<TaxID>" tipico di MIDORI2
# Es: "Vannellidae_95227"   -> "Vannellidae"
#     "Clydonella sawyeri_2201168" -> "Clydonella sawyeri"
#     "order_Vannellidae_95227"    -> "order_Vannellidae"  (solo ultimo _NUM rimosso)
strip_taxid <- function(x) {
    str_remove(x, "_\\d+$")
}


# ============================================================
# 2. PARSE FILE MIDORI2 .taxon
#
# Formato di ogni riga:
#   seq_id<TAB>rank1_ID;rank2_ID;...;rankN_ID;
#
# Differenze rispetto a PR2:
#   - ogni token contiene "_<TaxID>" da rimuovere
#   - il seq_id può contenere punti e angolari (MG559732.1.<1.>690)
# ============================================================

parse_midori2_taxon <- function(taxon_path) {

    stopifnot(file.exists(taxon_path))

    raw <- read_tsv(
        taxon_path,
        col_names      = c("seq_id", "taxonomy"),
        col_types      = cols(.default = col_character()),
        show_col_types = FALSE
    )

    # Splitta su ";" e rimuove token vuoti/NA
    tax_split <- str_split(raw$taxonomy, ";")
    tax_split <- lapply(tax_split, function(v) v[v != "" & !is.na(v)])

    # Rimuove il suffisso "_<TaxID>" da ogni token
    tax_split_clean <- lapply(tax_split, strip_taxid)

    # Numero massimo di livelli presenti nel file
    max_ranks <- max(lengths(tax_split_clean))
    rank_names <- paste0("rank", seq_len(max_ranks))

    # Costruisce data.frame con colonne rank1..rankN (NA se mancante)
    tax_df <- do.call(rbind, lapply(tax_split_clean, function(v) {
        length(v) <- max_ranks   # padding con NA
        v
    }))
    colnames(tax_df) <- rank_names

    # Ultimo valore non-NA per ogni riga -> colonna da usare nel match
    last_value <- apply(tax_df, 1, function(r) {
        vals <- r[!is.na(r)]
        if (length(vals) == 0) NA_character_ else tail(vals, 1)
    })

    result <- bind_cols(
        raw["seq_id"],
        as.data.frame(tax_df, stringsAsFactors = FALSE),
        tibble(match_query = last_value)
    )

    message("MIDORI2 parsed: ", nrow(result), " sequenze, ",
            max_ranks, " livelli (rank1-rank", max_ranks, ")")
    return(result)
}


# ============================================================
# 3. ANNOTAZIONE: match solo sul valore finale del .taxon
#    (identica alla versione PR2, riutilizzabile as-is)
# ============================================================

annotate_with_db <- function(midori2_df, db_name,
                             db_dir = "inst/db/") {

    db_path <- file.path(db_dir, paste0(db_name, ".rds"))
    if (!file.exists(db_path)) stop("Database non trovato: ", db_path)
    db <- readRDS(db_path)

    # Colonna di match nel db
    match_col <- if (db_name == "major") "Group" else "Species_Name"
    if (!match_col %in% names(db)) {
        stop("Colonna '", match_col, "' non trovata nel db '", db_name,
             "'. Colonne disponibili: ", paste(names(db), collapse = ", "))
    }

    # Normalizza la colonna di match nel db
    db$.match_value <- clean_taxon(db[[match_col]])
    db <- db[!duplicated(db$.match_value), ]

    # Colonne extra del db (esclude match_col e .match_value)
    db_extra_cols <- setdiff(names(db), c(match_col, ".match_value"))

    # Normalizza il valore di query
    query_clean <- clean_taxon(midori2_df$match_query)

    # Match vettorizzato
    idx <- match(query_clean, db$.match_value)

    n_matched <- sum(!is.na(idx))
    message(sprintf("  [%s] %d/%d sequenze matchate (%.1f%%)",
                    db_name, n_matched, nrow(midori2_df),
                    100 * n_matched / nrow(midori2_df)))

    db_matched <- db[idx, db_extra_cols, drop = FALSE]
    rownames(db_matched) <- NULL

    rank_cols <- grep("^rank", names(midori2_df), value = TRUE)

    final <- bind_cols(
        midori2_df[, c("seq_id", rank_cols)],
        db_matched,
        tibble(Matched = !is.na(idx))
    )

    return(final)
}


# ============================================================
# 4. FUNZIONE PRINCIPALE
# ============================================================

annotate_midori2_with_dbs <- function(taxon_path,
                                      db_names,
                                      db_dir     = "inst/db/",
                                      output_dir = "annotated_results/") {

    if (!file.exists(taxon_path)) {
        stop("File .taxon non trovato: ", taxon_path, call. = FALSE)
    }

    if (length(db_names) == 0) {
        stop("`db_names` deve contenere almeno un database.", call. = FALSE)
    }

    input_name <- file_path_sans_ext(basename(taxon_path))

    # Step 1: parsing
    message("\n=== Dataset: ", basename(taxon_path), " ===")
    message("=== Step 1: Parsing MIDORI2 .taxon ===")
    midori2 <- parse_midori2_taxon(taxon_path)

    # Step 2: annotazione per ogni db
    message("=== Step 2: Annotazione ===")

    results <- lapply(db_names, function(db_name) {

        status <- tryCatch({

            annotated <- annotate_with_db(
                midori2_df = midori2,
                db_name    = db_name,
                db_dir     = db_dir
            )

            out_dir <- file.path(output_dir, db_name)
            dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

            out_file <- file.path(
                out_dir,
                paste0(input_name, "_", db_name, "_annotated.tsv")
            )

            write_tsv(annotated, out_file)
            message("  Salvato: ", out_file)

            "ok"

        }, error = function(e) {
            message("  ERRORE [", db_name, "]: ", e$message)
            paste("error:", e$message)
        })

        tibble::tibble(
            dataset  = basename(taxon_path),
            database = db_name,
            status   = status
        )
    })

    dplyr::bind_rows(results)
}


# ============================================================
# 5. ANNOTAZIONE DI TUTTI I FILE *.taxon IN UNA CARTELLA
# ============================================================

annotate_all_midori2_taxon <- function(input_dir = ".",
                                       db_names,
                                       db_dir     = "inst/db/",
                                       output_dir = "annotated_results/",
                                       recursive  = FALSE) {

    if (!dir.exists(input_dir)) {
        stop("Cartella di input non trovata: ", input_dir, call. = FALSE)
    }

    taxon_files <- list.files(
        path       = input_dir,
        pattern    = "\\.taxon$",
        full.names = TRUE,
        recursive  = recursive,
        ignore.case = TRUE
    )

    if (length(taxon_files) == 0) {
        stop(
            "Nessun file *.taxon trovato nella cartella: ",
            normalizePath(input_dir, winslash = "/", mustWork = FALSE),
            call. = FALSE
        )
    }

    message(
        "\nTrovati ", length(taxon_files),
        " file *.taxon in: ",
        normalizePath(input_dir, winslash = "/", mustWork = FALSE)
    )

    all_results <- lapply(taxon_files, function(taxon_path) {

        tryCatch(
            annotate_midori2_with_dbs(
                taxon_path = taxon_path,
                db_names   = db_names,
                db_dir     = db_dir,
                output_dir = output_dir
            ),
            error = function(e) {
                message(
                    "\nERRORE durante l'elaborazione di ",
                    basename(taxon_path),
                    ": ",
                    e$message
                )

                tibble::tibble(
                    dataset  = basename(taxon_path),
                    database = NA_character_,
                    status   = paste("error:", e$message)
                )
            }
        )
    })

    summary <- dplyr::bind_rows(all_results)

    summary_file <- file.path(output_dir, "midori2_annotation_summary.tsv")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    write_tsv(summary, summary_file)

    message("\nRiepilogo salvato: ", summary_file)

    summary
}


# ============================================================
# ESECUZIONE
# ============================================================

result <- annotate_all_midori2_taxon(
    input_dir  = ".",
    db_names   = c("copepoda", "habs", "mixoplankton", "phytoplankton"),
    db_dir     = "inst/db/",
    output_dir = "annotated_results/",
    recursive  = FALSE
)

print(result)
