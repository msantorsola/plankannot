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

# Rimuove il suffisso "_EXT" tipico del formato MZGmothur
# Es: "Crustacea_EXT" -> "Crustacea"
#     "Aetideopsis_EXT" -> "Aetideopsis"
#     "Aetideopsis_multiserrata" -> "Aetideopsis_multiserrata"  (invariato)
strip_ext <- function(x) {
    str_remove(x, "_EXT$")
}


# ============================================================
# 2. PARSE FILE MZGmothur .txt
#
# Formato di ogni riga:
#   seq_id<TAB>rank1;rank2;...;rankN;
#
# Differenze rispetto a PR2 e MIDORI2:
#   - seq_id usa "__" come separatore interno (AccNum__Species)
#   - alcuni rank hanno suffisso "_EXT" (nodi interni espansi)
#   - nessun TaxID numerico
# ============================================================

parse_mzgmothur_txt <- function(txt_path) {

    stopifnot(file.exists(txt_path))

    raw <- read_tsv(
        txt_path,
        col_names      = c("seq_id", "taxonomy"),
        col_types      = cols(.default = col_character()),
        show_col_types = FALSE
    )

    # Splitta su ";" e rimuove token vuoti/NA
    tax_split <- str_split(raw$taxonomy, ";")
    tax_split <- lapply(tax_split, function(v) v[v != "" & !is.na(v)])

    # Rimuove il suffisso "_EXT" da ogni token
    tax_split_clean <- lapply(tax_split, strip_ext)

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
    # Nei nomi specie gli underscore sostituiscono gli spazi:
    # "Aetideopsis_multiserrata" -> corrisponde a "Aetideopsis multiserrata" nel db
    # La normalizzazione (spazio vs underscore) avviene in clean_taxon al momento del match
    last_value <- apply(tax_df, 1, function(r) {
        vals <- r[!is.na(r)]
        if (length(vals) == 0) NA_character_ else tail(vals, 1)
    })

    result <- bind_cols(
        raw["seq_id"],
        as.data.frame(tax_df, stringsAsFactors = FALSE),
        tibble(match_query = last_value)
    )

    message("MZGmothur parsed: ", nrow(result), " sequenze, ",
            max_ranks, " livelli (rank1-rank", max_ranks, ")")
    return(result)
}


# ============================================================
# 3. ANNOTAZIONE: match sul valore finale del .txt
#    Gestisce sia "Genus_species" (underscore) sia "Genus species" (spazio)
# ============================================================

annotate_with_db <- function(mzg_df, db_name,
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

    # Normalizza la colonna di match nel db:
    # lowercase + sostituisce underscore con spazio per uniformità
    db$.match_value <- clean_taxon(db[[match_col]]) %>%
        str_replace_all("_", " ")
    db <- db[!duplicated(db$.match_value), ]

    # Colonne extra del db (esclude match_col e .match_value)
    db_extra_cols <- setdiff(names(db), c(match_col, ".match_value"))

    # Normalizza il valore di query:
    # lowercase + sostituisce underscore con spazio
    query_clean <- clean_taxon(mzg_df$match_query) %>%
        str_replace_all("_", " ")

    # Match vettorizzato
    idx <- match(query_clean, db$.match_value)

    n_matched <- sum(!is.na(idx))
    message(sprintf("  [%s] %d/%d sequenze matchate (%.1f%%)",
                    db_name, n_matched, nrow(mzg_df),
                    100 * n_matched / nrow(mzg_df)))

    db_matched <- db[idx, db_extra_cols, drop = FALSE]
    rownames(db_matched) <- NULL

    rank_cols <- grep("^rank", names(mzg_df), value = TRUE)

    final <- bind_cols(
        mzg_df[, c("seq_id", rank_cols)],
        db_matched,
        tibble(Matched = !is.na(idx))
    )

    return(final)
}


# ============================================================
# 4. FUNZIONE PRINCIPALE
# ============================================================

annotate_mzgmothur_with_dbs <- function(txt_path,
                                        db_names,
                                        db_dir     = "inst/db/",
                                        output_dir = "annotated_results/") {

    stopifnot(file.exists(txt_path))
    stopifnot(length(db_names) > 0)

    input_name <- file_path_sans_ext(basename(txt_path))

    # Step 1: parsing
    message("\n=== Step 1: Parsing MZGmothur .txt ===")
    mzg <- parse_mzgmothur_txt(txt_path)

    # Step 2: annotazione per ogni db
    message("\n=== Step 2: Annotazione ===")

    results <- lapply(db_names, function(db_name) {

        status <- tryCatch({

            annotated <- annotate_with_db(
                mzg_df  = mzg,
                db_name = db_name,
                db_dir  = db_dir
            )

            out_dir  <- file.path(output_dir, db_name)
            dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
            out_file <- file.path(out_dir,
                                  paste0(input_name, "_", db_name, "_annotated.tsv"))
            write_tsv(annotated, out_file)
            message("  Salvato: ", out_file)
            "ok"

        }, error = function(e) {
            message("  ERRORE [", db_name, "]: ", e$message)
            paste("error:", e$message)
        })

        tibble::tibble(
            dataset  = basename(txt_path),
            database = db_name,
            status   = status
        )
    })

    dplyr::bind_rows(results)
}


# ============================================================
# ESECUZIONE
# ============================================================

# Processa tutti i file *.txt nella directory corrente
txt_files <- list.files(pattern = "\\.txt$", full.names = TRUE)

if (length(txt_files) == 0) {
    stop("Nessun file .txt trovato nella directory corrente.")
}

message("File .txt trovati: ", length(txt_files))
for (f in txt_files) message("  - ", f)

results <- lapply(txt_files, function(txt_path) {
    message("\n>>> Elaborazione: ", basename(txt_path))
    annotate_mzgmothur_with_dbs(
        txt_path   = txt_path,
        db_names   = c("copepoda", "habs", "mixoplankton", "phytoplankton"),
        db_dir     = "inst/db/",
        output_dir = "annotated_results/"
    )
})

result <- dplyr::bind_rows(results)
print(result)
