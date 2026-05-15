library(readr)
library(dplyr)
library(stringr)
library(tools)

# ============================================================
# 1. FUNZIONI DI PULIZIA
# ============================================================

clean_taxon <- function(x) {
    # Normalizza solo case e spazi bianchi.
    # NON tocca gli underscore: nel .tax le specie sono "Isospora_sp."
    # e nel db devono essere identiche ("Isospora_sp.", "Isospora_spp.", ecc.)
    x %>%
        as.character() %>%
        str_squish() %>%
        tolower()
}

# ============================================================
# 2. PARSE FILE PR2 .tax
#
# Ogni riga:  seq_id<TAB>rank1;rank2;...;rankN;
# I livelli vengono chiamati rank1, rank2, ... rankN
# L'ultimo valore non vuoto è il candidato per il match.
# ============================================================

parse_pr2_tax <- function(tax_path) {

    stopifnot(file.exists(tax_path))

    raw <- read_tsv(
        tax_path,
        col_names      = c("seq_id", "taxonomy"),
        col_types      = cols(.default = col_character()),
        show_col_types = FALSE
    )

    # Splitta su ";" e rimuove l'eventuale token vuoto finale
    tax_split <- str_split(raw$taxonomy, ";")
    tax_split <- lapply(tax_split, function(v) v[v != "" & !is.na(v)])

    # Numero massimo di livelli presenti nel file
    max_ranks <- max(lengths(tax_split))
    rank_names <- paste0("rank", seq_len(max_ranks))

    # Costruisce data.frame con colonne rank1..rankN (NA se mancante)
    tax_df <- do.call(rbind, lapply(tax_split, function(v) {
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

    message("PR2 parsed: ", nrow(result), " sequenze, ", max_ranks, " livelli (rank1-rank", max_ranks, ")")
    return(result)
}


# ============================================================
# 3. ANNOTAZIONE: match solo sul valore finale del .tax
# ============================================================

annotate_with_db <- function(pr2_df, db_name,
                             db_dir = "inst/db/") {

    db_path <- file.path(db_dir, paste0(db_name, ".rds"))
    if (!file.exists(db_path)) stop("Database non trovato: ", db_path)
    db <- readRDS(db_path)

    # Colonna di match nel db (Species_Name per tutti tranne "major")
    match_col <- if (db_name == "major") "Group" else "Species_Name"
    if (!match_col %in% names(db)) {
        stop("Colonna '", match_col, "' non trovata nel db '", db_name,
             "'. Colonne disponibili: ", paste(names(db), collapse = ", "))
    }

    # Normalizza la colonna di match nel db
    db$.match_value <- clean_taxon(db[[match_col]])
    # Rimuove duplicati (tiene prima occorrenza)
    db <- db[!duplicated(db$.match_value), ]

    # Colonne del db da aggiungere all'output
    # Esclude match_col e .match_value (metadati interni)
    db_extra_cols <- setdiff(names(db), c(match_col, ".match_value"))

    # Normalizza il valore di query dal .tax
    query_clean <- clean_taxon(pr2_df$match_query)

    # Match vettorizzato
    idx <- match(query_clean, db$.match_value)

    n_matched <- sum(!is.na(idx))
    message(sprintf("  [%s] %d/%d sequenze matchate (%.1f%%)",
                    db_name, n_matched, nrow(pr2_df),
                    100 * n_matched / nrow(pr2_df)))

    # Costruisce output:
    # - colonne rank1..rankN del .tax (senza match_query che è interna)
    # - colonne funzionali dal db (esattamente come sono nel .rds)
    # - colonna Matched (TRUE/FALSE)
    db_matched <- db[idx, db_extra_cols, drop = FALSE]
    rownames(db_matched) <- NULL

    rank_cols <- grep("^rank", names(pr2_df), value = TRUE)

    final <- bind_cols(
        pr2_df[, c("seq_id", rank_cols)],
        db_matched,
        tibble(Matched = !is.na(idx))
    )

    return(final)
}


# ============================================================
# 4. FUNZIONE PRINCIPALE
# ============================================================

annotate_pr2_with_dbs <- function(tax_path,
                                  db_names,
                                  db_dir     = "inst/db/",
                                  output_dir = "annotated_results/") {

    stopifnot(file.exists(tax_path))
    stopifnot(length(db_names) > 0)

    input_name <- file_path_sans_ext(basename(tax_path))

    # Step 1: parsing
    message("\n=== Step 1: Parsing PR2 .tax ===")
    pr2 <- parse_pr2_tax(tax_path)

    # Step 2: annotazione per ogni db
    message("\n=== Step 2: Annotazione ===")

    results <- lapply(db_names, function(db_name) {

        status <- tryCatch({

            annotated <- annotate_with_db(
                pr2_df  = pr2,
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
