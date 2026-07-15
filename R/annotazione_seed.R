# Simulazione tassonomica casuale per test_dataset.xlsx
# Compatibile con il repository plankannot (branch test)
#
# Struttura attesa:
# plankannot/
# ├── DESCRIPTION
# ├── annotated_results/
# │   ├── mixoplankton/pr2_version_5.0.0_SSU_mothur_mixoplankton_annotated.tsv
# │   ├── phytoplankton/pr2_version_5.0.0_SSU_mothur_phytoplankton_annotated.tsv
# │   ├── habs/pr2_version_5.0.0_SSU_mothur_habs_annotated.tsv
# │   └── copepoda/pr2_version_5.0.0_SSU_mothur_copepoda_annotated.tsv
# └── test_dataset.xlsx  (oppure indicare un percorso diverso sotto)

required_packages <- c("readxl", "readr", "dplyr", "purrr", "writexl")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
    install.packages(missing_packages)
}

library(readxl)
library(readr)
library(dplyr)
library(purrr)
library(writexl)

# -----------------------------------------------------------------------------
# PARAMETRI MODIFICABILI
# -----------------------------------------------------------------------------

# Se il file Excel è nella radice del repository, lascia questo valore.
input_file <- "test_dataset.xlsx"

# File prodotto dallo script.
output_file <- "test_dataset_annotated_PR2.xlsx"

# Seed per rendere l'estrazione casuale riproducibile.
random_seed <- 123

# FALSE = una sequenza PR2 non viene assegnata a più ASV.
# TRUE = la stessa sequenza può essere estratta più volte.
sample_with_replacement <- FALSE

# TRUE = usa soltanto righe con almeno un'annotazione funzionale (Matched == TRUE)
# nei database mixoplankton/phytoplankton/HABs/copepoda.
only_functionally_annotated <- TRUE

# -----------------------------------------------------------------------------
# FUNZIONI DI SUPPORTO
# -----------------------------------------------------------------------------

find_project_root <- function(start = getwd()) {
    current <- normalizePath(start, winslash = "/", mustWork = TRUE)

    repeat {
        description_file <- file.path(current, "DESCRIPTION")
        annotated_dir <- file.path(current, "annotated_results")

        if (file.exists(description_file) && dir.exists(annotated_dir)) {
            description_lines <- readLines(description_file, warn = FALSE)
            if (any(grepl("^Package:\\s*plankannot\\s*$", description_lines))) {
                return(current)
            }
        }

        parent <- dirname(current)
        if (identical(parent, current)) break
        current <- parent
    }

    stop(
        paste0(
            "Non trovo la radice del repository plankannot. ",
            "Apri plankannot.Rproj oppure imposta la working directory nella cartella del repository."
        )
    )
}

read_pr2_file <- function(path) {
    if (!file.exists(path)) {
        stop("File PR2 non trovato: ", path)
    }

    # readr normalmente gestisce il file; il fallback Latin1 evita problemi
    # con alcuni caratteri presenti nel database HABs.
    tryCatch(
        read_tsv(path, show_col_types = FALSE, progress = FALSE),
        error = function(e) {
            read_delim(
                path,
                delim = "\t",
                locale = locale(encoding = "Latin1"),
                show_col_types = FALSE,
                progress = FALSE
            )
        }
    )
}

as_logical_matched <- function(x) {
    if (is.logical(x)) return(replace(x, is.na(x), FALSE))
    tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

prepare_annotation <- function(df, prefix, rank_cols) {
    required <- c("seq_id", rank_cols)
    missing <- setdiff(required, names(df))
    if (length(missing) > 0) {
        stop(prefix, ": colonne mancanti: ", paste(missing, collapse = ", "))
    }

    if (!"Matched" %in% names(df)) {
        df$Matched <- FALSE
    }

    annotation_cols <- setdiff(names(df), c("seq_id", rank_cols, "Matched"))

    result <- df %>%
        mutate(.matched_tmp = as_logical_matched(Matched))

    if (only_functionally_annotated) {
        result <- result %>% filter(.matched_tmp)
    }

    result %>%
        select(seq_id, all_of(annotation_cols), .matched_tmp) %>%
        rename_with(~ paste0(prefix, "_", .x), all_of(annotation_cols)) %>%
        rename(!!paste0(prefix, "_Matched") := .matched_tmp) %>%
        distinct(seq_id, .keep_all = TRUE)
}

# -----------------------------------------------------------------------------
# 1. INDIVIDUAZIONE DEL REPOSITORY E DEI DATABASE
# -----------------------------------------------------------------------------

project_root <- find_project_root()

pr2_files <- c(
    mixoplankton = file.path(
        project_root, "annotated_results", "mixoplankton",
        "pr2_version_5.0.0_SSU_mothur_mixoplankton_annotated.tsv"
    ),
    phytoplankton = file.path(
        project_root, "annotated_results", "phytoplankton",
        "pr2_version_5.0.0_SSU_mothur_phytoplankton_annotated.tsv"
    ),
    habs = file.path(
        project_root, "annotated_results", "habs",
        "pr2_version_5.0.0_SSU_mothur_habs_annotated.tsv"
    ),
    copepoda = file.path(
        project_root, "annotated_results", "copepoda",
        "pr2_version_5.0.0_SSU_mothur_copepoda_annotated.tsv"
    )
)

missing_files <- pr2_files[!file.exists(pr2_files)]
if (length(missing_files) > 0) {
    stop(
        "Mancano uno o più database PR2:\n",
        paste("-", missing_files, collapse = "\n")
    )
}

# Rende relativi alla radice del progetto i percorsi semplici.
resolve_project_path <- function(path) {
    if (file.exists(path) || grepl("^([A-Za-z]:|/|~)", path)) return(path)
    file.path(project_root, path)
}

input_path <- resolve_project_path(input_file)
output_path <- if (grepl("^([A-Za-z]:|/|~)", output_file)) {
    output_file
} else {
    file.path(project_root, output_file)
}

if (!file.exists(input_path)) {
    stop("Dataset Excel non trovato: ", input_path)
}

# -----------------------------------------------------------------------------
# 2. LETTURA DEL DATASET DI CAMPIONAMENTO
# -----------------------------------------------------------------------------

samples <- read_excel(input_path)

sample_cols <- grep("^Sample[ _]?[1-6]$", names(samples), value = TRUE)
if (length(sample_cols) != 6) {
    stop(
        "Non trovo esattamente le sei colonne Sample 1-Sample 6. Colonne trovate: ",
        paste(names(samples), collapse = ", ")
    )
}

# Individua la colonna identificativa; se non esiste, la crea.
id_candidates <- c("ASV", "Component", "Component_ID", "component", "ID", "id")
id_col <- id_candidates[id_candidates %in% names(samples)][1]

if (is.na(id_col)) {
    samples <- samples %>% mutate(ASV = paste0("ASV_", row_number()), .before = 1)
    id_col <- "ASV"
}

# -----------------------------------------------------------------------------
# 3. LETTURA E UNIONE DEI QUATTRO DATABASE
# -----------------------------------------------------------------------------

pr2 <- map(pr2_files, read_pr2_file)
rank_cols <- paste0("rank", 1:9)

# La tassonomia PR2 è comune ai quattro file: la prendiamo dal primo.
taxonomy <- pr2$mixoplankton %>%
    select(seq_id, all_of(rank_cols)) %>%
    distinct(seq_id, .keep_all = TRUE)

annotations <- imap(
    pr2,
    ~ prepare_annotation(.x, .y, rank_cols)
)

pr2_complete <- reduce(
    annotations,
    full_join,
    by = "seq_id"
) %>%
    left_join(taxonomy, by = "seq_id") %>%
    relocate(seq_id, all_of(rank_cols))

matched_cols <- grep("_Matched$", names(pr2_complete), value = TRUE)

if (only_functionally_annotated) {
    sampling_pool <- pr2_complete %>%
        filter(if_any(all_of(matched_cols), ~ .x %in% TRUE))
} else {
    sampling_pool <- pr2_complete
}

sampling_pool <- sampling_pool %>%
    filter(!is.na(seq_id)) %>%
    distinct(seq_id, .keep_all = TRUE)

if (nrow(sampling_pool) == 0) {
    stop("Il pool PR2 risultante è vuoto.")
}

if (!sample_with_replacement && nrow(sampling_pool) < nrow(samples)) {
    stop(
        "Il pool contiene ", nrow(sampling_pool),
        " righe, ma il dataset ne contiene ", nrow(samples),
        ". Imposta sample_with_replacement <- TRUE."
    )
}

# -----------------------------------------------------------------------------
# 4. ASSEGNAZIONE CASUALE RIPRODUCIBILE
# -----------------------------------------------------------------------------

set.seed(random_seed)
selected_indices <- sample(
    seq_len(nrow(sampling_pool)),
    size = nrow(samples),
    replace = sample_with_replacement
)

random_annotations <- sampling_pool[selected_indices, , drop = FALSE] %>%
    mutate(
        Simulation_seed = random_seed,
        PR2_assignment = "random"
    )

result <- bind_cols(samples, random_annotations)

# -----------------------------------------------------------------------------
# 5. PULIZIA TESTI ED ESPORTAZIONE
# -----------------------------------------------------------------------------

# Pulisce caratteri non validi o problematici per XML/Excel.
clean_excel_text <- function(x) {
    if (!is.character(x)) {
        return(x)
    }

    clean_one <- function(value) {
        if (is.na(value)) {
            return(NA_character_)
        }

        # Converte nella codifica UTF-8.
        value <- iconv(value, from = "", to = "UTF-8", sub = "")

        # Sostituisce gli spazi Unicode speciali con uno spazio normale.
        value <- gsub("\u00A0", " ", value, fixed = TRUE)
        value <- gsub("[\u2007\u202F]", " ", value, perl = TRUE)

        # Conserva soltanto i caratteri ammessi da XML 1.0:
        # tab (9), newline (10), carriage return (13), e i normali intervalli Unicode.
        code_points <- utf8ToInt(value)
        valid <- code_points %in% c(9L, 10L, 13L) |
            (code_points >= 32L & code_points <= 55295L) |
            (code_points >= 57344L & code_points <= 65533L) |
            (code_points >= 65536L & code_points <= 1114111L)

        value <- intToUtf8(code_points[valid])
        trimws(value)
    }

    vapply(x, clean_one, character(1), USE.NAMES = FALSE)
}

result <- result %>%
    mutate(across(where(is.character), clean_excel_text))

simulation_info <- tibble(
    parameter = c(
        "input_file", "output_file", "number_of_components",
        "PR2_pool_size", "random_seed", "replacement",
        "only_functionally_annotated"
    ),
    value = as.character(c(
        input_path, output_path, nrow(samples), nrow(sampling_pool),
        random_seed, sample_with_replacement, only_functionally_annotated
    ))
) %>%
    mutate(across(where(is.character), clean_excel_text))

write_xlsx(
    list(
        annotated_dataset = result,
        simulation_info = simulation_info
    ),
    output_path
)

message("Operazione completata.")
message("Repository: ", project_root)
message("Componenti annotate: ", nrow(result))
message("Pool PR2 disponibile: ", nrow(sampling_pool))
message("File creato: ", output_path)
