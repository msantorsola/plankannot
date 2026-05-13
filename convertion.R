library(readr)

# Cartella dove sono i tuoi TSV
tsv_dir <- "inst/db"         # cambia se i TSV sono altrove
rds_dir <- "inst/db"         # dove salvare gli RDS (può essere la stessa)

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

tsv_files <- list.files(tsv_dir, pattern = "\\.tsv$", full.names = TRUE)

if (length(tsv_files) == 0) stop("Nessun file .tsv trovato in: ", tsv_dir)

for (f in tsv_files) {
    nome <- tools::file_path_sans_ext(basename(f))
    df   <- read_tsv(f, show_col_types = FALSE)
    out  <- file.path(rds_dir, paste0(nome, ".rds"))
    saveRDS(df, out, compress = "xz")
    message("✓ ", nome, ".tsv → ", nome, ".rds  (", nrow(df), " righe)")
}

message("\nFatto! ", length(tsv_files), " file convertiti.")
