#' Annotate a single dataset using one or more functional databases
#'
#' This function applies annotation to a taxonomic dataset using a list of specified databases.
#'
#' @param dataset_path Path to the input dataset (.tsv)
#' @param db_names Character vector of database names (e.g., c("copepoda", "mixo"))
#' @param db_dir Path to the directory containing .rds DB files (default = "data/")
#' @param output_dir Path to the output directory (default = "annotated_results/")
#'
#' @return A data frame with database name and annotation status ("ok" or error message)
#' @export
annotate_dataset_with_dbs <- function(dataset_path, db_names,
                                      db_dir = "data/",
                                      output_dir = "annotated_results/") {
    stopifnot(file.exists(dataset_path))
    stopifnot(length(db_names) > 0)

    results <- lapply(db_names, function(db_name) {
        status <- tryCatch({
            annotate_user_dataset(
                dataset_path = dataset_path,
                db_name = db_name,
                db_dir = db_dir,
                base_output_dir = output_dir
            )
            "ok"
        }, error = function(e) {
            paste("error:", e$message)
        })

        tibble::tibble(
            dataset = basename(dataset_path),
            database = db_name,
            status = status
        )
    })

    dplyr::bind_rows(results)
}
