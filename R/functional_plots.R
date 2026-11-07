#' Generate functional stacked-bar plots and per-trait heatmaps
#'
#' Reads one database-specific annotation file and summarizes only the selected
#' categorical functional traits using sample abundances. Each trait gets one
#' stacked bar plot and one heatmap.
#'
#' @param annotated_data Data frame, or path to an `.xlsx`/`.tsv` annotation file.
#' @param db_name One of `mixoplankton`, `copepoda`, `phytoplankton`, `habs`.
#' @param output_dir Directory for plots.
#' @param sample_cols Optional sample columns. If `NULL`, all Sample_N columns are detected.
#' @param trait_map Optional named character vector/list defining the traits to plot.
#' @return A list containing the detected traits and plot data.
#' @export
plot_functional_profiles <- function(
    annotated_data,
    db_name,
    output_dir = "functional_plots",
    sample_cols = NULL,
    trait_map = NULL
) {
    if (!db_name %in% c("mixoplankton", "copepoda", "phytoplankton", "habs")) {
        stop("Unknown database: ", db_name, call. = FALSE)
    }

    df <- plankannot_read_annotation_input(annotated_data)
    names(df) <- plankannot_safe_utf8(names(df))
    if (is.null(sample_cols)) sample_cols <- plankannot_detect_sample_columns(names(df))
    if (length(sample_cols) == 0) stop("No Sample_N columns found.", call. = FALSE)

    missing_samples <- setdiff(sample_cols, names(df))
    if (length(missing_samples) > 0) stop("Missing sample columns: ", paste(missing_samples, collapse = ", "), call. = FALSE)

    if (is.null(trait_map)) trait_map <- plankannot_default_trait_map()[[db_name]]
    if (is.null(trait_map)) stop("No trait map available for ", db_name, call. = FALSE)

    trait_cols <- plankannot_resolve_traits(names(df), trait_map)
    trait_cols <- trait_cols[!is.na(trait_cols)]
    if (length(trait_cols) == 0) {
        stop("None of the requested functional traits were found in the annotation file.", call. = FALSE)
    }

    out_db <- file.path(output_dir, db_name)
    bar_dir <- file.path(out_db, "stacked_bar")
    dir.create(bar_dir, recursive = TRUE, showWarnings = FALSE)

    bar_data <- list()
    heat_data <- list()

    for (trait in names(trait_cols)) {
        col <- unname(trait_cols[[trait]])
        long <- plankannot_trait_long_data(df, col, trait, sample_cols)
        if (nrow(long) == 0) next

        summary <- long |>
            dplyr::group_by(Sample, Category) |>
            dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") |>
            dplyr::group_by(Sample) |>
            dplyr::mutate(
                Relative_abundance = {
                    total_abundance <- sum(Abundance, na.rm = TRUE)
                    if (is.finite(total_abundance) && total_abundance > 0) {
                        Abundance / total_abundance
                    } else {
                        rep(0, dplyr::n())
                    }
                }
            ) |>
            dplyr::ungroup()

        bar_data[[trait]] <- summary
        heat_data[[trait]] <- summary |>
            dplyr::mutate(Trait = trait, Feature = paste(trait, Category, sep = ": ")) |>
            dplyr::select(Trait, Feature, Sample, Relative_abundance)

        p <- ggplot2::ggplot(summary, ggplot2::aes(x = Sample, y = Relative_abundance, fill = Category)) +
            ggplot2::geom_col(width = 0.85) +
            ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
            ggplot2::labs(
                title = paste0(db_name, " — ", trait),
                x = NULL, y = "Relative abundance", fill = trait
            ) +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

        safe_name <- gsub("[^A-Za-z0-9]+", "_", trait)
        ggplot2::ggsave(file.path(bar_dir, paste0(safe_name, ".png")), p, width = 9, height = 5, dpi = 300)

        # One heatmap per selected trait, as requested (not one combined heatmap per DB).
        # Build the full Category x Sample grid. A category not observed in a
        # sample is a true zero, not an NA; this prevents blank white cells.
        category_levels <- unique(as.character(summary$Category))
        heat_trait <- tidyr::expand_grid(
            Category = category_levels,
            Sample = sample_cols
        ) |>
            dplyr::left_join(
                summary |>
                    dplyr::select(Sample, Category, Relative_abundance),
                by = c("Sample", "Category")
            ) |>
            dplyr::mutate(
                Relative_abundance = dplyr::coalesce(Relative_abundance, 0),
                Category = factor(Category, levels = rev(category_levels)),
                Sample = factor(Sample, levels = sample_cols)
            )

        heat_plot <- ggplot2::ggplot(
            heat_trait,
            ggplot2::aes(x = Sample, y = Category, fill = Relative_abundance)
        ) +
            ggplot2::geom_tile(color = "grey88", linewidth = 0.18) +
            ggplot2::scale_fill_viridis_c(
                option = "D",
                direction = 1,
                limits = c(0, 1),
                breaks = c(0, 0.25, 0.5, 0.75, 1),
                labels = function(x) paste0(round(100 * x), "%"),
                na.value = "grey75",
                name = "Abbondanza\nrelativa"
            ) +
            ggplot2::labs(
                title = paste0(db_name, " — ", trait),
                x = NULL,
                y = trait,
                subtitle = "Le celle senza rilevamento sono rappresentate come zero"
            ) +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "grey20"),
                axis.text.y = ggplot2::element_text(color = "grey20"),
                plot.title = ggplot2::element_text(face = "bold"),
                legend.title = ggplot2::element_text(size = 9),
                legend.text = ggplot2::element_text(size = 8)
            )
        heat_dir <- file.path(out_db, "heatmap")
        dir.create(heat_dir, recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave(file.path(heat_dir, paste0(safe_name, ".png")), heat_plot, width = 8, height = max(4, 0.35 * length(unique(heat_trait$Category))), dpi = 300)
    }

    heatmap_data <- if (length(heat_data)) dplyr::bind_rows(heat_data) else data.frame()
    heatmap_path <- file.path(out_db, "heatmap")

    invisible(list(
        database = db_name,
        sample_columns = sample_cols,
        traits = trait_cols,
        stacked_bar_data = bar_data,
        heatmap_data = heatmap_data,
        output_dir = out_db,
        heatmap = heatmap_path
    ))
}

#' Generate all four database-specific functional plot sets
#'
#' @param annotated_files Named character vector/list of four annotation files,
#'   or the directory produced by [annotate_functional_databases()].
#' @param output_dir Plot output directory.
#' @param trait_map Optional trait configuration.
#' @return Named list of plot results.
#' @export
plot_all_functional_profiles <- function(
    annotated_files,
    output_dir = "functional_plots",
    trait_map = NULL
) {
    files <- plankannot_resolve_annotation_files(annotated_files)
    results <- lapply(names(files), function(db) {
        plot_functional_profiles(files[[db]], db_name = db, output_dir = output_dir, trait_map = trait_map)
    })
    names(results) <- names(files)
    invisible(results)
}

#' @keywords internal
plankannot_default_trait_map <- function() {
    # Only traits explicitly selected in the stakeholder feedback are plotted.
    list(
        copepoda = c(
            "Functional group" = "Functional group",
            "Myelination" = "Myelination",
            "Spawning mode" = "Spawning mode",
            "Trophic group" = "Trophic group",
            "Feeding mode" = "Feeding mode"
        ),
        mixoplankton = c(
            "MFT" = "MFT (MIXOPLANKTON FUNCTIONAL CLASSIFICATION)|MFT (Functional Classification)|MFT"
        ),
        phytoplankton = c(
            "Shape" = "Shape",
            "Elongation" = "Elongation",
            "Group" = "Group"
        ),
        habs = c(
            "Definition" = "Definition"
        )
    )
}

#' @keywords internal
plankannot_normalize_name <- function(x) {
    x <- plankannot_safe_utf8(as.character(x))
    x <- tolower(x)
    gsub("[^a-z0-9]+", "", x)
}

#' @keywords internal
plankannot_resolve_traits <- function(columns, trait_map) {
    columns <- plankannot_safe_utf8(columns)
    out <- rep(NA_character_, length(trait_map))
    names(out) <- names(trait_map)
    norm_cols <- plankannot_normalize_name(columns)

    for (i in seq_along(trait_map)) {
        aliases <- strsplit(as.character(trait_map[[i]]), "\\|", fixed = FALSE)[[1]]
        aliases <- aliases[nzchar(trimws(aliases))]
        norm_aliases <- unique(plankannot_normalize_name(aliases))
        norm_aliases <- norm_aliases[nzchar(norm_aliases)]
        if (!length(norm_aliases)) next

        # Prefer the longest alias to avoid matching a generic suffix first.
        norm_aliases <- norm_aliases[order(nchar(norm_aliases), decreasing = TRUE)]

        # 1) Exact normalized match.
        hit <- which(norm_cols %in% norm_aliases)

        # 2) Database-prefixed/suffixed columns, e.g. copepoda_Functional group.
        if (!length(hit)) {
            hit <- which(vapply(norm_cols, function(n) {
                any(vapply(norm_aliases, function(a) {
                    identical(n, a) || (nchar(n) > nchar(a) && endsWith(n, a))
                }, logical(1)))
            }, logical(1)))
        }

        out[[i]] <- if (length(hit)) columns[hit[[1]]] else NA_character_
    }
    out
}

#' @keywords internal
plankannot_trait_long_data <- function(df, trait_col, trait_name, sample_cols) {
    categories <- plankannot_safe_utf8(as.character(df[[trait_col]]))
    categories[is.na(categories) | trimws(categories) == ""] <- NA_character_

    pieces <- lapply(sample_cols, function(s) {
        abundance <- suppressWarnings(as.numeric(df[[s]]))
        abundance[is.na(abundance)] <- 0
        tibble::tibble(
            Sample = s,
            Category = categories,
            Abundance = abundance
        ) |>
            dplyr::filter(!is.na(Category), is.finite(Abundance), Abundance > 0)
    })

    dplyr::bind_rows(pieces)
}

#' @keywords internal
plankannot_safe_utf8 <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
        bad <- !is.na(x)
        x[bad] <- iconv(x[bad], from = "", to = "UTF-8", sub = "")
        x[bad] <- gsub("[\\x00-\\x08\\x0B\\x0C\\x0E-\\x1F]", "", x[bad], perl = TRUE)
    }
    x
}

plankannot_read_annotation_input <- function(x) {
    if (is.data.frame(x)) return(x)
    if (!is.character(x) || length(x) != 1 || !file.exists(x)) stop("Annotation input not found.", call. = FALSE)
    ext <- tolower(tools::file_ext(x))
    if (ext == "xlsx") return(readxl::read_excel(x))
    if (ext %in% c("tsv", "txt")) return(readr::read_tsv(x, show_col_types = FALSE))
    stop("Unsupported annotation file: ", x, call. = FALSE)
}

#' @keywords internal
plankannot_resolve_annotation_files <- function(x) {
    if (is.character(x) && length(x) == 1 && dir.exists(x)) {
        out <- character()
        for (db in c("copepoda", "habs", "mixoplankton", "phytoplankton")) {
            candidates <- list.files(file.path(x, db), pattern = "_annotated\\.(xlsx|tsv)$", full.names = TRUE, ignore.case = TRUE)
            if (length(candidates)) out[db] <- candidates[[1]]
        }
        x <- out
    }
    if (is.null(names(x)) || any(!names(x) %in% c("copepoda", "habs", "mixoplankton", "phytoplankton"))) {
        stop("annotated_files must be a named vector/list using the four database names.", call. = FALSE)
    }
    x <- x[names(x) %in% c("copepoda", "habs", "mixoplankton", "phytoplankton")]
    missing <- names(x)[!file.exists(unname(x))]
    if (length(missing)) stop("Missing annotation files for: ", paste(missing, collapse = ", "), call. = FALSE)
    x
}
