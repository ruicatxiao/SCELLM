create_seurat_llm_agent <- function(
  base_url = Sys.getenv("LLM_BASE_URL", "http://127.0.0.1:1234/v1/"),
  api_key = Sys.getenv("LLM_API_KEY", "lm-studio"),
  model = Sys.getenv("LLM_MODEL", "qwen/qwen3-vl-8b"),
  data_root = ".",
  default_object_name = ".active_seurat",
  system_prompt_suffix = NULL
) {
  required_packages <- c("ellmer", "Seurat", "ggplot2", "grid")

  for (package_name in required_packages) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      stop(
        sprintf("Package '%s' is required but not installed.", package_name),
        call. = FALSE
      )
    }
  }

  state_environment <- new.env(parent = emptyenv())
  state_environment$last_plot_path <- NULL
  state_environment$default_object_name <- default_object_name
  state_environment$data_root <- normalizePath(data_root, mustWork = FALSE)

  system_prompt <- paste(
    "You are a Seurat Analysis Agent.",
    "",
    "Tools and behavior:",
    "1. To discover available RDS files, call the tool 'list_rds_files'.",
    "2. To load a Seurat object, call 'load_seurat_object' before plotting.",
    "3. To inspect metadata and cluster identities, call 'inspect_active_seurat'.",
    "4. To search for genes, call 'list_genes' with an optional pattern (e.g. 'CD', 'IL', 'MT-').",
    "5. When asked to plot, always use 'plot_seurat'.",
    "6. You cannot see the plot immediately.",
    "7. After plotting you MUST respond with:",
    "   'I have displayed the plot. If you want me to analyze it, please run analyze_last_plot()'.",
    "",
    if (!is.null(system_prompt_suffix)) system_prompt_suffix else "",
    sep = "\n"
  )

  chat <- ellmer::chat_vllm(
    base_url = base_url,
    api_key = api_key,
    model = model,
    system_prompt = system_prompt
  )

  # TOOL 1: List RDS files
  list_rds_files_tool <- ellmer::tool(
    function(directory_path = state_environment$data_root,
             pattern = "\\.rds$") {

      if (is.list(directory_path)) directory_path <- unlist(directory_path)[1]
      if (is.list(pattern)) pattern <- unlist(pattern)[1]

      if (!dir.exists(directory_path)) {
        return(list(
          status = "error",
          message = paste("Directory not found:", directory_path),
          files = list()
        ))
      }

      file_paths <- list.files(
        path = directory_path,
        pattern = pattern,
        full.names = TRUE,
        ignore.case = TRUE
      )

      if (length(file_paths) == 0) {
        return(list(
          status = "ok",
          message = "No RDS files found in this directory.",
          files = list()
        ))
      }

      file_info <- file.info(file_paths)

      files_table <- data.frame(
        file_name = basename(file_paths),
        file_path = normalizePath(file_paths, winslash = "/", mustWork = FALSE),
        size_megabytes = round(file_info$size / 1024^2, 2),
        modification_time = file_info$mtime,
        stringsAsFactors = FALSE
      )

      summary_text <- paste(
        apply(
          files_table,
          1,
          function(row_values) {
            paste0(
              row_values[["file_name"]],
              " (",
              row_values[["size_megabytes"]],
              " MB, modified ",
              row_values[["modification_time"]],
              ")"
            )
          }
        ),
        collapse = "\n"
      )

      list(
        status = "ok",
        message = paste("Found the following RDS files:\n", summary_text),
        files = files_table
      )
    },
    name = "list_rds_files",
    description = paste(
      "List .rds files in a directory. Use this to discover valid RDS files",
      "before calling 'load_seurat_object'."
    ),
    arguments = list(
      directory_path = ellmer::type_string(
        "Directory path. Default: the agent's data_root.",
        required = FALSE
      ),
      pattern = ellmer::type_string(
        "Regular expression pattern used to match RDS files (default '\\\\.rds$').",
        required = FALSE
      )
    )
  )

  # TOOL 2: Load Seurat Object
  load_seurat_object_tool <- ellmer::tool(
    function(file_path, object_name = state_environment$default_object_name) {
      file_path <- if (is.list(file_path)) unlist(file_path)[1] else file_path
      object_name <- if (is.list(object_name)) unlist(object_name)[1] else object_name

      if (!file.exists(file_path)) {
        candidate_file_path <- paste0(file_path, ".rds")
        if (file.exists(candidate_file_path)) {
          file_path <- candidate_file_path
        } else {
          stop(paste("File not found:", file_path))
        }
      }

      message("Loading Seurat object from: ", normalizePath(file_path))

      tryCatch({
        seurat_object <- readRDS(file_path)

        if (!inherits(seurat_object, "Seurat")) {
          stop("Loaded object is not a Seurat object.")
        }

        assign(
          x = object_name,
          value = seurat_object,
          envir = .GlobalEnv
        )

        list(
          status = "ok",
          message = paste(
            "Seurat object loaded into .GlobalEnv as",
            shQuote(object_name)
          ),
          object_name = object_name,
          cells = ncol(seurat_object),
          features = nrow(seurat_object),
          assays = names(seurat_object@assays),
          active_assay = Seurat::DefaultAssay(seurat_object),
          reductions = names(seurat_object@reductions)
        )
      }, error = function(error_condition) {
        list(
          status = "error",
          message = paste("Error loading Seurat object:", error_condition$message)
        )
      })
    },
    name = "load_seurat_object",
    description = paste(
      "Load an RDS file containing a Seurat object into memory.",
      "Must be called before any plotting or metadata inspection."
    ),
    arguments = list(
      file_path = ellmer::type_string(
        "Path to the .rds file containing a Seurat object.",
        required = TRUE
      ),
      object_name = ellmer::type_string(
        "Symbol name under which the Seurat object will be stored in .GlobalEnv.",
        required = FALSE
      )
    )
  )

  # TOOL 3: Inspect Metadata (Fuzzy Match)
  inspect_active_seurat_tool <- ellmer::tool(
    function(object_name = state_environment$default_object_name) {
      object_name <- if (is.list(object_name)) unlist(object_name)[1] else object_name

      if (!exists(object_name, envir = .GlobalEnv)) {
        return(list(
          status = "error",
          message = paste(
            "No object found named",
            object_name,
            "- call 'load_seurat_object' first."
          )
        ))
      }

      seurat_object <- get(object_name, envir = .GlobalEnv)

      metadata_column_names <- names(seurat_object@meta.data)

      metadata_summary <- lapply(metadata_column_names, function(column_name) {
        column_values <- seurat_object@meta.data[[column_name]]

        if (is.numeric(column_values)) {
          paste0(
            column_name,
            " (Numeric: ",
            round(min(column_values, na.rm = TRUE), 2),
            " to ",
            round(max(column_values, na.rm = TRUE), 2),
            ")"
          )
        } else {
          level_values <- unique(as.character(column_values))
          if (length(level_values) > 20) {
            paste0(
              column_name,
              " (High cardinality: ",
              length(level_values),
              " unique values)"
            )
          } else {
            paste0(
              column_name,
              " (Levels: ",
              paste(level_values, collapse = ", "),
              ")"
            )
          }
        }
      })

      list(
        status = "ok",
        message = "Seurat object inspected successfully.",
        object_name = object_name,
        num_cells = ncol(seurat_object),
        num_features = nrow(seurat_object),
        active_assay = Seurat::DefaultAssay(seurat_object),
        available_assays = names(seurat_object@assays),
        reductions = names(seurat_object@reductions),
        metadata_details = unlist(metadata_summary),
        active_identity_levels = levels(Seurat::Idents(seurat_object))
      )
    },
    name = "inspect_active_seurat",
    description = paste(
      "Inspect the currently loaded Seurat object.",
      "Returns cell and feature counts, assay information, reductions,",
      "metadata column summaries and cluster identities."
    ),
    arguments = list(
      object_name = ellmer::type_string(
        "Name of the Seurat object in .GlobalEnv. Default: agent's default_object_name.",
        required = FALSE
      )
    )
  )

  # TOOL 4: List Available Genes
  list_genes_tool <- ellmer::tool(
    function(pattern = NULL,
             object_name = state_environment$default_object_name,
             max_results = 100) {

      pattern <- if (is.list(pattern)) unlist(pattern)[1] else pattern
      object_name <- if (is.list(object_name)) unlist(object_name)[1] else object_name
      max_results <- if (is.list(max_results)) unlist(max_results)[1] else as.numeric(max_results)

      if (!exists(object_name, envir = .GlobalEnv)) {
        return(list(
          status = "error",
          message = paste(
            "No object found named",
            object_name,
            "- call 'load_seurat_object' first."
          )
        ))
      }

      seurat_object <- get(object_name, envir = .GlobalEnv)
      all_genes <- rownames(seurat_object)

      if (!is.null(pattern) && nchar(pattern) > 0) {
        matched_genes <- grep(
          pattern = pattern,
          x = all_genes,
          ignore.case = TRUE,
          value = TRUE
        )

        if (length(matched_genes) == 0) {
          return(list(
            status = "ok",
            message = paste(
              "No genes matched pattern:",
              shQuote(pattern),
              "- Total genes available:",
              length(all_genes)
            ),
            matched_genes = list(),
            total_genes = length(all_genes)
          ))
        }

        display_genes <- head(matched_genes, max_results)
        more_available <- length(matched_genes) > max_results

        list(
          status = "ok",
          message = paste(
            "Found",
            length(matched_genes),
            "genes matching pattern",
            shQuote(pattern),
            if (more_available) paste("(showing first", max_results, ")") else ""
          ),
          matched_genes = display_genes,
          total_matched = length(matched_genes),
          total_genes = length(all_genes),
          more_available = more_available
        )
      } else {
        display_genes <- head(all_genes, max_results)
        more_available <- length(all_genes) > max_results

        list(
          status = "ok",
          message = paste(
            "Total genes available:",
            length(all_genes),
            if (more_available) paste("(showing first", max_results, ")") else ""
          ),
          genes = display_genes,
          total_genes = length(all_genes),
          more_available = more_available
        )
      }
    },
    name = "list_genes",
    description = paste(
      "List all available genes in the Seurat object.",
      "Optionally filter by pattern (case-insensitive).",
      "Use this to discover valid gene names before plotting."
    ),
    arguments = list(
      pattern = ellmer::type_string(
        "Optional pattern to search for genes (case-insensitive, e.g. 'CD3', 'IL', 'MT-').",
        required = FALSE
      ),
      object_name = ellmer::type_string(
        "Name of the Seurat object in .GlobalEnv. Default: agent's default_object_name.",
        required = FALSE
      ),
      max_results = ellmer::type_integer(
        "Maximum number of genes to return (default 100).",
        required = FALSE
      )
    )
  )

  # TOOL 5: Universal Seurat Plotter
  plot_seurat_tool <- ellmer::tool(
    function(feature,
             type = "auto",
             object_name = state_environment$default_object_name) {

      feature <- if (is.list(feature)) unlist(feature)[1] else feature
      type <- if (is.list(type)) unlist(type)[1] else type
      object_name <- if (is.list(object_name)) unlist(object_name)[1] else object_name

      if (!exists(object_name, envir = .GlobalEnv)) {
        stop("No Seurat object loaded. Call 'load_seurat_object' first.")
      }

      seurat_object <- get(object_name, envir = .GlobalEnv)

      target_feature <- feature
      feature_mode <- "unknown"

      if (feature %in% colnames(seurat_object@meta.data)) {
        feature_mode <- "metadata"
      } else {
        matched_metadata_columns <- grep(
          pattern = feature,
          x = colnames(seurat_object@meta.data),
          ignore.case = TRUE,
          value = TRUE
        )
        if (length(matched_metadata_columns) > 0) {
          target_feature <- matched_metadata_columns[1]
          feature_mode <- "metadata"
        }
      }

      if (feature_mode == "unknown") {
        if (feature %in% rownames(seurat_object)) {
          feature_mode <- "gene"
        } else {
          stop(
            paste(
              "Could not find feature",
              shQuote(feature),
              "in metadata or genes."
            )
          )
        }
      }

      plot_object <- NULL

      if (feature_mode == "metadata") {
        column_values <- seurat_object@meta.data[[target_feature]]
        is_numeric_metadata <- is.numeric(column_values) &&
          length(unique(column_values)) > 20

        if (is_numeric_metadata) {
          plot_object <- Seurat::FeaturePlot(
            seurat_object,
            features = target_feature,
            order = TRUE
          ) +
            ggplot2::scale_color_viridis_c()
        } else {
          plot_object <- Seurat::DimPlot(
            seurat_object,
            group.by = target_feature,
            label = TRUE,
            label.size = 5
          )
        }
      } else {
        plot_object <- Seurat::FeaturePlot(
          seurat_object,
          features = target_feature,
          order = TRUE,
          pt.size = 0.8
        ) +
          ggplot2::scale_color_viridis_c(option = "magma")
      }

      plot_object <- plot_object +
        ggplot2::ggtitle(
          paste(
            target_feature,
            ifelse(feature_mode == "gene", "Expression", "")
          )
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 18, face = "bold"),
          legend.position = "bottom"
        )

      temp_file <- tempfile(pattern = "llm_plot_", fileext = ".png")

      ggplot2::ggsave(
        filename = temp_file,
        plot     = plot_object,
        width    = 4,
        height   = 4,
        units    = "in",
        dpi      = 120,
        bg       = "white"
      )

      state_environment$last_plot_path <- temp_file

      print(plot_object)

      list(
        status = "ok",
        message = paste(
          "Plot generated and saved to",
          normalizePath(temp_file)
        ),
        file_path = temp_file,
        feature_plotted = target_feature,
        mode = feature_mode,
        instruction = paste(
          "The plot is visible to the user.",
          "Tell the user:",
          "'I have displayed the plot. If you want me to analyze it, please run analyze_last_plot()'."
        )
      )
    },
    name = "plot_seurat",
    description = paste(
      "Plot a gene or metadata column on the active Seurat object.",
      "Automatically detects whether the requested feature is a gene or metadata.",
      "For numeric metadata, uses FeaturePlot; for categorical metadata, uses DimPlot."
    ),
    arguments = list(
      feature = ellmer::type_string(
        "Gene name (e.g. 'CD14') or metadata column (e.g. 'seurat_clusters').",
        required = TRUE
      ),
      type = ellmer::type_string(
        "Manual override for plot type (currently advisory only).",
        required = FALSE
      ),
      object_name = ellmer::type_string(
        "Name of the Seurat object in .GlobalEnv. Default: agent's default_object_name.",
        required = FALSE
      )
    )
  )

  chat$register_tool(list_rds_files_tool)
  chat$register_tool(load_seurat_object_tool)
  chat$register_tool(inspect_active_seurat_tool)
  chat$register_tool(list_genes_tool)
  chat$register_tool(plot_seurat_tool)

  # HELPER: The Vision Bridge

  analyze_last_plot <- function(
    question = "Describe the biological patterns in this plot."
  ) {
    if (is.null(state_environment$last_plot_path) ||
        !file.exists(state_environment$last_plot_path)) {
      stop("No plot has been generated yet. Call 'plot_seurat' via the agent first.")
    }

    message("Sending plot to LLM for analysis...")

    response <- chat$chat(
      question,
      ellmer::content_image_file(state_environment$last_plot_path)
    )

    response
  }

  structure(
    list(
      chat = chat,
      analyze_last_plot = analyze_last_plot,
      state = state_environment
    ),
    class = "SeuratLlmAgent"
  )
}
