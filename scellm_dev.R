# Additional dev for LLM calling on fucntions with ellmer
library(ellmer)
library(Seurat)
library(ggplot2)
library(grid)

# LLM SETUP

# Initialize Chat
chat <- chat_vllm(
  base_url = "http://127.0.0.1:1234/v1/",
  api_key = "lm-studio",
  model = "qwen/qwen3-vl-8b", 
  system_prompt = "You are a Seurat Analysis Agent. 
  1. When asked to plot, use 'plot_seurat'. 
  2. You cannot see the plot immediately. 
  3. After plotting, tell the user: 'I have displayed the plot. If you want me to analyze it, please run analyze_last_plot()'."
)


# Global variable to track the last generated plot for easy "follow-up" analysis
.last_plot_path <- NULL


# TOOL 1: List RDS Files
list_rds_files_tool <- tool(
  function(path = ".") {
    if (is.list(path)) path <- unlist(path)[1]
    
    # Check if path exists
    if (!dir.exists(path)) return(paste("Directory not found:", path))
    
    files <- list.files(path, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
    
    if (length(files) == 0) {
      return("No RDS files found in this directory.")
    }
    
    info <- file.info(files)
    file_list <- paste(
      basename(files), 
      paste0("(", round(info$size / 1024^2, 2), " MB)"),
      collapse = "\n"
    )
    return(paste("Found the following files:\n", file_list))
  },
  name = "list_rds_files",
  description = "List .rds files in the directory. Use this to find the filename before loading.",
  arguments = list(
    path = type_string("Directory path. Default: '.'", required = FALSE)
  )
)

# TOOL 2: Load Seurat Object
load_seurat_object_tool <- tool(
  function(file_path, object_name = ".active_seurat") {
    file_path <- if (is.list(file_path)) unlist(file_path)[1] else file_path
    object_name <- if (is.list(object_name)) unlist(object_name)[1] else object_name
    
    # 1. Check if file exists
    if (!file.exists(file_path)) {
      # Try to help the LLM if it missed the extension or path
      if (file.exists(paste0(file_path, ".rds"))) file_path <- paste0(file_path, ".rds")
      else stop(paste("File not found:", file_path))
    }
    
    cat("Loading Seurat object...\n")
    
    # 2. Load and Validate
    tryCatch({
      obj <- readRDS(file_path)
      if (!inherits(obj, "Seurat")) stop("File is not a Seurat object.")
      
      # 3. Assign to Global Env
      assign(object_name, obj, envir = .GlobalEnv)
      
      # 4. Return Summary
      list(
        status = "Loaded",
        cells = ncol(obj),
        genes = nrow(obj),
        assays = names(obj@assays),
        active.assay = DefaultAssay(obj),
        reductions = names(obj@reductions),
        message = "Object loaded. You can now inspect metadata or plot."
      )
    }, error = function(e) {
      paste("Error loading file:", e$message)
    })
  },
  name = "load_seurat_object",
  description = "Load an RDS file into memory. Must be done before any plotting.",
  arguments = list(
    file_path = type_string("Path to the .rds file", required = TRUE),
    object_name = type_string("Internal var name (default: .active_seurat)", required = FALSE)
  )
)

# TOOL 3: Inspect Metadata (Fuzzy Match)
inspect_active_seurat_tool <- tool(
  function(object_name = ".active_seurat") {
    # 1. Input sanitization
    object_name <- if (is.list(object_name)) unlist(object_name)[1] else object_name
    
    # 2. Check for cached object
    if (!exists(object_name, envir = .GlobalEnv)) {
      return(paste("❌ No object found named", object_name, "- Please load a file first."))
    }
    
    obj <- get(object_name, envir = .GlobalEnv)
    
    # 3. Smart Metadata Summary (The upgrade)
    # Instead of just counting levels, we show the first few to help the LLM.
    meta_summary <- lapply(names(obj@meta.data), function(col) {
      val <- obj@meta.data[[col]]
      
      if (is.numeric(val)) {
        # For numbers: Show Range
        paste0(col, " (Numeric: ", round(min(val, na.rm=T), 2), " to ", round(max(val, na.rm=T), 2), ")")
      } else {
        # For categories: Show specific levels
        levs <- unique(as.character(val))
        if (length(levs) > 20) {
          # Too many? Truncate
          paste0(col, " (High cardinality: ", length(levs), " unique values)")
        } else {
          # Show the actual names (Crucial for LLM to know what to plot!)
          paste0(col, " (Levels: ", paste(levs, collapse = ", "), ")")
        }
      }
    })
    
    # 4. Return comprehensive list
    list(
      status = "✅ Object Inspected",
      num_cells = ncol(obj),
      num_features = nrow(obj),
      active_assay = DefaultAssay(obj),
      available_assays = names(obj@assays),
      reductions = names(obj@reductions),
      # This returns your full list of metadata columns + their content
      metadata_details = unlist(meta_summary), 
      # Explicitly grab active identity levels (Cluster names)
      active_identity_levels = levels(Idents(obj))
    )
  },
  name = "inspect_active_seurat",
  description = "Inspects the loaded Seurat object. Returns cell counts, active assay, specific metadata column names, and cluster identities.",
  # FIXED: Arguments now match the function exactly
  arguments = list(
    object_name = type_string("Internal variable name. Default: '.active_seurat'", required = FALSE)
  )
)

# TOOL 4: Universal Plotter
plot_seurat_tool <- tool(
  function(feature, type = "auto", object_name = ".active_seurat") {
    # Standardize inputs from LLM
    feature <- if (is.list(feature)) unlist(feature)[1] else feature
    type <- if (is.list(type)) unlist(type)[1] else type
    object_name <- if (is.list(object_name)) unlist(object_name)[1] else object_name
    
    if (!exists(object_name, envir = .GlobalEnv)) stop("No object loaded.")
    obj <- get(object_name, envir = .GlobalEnv)
    
    target <- feature
    mode <- "unknown"
    
    if (feature %in% colnames(obj@meta.data)) {
      mode <- "meta"
    } else {
      # Fuzzy match metadata
      matches <- grep(feature, colnames(obj@meta.data), ignore.case = TRUE, value = TRUE)
      if (length(matches) > 0) {
        target <- matches[1]
        mode <- "meta"
      }
    }
    
    if (mode == "unknown") {
      # Check against features in the object
      if (feature %in% rownames(obj)) {
        mode <- "gene"
      } else {
        stop(paste("Could not find '", feature, "' in metadata or genes."))
      }
    }
    
    p <- NULL
    if (mode == "meta") {
      is_numeric <- is.numeric(obj@meta.data[[target]]) && length(unique(obj@meta.data[[target]])) > 20
      if (is_numeric) {
        p <- FeaturePlot(obj, features = target, order = TRUE) + 
             scale_color_viridis_c()
      } else {
        p <- DimPlot(obj, group.by = target, label = TRUE, label.size = 5)
      }
    } else {
      # Gene Expression
      p <- FeaturePlot(obj, features = target, order = TRUE, pt.size = 0.8) + 
           scale_color_viridis_c(option = "magma")
    }
    
    p <- p + 
      ggtitle(paste(target, ifelse(mode=="gene", "Expression", ""))) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        legend.position = "bottom"
      )
    

    temp_file <- tempfile(pattern = "llm_plot_", fileext = ".png")
    ggsave(temp_file, p, width = 8, height = 8, bg = "white")
    
    assign(".last_plot_path", temp_file, envir = .GlobalEnv)
    
    print(p) 
    
    return(list(
      status = "Plot Displayed",
      file_path = temp_file,
      feature_plotted = target,
      mode = mode,
      instruction = "The plot is visible to the user. If you need to analyze it, ask the user to 'Use the vision helper'."
    ))
  },
  name = "plot_seurat",
  description = "Plots a gene OR metadata column. Auto-detects which one it is.",
  arguments = list(
    feature = type_string("Gene name (e.g. 'CD14') or Metadata column (e.g. 'clusters')", required = TRUE),
    type = type_string("Manual override type (optional)", required = FALSE),
    object_name = type_string("Internal variable name (optional)", required = FALSE)
  )
)




chat$register_tool(list_rds_files_tool)
chat$register_tool(load_seurat_object_tool)
chat$register_tool(inspect_active_seurat_tool)
chat$register_tool(plot_seurat_tool)

# HELPER: The Vision Bridge

analyze_last_plot <- function(question = "Describe the biological patterns in this plot.") {
  if (is.null(.last_plot_path) || !file.exists(.last_plot_path)) {
    stop("No plot has been generated yet!")
  }
  
  cat("Sending plot to LLM for analysis...\n")
  
  response <- chat$chat(question, content_image_file(.last_plot_path))
  
  return(response)
}



chat$chat("List rds files here and load the one called all.rds.")
# > LLM calls list_rds_files
# > LLM calls load_seurat_object

chat$chat("Plot the UMAP colored by clusters.")
# > LLM calls plot_seurat(feature="seurat_clusters")
# > RStudio displays the plot in the pane.
# > LLM responds: "I've plotted the clusters. Run analyze_last_plot() to have me look at it."

chat$chat("Plot gene cgd7-1900 expression.")
# > RStudio updates plot.

analyze_last_plot("Analyze the generated plot")
# > Function grabs .last_plot_path automatically.
# > LLM sees the image and responds.