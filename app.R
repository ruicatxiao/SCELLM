options(shiny.maxRequestSize = 1024^3)

library(shiny)
library(bslib)
library(ellmer)
library(Seurat)
library(shinychat)
library(base64enc)

source("seurat_llm_agent.R")

`%+%` <- function(lhs, rhs) paste0(lhs, rhs)

agent <- create_seurat_llm_agent(
  base_url = "http://127.0.0.1:1234/v1/",
  api_key  = "lm-studio",
  model    = "qwen/qwen3-vl-8b",
  data_root = ".",
  system_prompt_suffix =
    "You are running inside an R Shiny app.\n" %+%
    "- The Seurat object is loaded for you into '.active_seurat' when the user uploads an .rds file.\n" %+%
    "- Do NOT try to list or load files. Assume data is already loaded.\n" %+%
    "- You may call 'inspect_active_seurat', 'list_genes', and 'plot_seurat'.\n" %+%
    "- Use 'list_genes' to search for available genes before plotting.\n" %+%
    "- After plotting, always say exactly:\n" %+%
    "  'I have displayed the plot. Use the \"Analyze last plot\" button in the app if you want me to analyze it.'\n"
)

chat_client <- agent$chat
agent_state <- agent$state

ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),

  titlePanel("Seurat LLM Assistant"),

  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("Data"),
      fileInput(
        inputId = "seurat_rds",
        label   = "Upload Seurat .rds",
        accept  = ".rds"
      ),
      uiOutput("object_summary"),
      hr(),
      tags$div(
        style = "cursor: pointer;",
        imageOutput("plot_thumbnail", click = "plot_click")
      )
    ),

    mainPanel(
      width = 8,
      shinychat::chat_mod_ui(
        id = "llm_chat",
        fill = TRUE,
        messages = list(
          list(
            role = "assistant",
            content = paste(
              "Upload a Seurat .rds file on the left, or specify and S3 bucket and key"
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  current_plot <- reactiveVal(NULL)

  observe({
    invalidateLater(500)

    plot_path <- agent_state$last_plot_path

    if (!is.null(plot_path) && file.exists(plot_path)) {
      if (is.null(current_plot()) || current_plot() != plot_path) {
        current_plot(plot_path)
      }
    }
  })

  observeEvent(input$seurat_rds, {
    req(input$seurat_rds)

    rds_path <- input$seurat_rds$datapath

    seurat_object <- tryCatch(
      readRDS(rds_path),
      error = function(error_condition) {
        showNotification(
          paste("Failed to read RDS file:", error_condition$message),
          type = "error"
        )
        NULL
      }
    )

    validate(
      need(!is.null(seurat_object), "No valid object loaded from RDS."),
      need(inherits(seurat_object, "Seurat"), "Uploaded file is not a Seurat object.")
    )

    assign(
      x = agent_state$default_object_name,
      value = seurat_object,
      envir = .GlobalEnv
    )

    agent_state$last_plot_path <- NULL
    current_plot(NULL)

    showNotification(
      paste(
        "Seurat object loaded with",
        ncol(seurat_object), "cells and",
        nrow(seurat_object), "features."
      ),
      type = "message"
    )
  })

  output$object_summary <- renderUI({
    if (!exists(agent_state$default_object_name, envir = .GlobalEnv)) {
      return(helpText("No Seurat object loaded yet."))
    }

    seurat_object <- get(agent_state$default_object_name, envir = .GlobalEnv)

    tagList(
      h5("Active Seurat object"),
      tags$ul(
        tags$li(sprintf("Cells: %s", ncol(seurat_object))),
        tags$li(sprintf("Features: %s", nrow(seurat_object))),
        tags$li(sprintf("Assays: %s", paste(names(seurat_object@assays), collapse = ", "))),
        tags$li(sprintf("Reductions: %s", paste(names(seurat_object@reductions), collapse = ", ")))
      )
    )
  })

  shinychat::chat_mod_server("llm_chat", chat_client)

  observeEvent(input$analyze_plot, {
    plot_path <- current_plot()

    if (is.null(plot_path) || !file.exists(plot_path)) {
      showNotification(
        "No plot found. Ask the agent to generate a plot first.",
        type = "warning"
      )
      return(NULL)
    }

    analysis_stream <- chat_client$stream(
      "Analyze the biological patterns and cluster-level differences in this plot.",
      ellmer::content_image_file(plot_path)
    )

    shinychat::chat_append(
      id = "llm_chat",
      response = analysis_stream,
      role = "assistant"
    )
  })

  output$plot_thumbnail <- renderImage({
    plot_path <- current_plot()
    req(plot_path, file.exists(plot_path))

    list(
      src         = plot_path,
      contentType = "image/png",
      alt         = "LLM generated plot thumbnail",
      width       = "100%",
      style       = "border: 1px solid #ddd; border-radius: 4px; max-width: 100%;"
    )
  }, deleteFile = FALSE)

  observeEvent(input$plot_click, {
    plot_path <- current_plot()
    req(plot_path, file.exists(plot_path))

    img_data <- base64enc::base64encode(plot_path)

    showModal(
      modalDialog(
        title = "Plot",
        easyClose = TRUE,
        size = "l",
        footer = NULL,
        tags$img(
          src = paste0("data:image/png;base64,", img_data),
          style = "width: 100%; height: auto;"
        )
      )
    )
  })

}

shinyApp(ui, server)
