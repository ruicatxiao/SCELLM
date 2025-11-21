source("seurat_llm_agent.R")

agent <- create_seurat_llm_agent(
  model = "qwen/qwen3-vl-8b",
  base_url = "http://127.0.0.1:1234/v1/",
  api_key = "lm-studio",
  data_root = "."
)

agent$chat$chat("List RDS files here and load the one called all.rds.")
agent$chat$chat("Plot the UMAP colored by clusters.")
agent$chat$chat("Plot gene cgd7-1900 expression.")
agent$analyze_last_plot("Analyze the generated plot.")
