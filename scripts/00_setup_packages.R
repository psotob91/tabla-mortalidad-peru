# scripts/00_setup_packages.R
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::restore(prompt = FALSE)
message("✅ Entorno restaurado desde renv.lock")
