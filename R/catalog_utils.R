# R/catalog_utils.R
library(data.table)
library(digest)
library(here)

catalog_paths <- function() {
  list(
    catalog_dir = here("data", "_catalog"),
    artifacts_csv = here("data", "_catalog", "catalogo_artefactos.csv"),
    runs_csv = here("data", "_catalog", "provenance_runs.csv")
  )
}

ensure_catalog_files <- function() {
  p <- catalog_paths()
  dir.create(p$catalog_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!file.exists(p$artifacts_csv)) {
    fwrite(
      data.table(
        dataset_id = character(),
        table_name = character(),
        version = character(),
        run_id = character(),
        artifact_type = character(),      # final_dataset | dictionary_ext | qc | report | staging | spec | master
        artifact_path = character(),
        file_ext = character(),
        n_rows = integer(),
        n_cols = integer(),
        file_hash = character(),
        created_at = character(),
        notes = character()
      ),
      p$artifacts_csv
    )
  }
  
  if (!file.exists(p$runs_csv)) {
    fwrite(
      data.table(
        run_id = character(),
        dataset_id = character(),
        version = character(),
        started_at = character(),
        finished_at = character(),
        status = character(),             # success | failed | running
        message = character()
      ),
      p$runs_csv
    )
  }
  
  invisible(TRUE)
}

# Normaliza schema del runs catalog (por si existe viejo con tipos incorrectos)
normalize_runs_schema <- function(runs_dt) {
  wanted <- c("run_id","dataset_id","version","started_at","finished_at","status","message")
  
  # asegurar columnas
  for (nm in wanted) if (!nm %in% names(runs_dt)) runs_dt[, (nm) := NA_character_]
  
  # castear a character
  for (nm in wanted) runs_dt[, (nm) := as.character(get(nm))]
  
  runs_dt[, ..wanted]
}

file_hash_md5 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "md5")
}

register_artifact <- function(dataset_id, table_name, version, run_id,
                              artifact_type, artifact_path,
                              n_rows = NA_integer_, n_cols = NA_integer_,
                              notes = NA_character_) {
  ensure_catalog_files()
  p <- catalog_paths()
  
  ext <- tools::file_ext(artifact_path)
  h <- file_hash_md5(artifact_path)
  
  row <- data.table(
    dataset_id = dataset_id,
    table_name = table_name,
    version = version,
    run_id = run_id,
    artifact_type = artifact_type,
    artifact_path = normalizePath(artifact_path, winslash = "/", mustWork = FALSE),
    file_ext = ext,
    n_rows = as.integer(n_rows),
    n_cols = as.integer(n_cols),
    file_hash = h,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    notes = notes
  )
  
  cat_dt <- fread(p$artifacts_csv)
  cat_dt <- rbind(cat_dt, row, fill = TRUE)
  fwrite(cat_dt, p$artifacts_csv)
  invisible(row)
}

register_run_start <- function(run_id, dataset_id, version) {
  ensure_catalog_files()
  p <- catalog_paths()
  
  runs <- fread(p$runs_csv)
  runs <- normalize_runs_schema(runs)
  
  runs <- rbind(
    runs,
    data.table(
      run_id = as.character(run_id),
      dataset_id = as.character(dataset_id),
      version = as.character(version),
      started_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      finished_at = NA_character_,
      status = "running",
      message = NA_character_
    ),
    fill = TRUE
  )
  
  runs <- normalize_runs_schema(runs)
  fwrite(runs, p$runs_csv)
}

register_run_finish <- function(run_id, status = c("success","failed"), message = NA_character_) {
  status <- match.arg(status)
  p <- catalog_paths()
  
  runs <- fread(p$runs_csv)
  runs <- normalize_runs_schema(runs)
  
  idx <- which(runs$run_id == as.character(run_id))
  if (length(idx) == 0) return(invisible(FALSE))
  
  runs[idx, `:=`(
    finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    status = status,
    message = as.character(message)
  )]
  
  runs <- normalize_runs_schema(runs)
  fwrite(runs, p$runs_csv)
  invisible(TRUE)
}
