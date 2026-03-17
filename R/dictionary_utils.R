# R/dictionary_utils.R
library(data.table)
library(yaml)

`%||%` <- function(x, y) if (is.null(x)) y else x
read_spec <- function(path) yaml::read_yaml(path)

# build labels for location_id from maestro_location_dept.csv (si existe)
location_allowed_values_label <- function(config_dir) {
  fp <- file.path(config_dir, "maestro_location_dept.csv")
  if (!file.exists(fp)) return(NA_character_)
  m <- fread(fp)
  # espera: location_id, location_name, ubigeo_dept_str, level
  m <- m[order(location_id)]
  paste0(m$ubigeo_dept_str, "=", m$location_name, collapse = "; ")
}

dict_from_spec <- function(spec, dataset_version = "v1.0.0", run_id = NA_character_, config_dir = "config") {
  table_name <- spec$table_name
  dataset_id <- spec$dataset_id
  req_cols <- names(spec$required_columns)
  pk <- spec$primary_key
  cts <- spec$constraints
  dict_meta <- spec$dictionary
  
  get_meta <- function(col, field, default = NA_character_) {
    if (!is.null(dict_meta) && !is.null(dict_meta[[col]]) && !is.null(dict_meta[[col]][[field]])) {
      dict_meta[[col]][[field]]
    } else default
  }
  
  # allowed_values_label special-case: location_id
  loc_avl <- location_allowed_values_label(config_dir)
  
  dt <- rbindlist(lapply(req_cols, function(col) {
    dtype <- spec$required_columns[[col]]
    allow_na <- cts[[col]]$allow_na %||% TRUE
    minv <- cts[[col]]$min %||% NA
    maxv <- cts[[col]]$max %||% NA
    
    allowed_values <- if (!is.null(cts[[col]]$allowed_values)) paste(cts[[col]]$allowed_values, collapse = "|") else NA_character_
    
    av_label <- get_meta(col, "allowed_values_label", NA_character_)
    if (col == "location_id" && !is.na(loc_avl)) av_label <- loc_avl
    
    data.table(
      dataset_id = dataset_id,
      version = dataset_version,
      run_id = run_id,
      table_name = table_name,
      column_name = col,
      label = get_meta(col, "label", col),
      data_type = dtype,
      is_required = TRUE,
      is_key = col %in% pk,
      allow_na = allow_na,
      value_min = minv,
      value_max = maxv,
      allowed_values = allowed_values,
      allowed_values_label = av_label,
      units = get_meta(col, "units", NA_character_),
      omop_domain_id = get_meta(col, "omop_domain_id", NA_character_),
      is_omop_concept_id = as.logical(get_meta(col, "is_omop_concept_id", FALSE)),
      notes = NA_character_
    )
  }))
  
  meta_rows <- data.table(
    dataset_id = dataset_id,
    version = dataset_version,
    run_id = run_id,
    table_name = table_name,
    column_name = c("META:age_definition", "META:year_definition", "META:population_units",
                    "META:location_definition", "META:extrapolation_method"),
    label = c("Definición de edad", "Definición de año", "Unidades de población",
              "Definición de ubicación", "Método de extrapolación"),
    data_type = "character",
    is_required = FALSE,
    is_key = FALSE,
    allow_na = TRUE,
    value_min = NA,
    value_max = NA,
    allowed_values = NA_character_,
    allowed_values_label = NA_character_,
    units = NA_character_,
    omop_domain_id = "metadata",
    is_omop_concept_id = FALSE,
    notes = c(
      spec$policy$age_definition %||% NA_character_,
      spec$policy$year_definition %||% NA_character_,
      spec$policy$population_units %||% NA_character_,
      spec$policy$location_definition %||% NA_character_,
      spec$policy$extrapolation$method %||% NA_character_
    )
  )
  
  rbind(dt, meta_rows, fill = TRUE)
}

enrich_dict_with_stats <- function(dict, data) {
  dt <- as.data.table(data)
  dict2 <- copy(dict)
  
  cols <- dict2[!startsWith(column_name, "META:"), unique(column_name)]
  cols <- intersect(cols, names(dt))
  
  stats <- rbindlist(lapply(cols, function(col) {
    x <- dt[[col]]
    ex <- unique(x[!is.na(x)])
    ex <- head(ex, 5)
    ex_str <- paste(ex, collapse = "|")
    
    data.table(
      column_name = col,
      n_missing = sum(is.na(x)),
      pct_missing = round(100 * mean(is.na(x)), 4),
      observed_min = suppressWarnings(if (is.numeric(x) || is.integer(x)) min(x, na.rm = TRUE) else NA_real_),
      observed_max = suppressWarnings(if (is.numeric(x) || is.integer(x)) max(x, na.rm = TRUE) else NA_real_),
      example_values = ex_str
    )
  }), fill = TRUE)
  
  merge(dict2, stats, by = "column_name", all.x = TRUE, sort = FALSE)
}
