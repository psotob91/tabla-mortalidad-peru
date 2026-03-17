# R/spec_utils.R
library(data.table)
library(yaml)

read_spec <- function(path) yaml::read_yaml(path)

assert_cols_and_types <- function(dt, required_columns) {
  missing_cols <- setdiff(names(required_columns), names(dt))
  if (length(missing_cols) > 0)
    stop("Faltan columnas obligatorias: ", paste(missing_cols, collapse = ", "))
  
  for (nm in names(required_columns)) {
    want <- required_columns[[nm]]
    x <- dt[[nm]]
    
    ok <- switch(
      want,
      integer   = is.integer(x) || (is.numeric(x) && all(is.na(x) | x == as.integer(x))),
      numeric   = is.numeric(x),
      character = is.character(x),
      logical   = is.logical(x),
      TRUE
    )
    if (!ok) stop("Tipo inválido en '", nm, "'. Esperado: ", want)
  }
  invisible(TRUE)
}

assert_unique_pk <- function(dt, pk) {
  dups <- dt[, .N, by = pk][N > 1]
  if (nrow(dups) > 0) stop("Duplicados en PK lógica: ", paste(pk, collapse = ", "))
  invisible(TRUE)
}

assert_range <- function(dt, col, min = NULL, max = NULL, allow_na = TRUE) {
  x <- dt[[col]]
  if (!allow_na && anyNA(x)) stop("NA no permitido en: ", col)
  if (!is.null(min) && any(x < min, na.rm = TRUE)) stop("Valores < ", min, " en: ", col)
  if (!is.null(max) && any(x > max, na.rm = TRUE)) stop("Valores > ", max, " en: ", col)
  invisible(TRUE)
}

assert_allowed <- function(dt, col, allowed_values, allow_na = TRUE) {
  x <- dt[[col]]
  if (!allow_na && anyNA(x)) stop("NA no permitido en: ", col)
  bad <- x[!(x %in% allowed_values) & !is.na(x)]
  if (length(bad) > 0) stop("Valores no permitidos en ", col, ": ", paste(unique(bad), collapse = ", "))
  invisible(TRUE)
}

# ✅ NUEVO: Validador genérico basado 100% en spec
validate_by_spec <- function(dt, spec) {
  dt <- as.data.table(dt)
  
  assert_cols_and_types(dt, spec$required_columns)
  assert_unique_pk(dt, spec$primary_key)
  
  cts <- spec$constraints
  if (is.null(cts)) return(invisible(TRUE))
  
  for (col in intersect(names(cts), names(dt))) {
    rule <- cts[[col]]
    
    allow_na <- TRUE
    if (!is.null(rule$allow_na)) allow_na <- isTRUE(rule$allow_na)
    
    # allowed_values tiene prioridad
    if (!is.null(rule$allowed_values)) {
      assert_allowed(dt, col, rule$allowed_values, allow_na = allow_na)
    }
    
    # rango (si existe)
    minv <- if (!is.null(rule$min)) rule$min else NULL
    maxv <- if (!is.null(rule$max)) rule$max else NULL
    if (!is.null(minv) || !is.null(maxv) || !allow_na) {
      assert_range(dt, col, min = minv, max = maxv, allow_na = allow_na)
    }
  }
  
  invisible(TRUE)
}

# Mantener el antiguo para proyectos INEI (sin romper nada)
validate_population_result <- function(dt, spec) {
  dt <- as.data.table(dt)
  
  assert_cols_and_types(dt, spec$required_columns)
  assert_unique_pk(dt, spec$primary_key)
  
  cts <- spec$constraints
  
  assert_range(dt, "year_id", cts$year_id$min, cts$year_id$max, cts$year_id$allow_na)
  assert_range(dt, "age",     cts$age$min,     cts$age$max,     cts$age$allow_na)
  
  assert_allowed(dt, "sex_id", cts$sex_id$allowed_values, cts$sex_id$allow_na)
  assert_allowed(dt, "location_id", cts$location_id$allowed_values, cts$location_id$allow_na)
  
  assert_range(dt, "population", cts$population$min, NULL, cts$population$allow_na)
  
  invisible(TRUE)
}