# # R/qc_life_table_utils.R
# library(data.table)
# 
# qc_life_table_basic <- function(dt) {
#   x <- as.data.table(dt)
#   
#   out <- list()
#   
#   # 1) Rango probabilidades
#   out$qx_out_of_range   <- x[!is.na(qx) & (qx < 0 | qx > 1)]
#   out$px_out_of_range   <- x[!is.na(px_5y) & (px_5y < 0 | px_5y > 1)]
#   
#   # 2) No negativos (tasas y conteos)
#   nonneg_cols <- c("mx","lx","dx","Lx","Tx","ex")
#   for (cc in nonneg_cols) {
#     if (cc %in% names(x)) out[[paste0(cc, "_negative")]] <- x[get(cc) < 0]
#   }
#   
#   # 3) Chequeo simple: period_end_year > period_start_year
#   out$bad_period_order <- x[period_end_year <= period_start_year]
#   
#   # 4) age_interval_width positivo
#   out$bad_age_width <- x[age_interval_width <= 0]
#   
#   out
# }
# 
# write_qc_list <- function(qc_list, qc_dir) {
#   dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
#   for (nm in names(qc_list)) {
#     fp <- file.path(qc_dir, paste0("qc_", nm, ".csv"))
#     fwrite(qc_list[[nm]], fp)
#   }
#   invisible(TRUE)
# }

# R/qc_life_table_utils.R
library(data.table)

qc_life_table_basic <- function(dt) {
  x <- as.data.table(dt)
  
  out <- list()
  
  # 1) Rango probabilidades
  if ("qx" %in% names(x)) {
    out$qx_out_of_range <- x[!is.na(qx) & (qx < 0 | qx > 1)]
  }
  
  if ("px_5y" %in% names(x)) {
    out$px_out_of_range <- x[!is.na(px_5y) & (px_5y < 0 | px_5y > 1)]
  }
  
  # 2) No negativos
  nonneg_cols <- c("mx", "lx", "dx", "Lx", "Tx", "ex")
  for (cc in nonneg_cols) {
    if (cc %in% names(x)) {
      out[[paste0(cc, "_negative")]] <- x[get(cc) < 0]
    }
  }
  
  # 3) Orden temporal solo si existen columnas de periodo
  if (all(c("period_start_year", "period_end_year") %in% names(x))) {
    out$bad_period_order <- x[period_end_year <= period_start_year]
  }
  
  # 4) age_interval_width positivo, pero permitir NA si es intervalo abierto
  if ("age_interval_width" %in% names(x)) {
    if ("age_interval_open" %in% names(x)) {
      out$bad_age_width <- x[
        !is.na(age_interval_width) & age_interval_width <= 0 |
          (is.na(age_interval_width) & !age_interval_open)
      ]
    } else {
      out$bad_age_width <- x[!is.na(age_interval_width) & age_interval_width <= 0]
    }
  }
  
  out
}

write_qc_list <- function(qc_list, qc_dir) {
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (nm in names(qc_list)) {
    fp <- file.path(qc_dir, paste0("qc_", nm, ".csv"))
    fwrite(qc_list[[nm]], fp)
  }
  
  invisible(TRUE)
}