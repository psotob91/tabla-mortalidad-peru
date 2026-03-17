# R/qc_utils.R
library(data.table)

qc_pk_duplicates <- function(dt, pk) {
  dt[, .N, by = pk][N > 1]
}

qc_missing_required <- function(dt, req_cols) {
  data.table(
    column = req_cols,
    n_missing = sapply(req_cols, \(c) sum(is.na(dt[[c]]))),
    pct_missing = sapply(req_cols, \(c) round(100 * mean(is.na(dt[[c]])), 6))
  )
}

qc_nonnegative <- function(dt, col) {
  dt[get(col) < 0]
}

# FIX: setorder(x, by_cols, age_col) -> setorderv con vector de columnas
qc_tail_monotone_flag <- function(dt,
                                  age_col = "age",
                                  pop_col = "population",
                                  start_age = 70L,
                                  by_cols = c("year_id", "sex_id", "location_id")) {
  x <- as.data.table(dt)[get(age_col) >= start_age]
  if (nrow(x) == 0) return(data.table())
  
  # ordenar por columnas dinámicas
  data.table::setorderv(x, cols = c(by_cols, age_col))
  
  x[, pop_prev := shift(get(pop_col)), by = by_cols]
  x[, inc_flag := !is.na(pop_prev) & get(pop_col) > pop_prev]
  
  x[inc_flag == TRUE, .N, by = by_cols]
}

# QC detectivo: national (00) vs suma deptos (01–25)
qc_national_vs_dept_sum <- function(dt,
                                    location_col = "location_id",
                                    national_id = 0L,
                                    dept_ids = 1:25,
                                    by_cols = c("year_id", "sex_id", "age"),
                                    pop_col = "population",
                                    pct_tol = 0.001) {
  x <- as.data.table(dt)
  
  nat <- x[get(location_col) == national_id,
           .(pop_national = sum(get(pop_col), na.rm = TRUE)),
           by = by_cols]
  
  dep <- x[get(location_col) %in% dept_ids,
           .(pop_dept_sum = sum(get(pop_col), na.rm = TRUE)),
           by = by_cols]
  
  out <- merge(nat, dep, by = by_cols, all = TRUE)
  
  out[, diff_abs := pop_national - pop_dept_sum]
  out[, diff_pct := fifelse(!is.na(pop_dept_sum) & pop_dept_sum > 0,
                            100 * diff_abs / pop_dept_sum,
                            NA_real_)]
  
  out[, flag_diff := fifelse(!is.na(diff_pct) & abs(diff_pct) > (pct_tol * 100), TRUE, FALSE)]
  out[, flag_missing := is.na(pop_national) | is.na(pop_dept_sum)]
  
  data.table::setorderv(out, cols = by_cols)
  out[]
}

qc_hierarchical_national_additive_hard <- function(dt,
                                                   base_locations = 1:25,
                                                   national_additive_id = 9000L) {
  stopifnot(all(c("year_id","age","sex_id","location_id","population") %in% names(dt)))
  
  # 1) No permitir location_id = 0
  if (any(dt$location_id == 0L)) {
    stop("QC HARD FAIL: Vista jerárquica contiene location_id=0 (nacional original). Debe excluirse.")
  }
  
  # 2) Debe existir 9000
  if (!any(dt$location_id == national_additive_id)) {
    stop("QC HARD FAIL: Vista jerárquica no contiene location_id=9000 (nacional aditivo).")
  }
  
  # 3) Debe cubrir exactamente base_locations
  missing_base <- setdiff(base_locations, unique(dt$location_id))
  if (length(missing_base) > 0) {
    stop("QC HARD FAIL: faltan deptos base en vista jerárquica: ", paste(missing_base, collapse = ", "))
  }
  
  # 4) Chequeo exacto: 9000 == suma deptos por celda
  dt_dept <- dt[location_id %in% base_locations,
                .(dept_sum = sum(population)), by = .(year_id, age, sex_id)]
  
  dt_nat <- dt[location_id == national_additive_id,
               .(nat_pop = population), by = .(year_id, age, sex_id)]
  
  chk <- merge(dt_dept, dt_nat, by = c("year_id","age","sex_id"), all = TRUE)
  
  # NAs son errores de cobertura
  if (anyNA(chk$dept_sum) || anyNA(chk$nat_pop)) {
    stop("QC HARD FAIL: hay celdas donde falta dept_sum o nat_pop (cobertura incompleta).")
  }
  
  chk[, diff := nat_pop - dept_sum]
  bad <- chk[diff != 0]
  
  if (nrow(bad) > 0) {
    # no guardamos silenciosamente; fallamos mostrando un resumen
    top <- bad[order(-abs(diff))][1:min(20, .N)]
    stop(
      "QC HARD FAIL: nacional aditivo (9000) NO es suma exacta de deptos.\n",
      "Ejemplos (top 20 por |diff|):\n",
      paste(capture.output(print(top)), collapse = "\n")
    )
  }
  
  invisible(TRUE)
}
