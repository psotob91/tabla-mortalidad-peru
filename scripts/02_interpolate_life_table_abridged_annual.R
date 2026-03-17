# scripts/02_interpolate_life_table_abridged_annual.R
# ------------------------------------------------------------
# Script 02
# - Toma la tabla de vida quinquenal/por periodos del Script 01
# - Interpola años simples 2018:2024 usando period_start_year como ancla
# - Interpola qx en escala cloglog (recomendado)
# - Reconstruye tabla de vida ABRIDGED anual consistente
# - Respeta estructura del proyecto basada en data/, no datos/
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(dplyr)
})

# ------------------------------------------------------------
# Sources del proyecto
# ------------------------------------------------------------
source(here("R/io_utils.R"))
source(here("R/spec_utils.R"))
source(here("R/dictionary_utils.R"))
source(here("R/qc_life_table_utils.R"))
source(here("config/parametros-proyecto.R"))

# ------------------------------------------------------------
# Config basada en estructura real del proyecto
# ------------------------------------------------------------
P <- paths_mort()

CFG <- list(
  in_file   = file.path(P$FINAL_DIR, "ref_life_table_mortality.csv"),
  out_dir   = file.path(P$FINAL_DIR, "annual_abridged"),
  qc_dir    = file.path(P$QC_DIR, "annual_abridged"),
  out_file  = "ref_life_table_mortality_annual_abridged.csv",
  out_dict  = "ref_life_table_mortality_annual_abridged_dictionary_ext.csv",
  spec_file = file.path(P$CONFIG_DIR, "spec_life_table_mortality_annual_abridged.yml"),
  interp_scale = "cloglog",
  radix = 100000,
  eps = 1e-12
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CFG$qc_dir, recursive = TRUE, showWarnings = FALSE)

# Blindaje contra rutas mal puestas
bad_paths <- c(CFG$in_file, CFG$out_dir, CFG$qc_dir, CFG$spec_file)
if (any(grepl("(^|[/\\\\])datos([/\\\\]|$)", bad_paths))) {
  stop("Se detectó una ruta con 'datos/'. Este proyecto usa 'data/'.")
}

if (!file.exists(CFG$in_file)) {
  stop("No existe el archivo de entrada: ", CFG$in_file)
}

if (!file.exists(CFG$spec_file)) {
  stop("No existe el spec requerido: ", CFG$spec_file)
}

# Años objetivo
if (exists("TARGET_YEARS_LIFE_TABLE")) {
  target_years <- TARGET_YEARS_LIFE_TABLE
} else if (exists("YEARS_ANALISIS")) {
  target_years <- YEARS_ANALISIS
} else {
  target_years <- 2018:2024
}

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
safe_num <- function(x) suppressWarnings(as.numeric(x))

cloglog <- function(q, eps = 1e-12) {
  q <- pmin(1 - eps, pmax(eps, q))
  log(-log(1 - q))
}

inv_cloglog <- function(z) {
  1 - exp(-exp(z))
}

find_bracketing_years <- function(y, anchors) {
  anchors <- sort(unique(anchors))
  
  if (y < min(anchors) || y > max(anchors)) {
    stop(
      "Año ", y, " fuera del rango cubierto por anclas: ",
      paste(anchors, collapse = ", ")
    )
  }
  
  if (y %in% anchors) return(c(y, y))
  
  left  <- max(anchors[anchors < y])
  right <- min(anchors[anchors > y])
  
  c(left, right)
}

# Reconstrucción de tabla abreviada anual a partir de qx anual interpolado
reconstruct_life_table_abridged <- function(dt_group, radix = 100000, eps = 1e-12) {
  g <- as.data.table(copy(dt_group))
  setorder(g, age_start)
  
  # asegurar qx válido
  g[, qx := pmin(1 - eps, pmax(eps, qx))]
  
  # prob. de sobrevivir del intervalo
  g[, px_5y := 1 - qx]
  
  # lx
  g[, lx := NA_real_]
  g[1, lx := radix]
  
  if (nrow(g) > 1) {
    for (i in 2:nrow(g)) {
      g[i, lx := g[i - 1, lx] * (1 - g[i - 1, qx])]
    }
  }
  
  # dx
  g[, dx := lx * qx]
  
  # hazard medio del intervalo cerrado
  g[, mx := fifelse(
    !is.na(age_interval_width) & age_interval_width > 0,
    -log(1 - qx) / age_interval_width,
    NA_real_
  )]
  
  # Lx cerrado bajo hazard constante
  g[, Lx := fifelse(
    !is.na(age_interval_width) & age_interval_width > 0 & mx > 0,
    dx / mx,
    NA_real_
  )]
  
  # respaldo cuando mx ~ 0
  g[is.na(Lx) & !is.na(age_interval_width) & age_interval_width > 0,
    Lx := lx * age_interval_width]
  
  # intervalo abierto
  if ("age_interval_open" %in% names(g) && any(g$age_interval_open, na.rm = TRUE)) {
    idx_open <- which(g$age_interval_open %in% TRUE)
    
    # usar hazard implícito positivo derivado de qx anualizada del intervalo abierto
    g[idx_open, mx := -log(1 - qx)]
    g[idx_open, Lx := lx / pmax(mx, eps)]
  }
  
  # Tx y ex
  g[, Tx := rev(cumsum(rev(Lx)))]
  g[, ex := Tx / pmax(lx, eps)]
  
  g[]
}

# ------------------------------------------------------------
# 1) Leer input del script 01
# ------------------------------------------------------------
cat("==> Leyendo input abridged: ", CFG$in_file, "\n", sep = "")
dt <- fread(CFG$in_file)

# Tipos esperados
if ("age_start" %in% names(dt)) dt[, age_start := as.integer(age_start)]
if ("age_interval_width" %in% names(dt)) dt[, age_interval_width := as.integer(age_interval_width)]
if ("period_start_year" %in% names(dt)) dt[, period_start_year := as.integer(period_start_year)]
if ("period_end_year" %in% names(dt)) dt[, period_end_year := as.integer(period_end_year)]
if ("sex_id" %in% names(dt)) dt[, sex_id := as.integer(sex_id)]
if ("location_id" %in% names(dt)) dt[, location_id := as.integer(location_id)]
if ("qx" %in% names(dt)) dt[, qx := safe_num(qx)]

# age_interval_open puede no existir en algunas corridas antiguas
if (!"age_interval_open" %in% names(dt)) {
  dt[, age_interval_open := FALSE]
  
  dt[
    ,
    age_interval_open := age_start == max(age_start, na.rm = TRUE),
    by = .(period_start_year, period_end_year, location_id, sex_id, life_table_label)
  ]
}

dt[, age_interval_open := as.logical(age_interval_open)]

needed <- c(
  "period_start_year", "period_end_year",
  "location_id", "location_name",
  "sex_id", "sex_source_value",
  "life_table_label",
  "age_start", "age_interval_width", "age_interval_open",
  "qx"
)

miss <- setdiff(needed, names(dt))
if (length(miss) > 0) {
  stop("Faltan columnas requeridas en input: ", paste(miss, collapse = ", "))
}

# ------------------------------------------------------------
# 2) Validación básica del input
# ------------------------------------------------------------
anchors <- sort(unique(dt$period_start_year))

cat("==> Años ancla detectados: ", paste(anchors, collapse = ", "), "\n", sep = "")
cat("==> Años objetivo: ", paste(target_years, collapse = ", "), "\n", sep = "")

if (min(target_years) < min(anchors) || max(target_years) > max(anchors)) {
  stop(
    "Los años objetivo salen fuera del rango de anclas. ",
    "Anclas: ", paste(anchors, collapse = ", "),
    ". Objetivo: ", paste(target_years, collapse = ", ")
  )
}

# unicidad por celda y ancla
anchor_pk <- c(
  "period_start_year",
  "location_id", "location_name",
  "sex_id", "sex_source_value",
  "life_table_label",
  "age_start", "age_interval_width", "age_interval_open"
)

dup_chk <- dt[, .N, by = anchor_pk][N > 1]
fwrite(dup_chk, file.path(CFG$qc_dir, "qc_duplicates_by_anchor.csv"))
if (nrow(dup_chk) > 0) {
  stop("Hay duplicados por celda etaria y año ancla. Revisa qc_duplicates_by_anchor.csv")
}

# ------------------------------------------------------------
# 3) Interpolación anual de qx por celda etaria
# ------------------------------------------------------------
cat("==> Interpolando qx anual en escala ", CFG$interp_scale, "...\n", sep = "")

interp_list <- lapply(target_years, function(yy) {
  
  br <- find_bracketing_years(yy, anchors)
  y0 <- br[1]
  y1 <- br[2]
  
  # caso exacto en ancla
  if (y0 == y1) {
    ans <- dt[period_start_year == y0][
      ,
      .(
        year_id = as.integer(yy),
        source_year_left = as.integer(y0),
        source_year_right = as.integer(y1),
        interp_weight = 0,
        interp_scale = CFG$interp_scale,
        
        location_id, location_name,
        sex_id, sex_source_value,
        life_table_label,
        age_start, age_interval_width, age_interval_open,
        qx
      )
    ]
    return(ans)
  }
  
  w <- (yy - y0) / (y1 - y0)
  
  left_dt <- dt[period_start_year == y0][
    ,
    .(
      location_id, location_name,
      sex_id, sex_source_value,
      life_table_label,
      age_start, age_interval_width, age_interval_open,
      qx_left = qx
    )
  ]
  
  right_dt <- dt[period_start_year == y1][
    ,
    .(
      location_id, location_name,
      sex_id, sex_source_value,
      life_table_label,
      age_start, age_interval_width, age_interval_open,
      qx_right = qx
    )
  ]
  
  m <- merge(
    left_dt,
    right_dt,
    by = c(
      "location_id", "location_name",
      "sex_id", "sex_source_value",
      "life_table_label",
      "age_start", "age_interval_width", "age_interval_open"
    ),
    all = FALSE
  )
  
  if (CFG$interp_scale == "cloglog") {
    m[, qx := inv_cloglog(cloglog(qx_left, eps = CFG$eps) +
                            w * (cloglog(qx_right, eps = CFG$eps) - cloglog(qx_left, eps = CFG$eps)))]
  } else if (CFG$interp_scale == "linear") {
    m[, qx := qx_left + w * (qx_right - qx_left)]
  } else {
    stop("Escala de interpolación no soportada: ", CFG$interp_scale)
  }
  
  m[, `:=`(
    year_id = as.integer(yy),
    source_year_left = as.integer(y0),
    source_year_right = as.integer(y1),
    interp_weight = as.numeric(w),
    interp_scale = CFG$interp_scale
  )]
  
  m[, .(
    year_id,
    source_year_left,
    source_year_right,
    interp_weight,
    interp_scale,
    
    location_id, location_name,
    sex_id, sex_source_value,
    life_table_label,
    age_start, age_interval_width, age_interval_open,
    qx
  )]
})

dt_interp_qx <- rbindlist(interp_list, use.names = TRUE)

# ------------------------------------------------------------
# 4) Reconstrucción de tabla abridged anual
# ------------------------------------------------------------
cat("==> Reconstruyendo tabla de vida abridged anual...\n")

lt_group_vars <- c(
  "year_id",
  "source_year_left", "source_year_right", "interp_weight", "interp_scale",
  "location_id", "location_name",
  "sex_id", "sex_source_value",
  "life_table_label"
)

setorderv(
  dt_interp_qx,
  c("year_id", "location_id", "sex_id", "life_table_label", "age_start")
)

dt_annual <- dt_interp_qx[
  ,
  reconstruct_life_table_abridged(.SD, radix = CFG$radix, eps = CFG$eps),
  by = lt_group_vars
]

setcolorder(
  dt_annual,
  c(
    "year_id",
    "source_year_left", "source_year_right", "interp_weight", "interp_scale",
    "location_id", "location_name",
    "sex_id", "sex_source_value",
    "life_table_label",
    "age_start", "age_interval_width", "age_interval_open",
    "mx", "qx", "px_5y", "lx", "dx", "Lx", "Tx", "ex"
  )
)

# ------------------------------------------------------------
# 5) QC
# ------------------------------------------------------------
cat("==> Ejecutando QC...\n")

# PK lógica
pk_chk <- dt_annual[
  ,
  .N,
  by = .(year_id, age_start, sex_id, location_id, life_table_label)
][N > 1]
fwrite(pk_chk, file.path(CFG$qc_dir, "qc_pk_duplicates.csv"))

# qx
qc_bad_qx <- dt_annual[qx < 0 | qx > 1 | is.na(qx)]
fwrite(qc_bad_qx, file.path(CFG$qc_dir, "qc_bad_qx.csv"))

# intervalos abiertos inconsistentes
qc_bad_open_flag <- dt_annual[is.na(age_interval_width) & age_interval_open == FALSE]
fwrite(qc_bad_open_flag, file.path(CFG$qc_dir, "qc_bad_open_flag.csv"))

# más de un intervalo abierto por estrato
qc_open_count <- dt_annual[
  ,
  .(n_open = sum(age_interval_open %in% TRUE, na.rm = TRUE)),
  by = .(year_id, location_id, sex_id, life_table_label)
]
fwrite(qc_open_count, file.path(CFG$qc_dir, "qc_open_interval_count.csv"))

# monotonicidad lx
qc_lx <- dt_annual[
  order(year_id, location_id, sex_id, life_table_label, age_start)
][
  ,
  .(
    n_bad_lx = sum(lx > shift(lx, type = "lag", fill = Inf), na.rm = TRUE)
  ),
  by = .(year_id, location_id, sex_id, life_table_label)
]
fwrite(qc_lx, file.path(CFG$qc_dir, "qc_lx_monotonicity_summary.csv"))

# monotonicidad Tx
qc_tx <- dt_annual[
  order(year_id, location_id, sex_id, life_table_label, age_start)
][
  ,
  .(
    n_bad_tx = sum(Tx > shift(Tx, type = "lag", fill = Inf), na.rm = TRUE)
  ),
  by = .(year_id, location_id, sex_id, life_table_label)
]
fwrite(qc_tx, file.path(CFG$qc_dir, "qc_tx_monotonicity_summary.csv"))

# qc básico reutilizando helper del proyecto
qc_basic <- qc_life_table_basic(as.data.table(dt_annual))
write_qc_list(qc_basic, file.path(CFG$qc_dir, "qc_basic"))

# fallas duras mínimas antes del spec
if (nrow(pk_chk) > 0) stop("FAIL QC: duplicados en PK lógica. Revisa qc_pk_duplicates.csv")
if (nrow(qc_bad_qx) > 0) stop("FAIL QC: qx fuera de rango o NA. Revisa qc_bad_qx.csv")
if (nrow(qc_bad_open_flag) > 0) stop("FAIL QC: intervalos abiertos mal marcados. Revisa qc_bad_open_flag.csv")
if (any(qc_open_count$n_open != 1)) {
  fwrite(qc_open_count[n_open != 1], file.path(CFG$qc_dir, "qc_bad_open_interval_count.csv"))
  stop("FAIL QC: cada estrato debe tener exactamente un intervalo abierto. Revisa qc_bad_open_interval_count.csv")
}

# ------------------------------------------------------------
# 6) Validación por spec
# ------------------------------------------------------------
cat("==> Validando contra spec annual abridged...\n")

spec <- read_spec(CFG$spec_file)
validate_by_spec(as.data.table(dt_annual), spec)

# ------------------------------------------------------------
# 7) Diccionario extendido
# ------------------------------------------------------------
cat("==> Construyendo diccionario extendido...\n")

dict <- dict_from_spec(
  spec = spec,
  dataset_version = "v1.0.0",
  run_id = NA_character_,
  config_dir = P$CONFIG_DIR
)

dict_ext <- enrich_dict_with_stats(dict, as.data.table(dt_annual))

# ------------------------------------------------------------
# 8) Export
# ------------------------------------------------------------
out_csv  <- file.path(CFG$out_dir, CFG$out_file)
out_dict <- file.path(CFG$out_dir, CFG$out_dict)

fwrite(as.data.table(dt_annual), out_csv)
fwrite(as.data.table(dict_ext), out_dict)

cat("==> OK export dataset: ", out_csv, "\n", sep = "")
cat("==> OK export dictionary: ", out_dict, "\n", sep = "")
cat("==> OK QC dir: ", CFG$qc_dir, "\n", sep = "")
cat("==> Script 02 finalizado correctamente.\n")