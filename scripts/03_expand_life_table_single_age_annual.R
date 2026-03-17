# scripts/03_expand_life_table_single_age_annual.R
# ------------------------------------------------------------
# Script 03
# - Toma la tabla de vida annual abridged del Script 02
# - Expande a edad simple con un esquema demográfico por tramos:
#     * edad 0: fija
#     * 1-14: spline suave sobre log(mx) con ancla de continuidad en 15
#     * 15-84: Gompertz sobre log(mx)
#     * 85-109: Kannisto-like sobre qx (logit)
#     * 110+: intervalo abierto real
# - Recalibra qx SOLO en intervalos abridged cerrados
# - Reconstruye tabla de vida single-age consistente
# - Valida conservación de qx SOLO para intervalos cerrados
# - Añade QC específico para el intervalo abierto 110+
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
# Config
# ------------------------------------------------------------
P <- paths_mort()

CFG <- list(
  in_file          = file.path(P$FINAL_DIR, "annual_abridged", "ref_life_table_mortality_annual_abridged.csv"),
  out_dir          = file.path(P$FINAL_DIR, "single_age"),
  qc_dir           = file.path(P$QC_DIR, "single_age"),
  out_file         = "ref_life_table_mortality_single_age.csv",
  out_dict         = "ref_life_table_mortality_single_age_dictionary_ext.csv",
  spec_file        = file.path(P$CONFIG_DIR, "spec_life_table_mortality_single_age.yml"),
  age_open_start   = 110L,   # 110+
  radix            = 100000,
  eps_q            = 1e-12,
  qx_gold_tol      = 1e-10,
  old_age_cut      = 85L,
  adult_age_cut    = 15L
)

CFG$age_last_closed <- CFG$age_open_start - 1L

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CFG$qc_dir, recursive = TRUE, showWarnings = FALSE)

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

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
safe_num <- function(x) suppressWarnings(as.numeric(x))

qc_print <- function(dt, n = 20) {
  print(dt[1:min(.N, n)])
  invisible(dt)
}

clip_qx <- function(q, eps = 1e-12) {
  pmin(1 - eps, pmax(0, safe_num(q)))
}

qx_to_mu <- function(q, width = 1, eps = 1e-12) {
  q <- clip_qx(q, eps = eps)
  width <- pmax(1, safe_num(width))
  -log(1 - pmin(1 - eps, q)) / width
}

mu_to_qx <- function(mu, width = 1) {
  mu <- pmax(0, safe_num(mu))
  width <- pmax(1, safe_num(width))
  1 - exp(-mu * width)
}

interval_midpoint <- function(age_start, width, age_interval_open = FALSE, fallback_open_age = NULL) {
  out <- age_start + width / 2
  if (any(age_interval_open %in% TRUE, na.rm = TRUE)) {
    idx <- which(age_interval_open %in% TRUE)
    if (length(idx) > 0) {
      if (is.null(fallback_open_age)) {
        out[idx] <- age_start[idx]
      } else {
        out[idx] <- (age_start[idx] + fallback_open_age) / 2
      }
    }
  }
  out
}

# ------------------------------------------------------------
# Ajustadores por tramo
# ------------------------------------------------------------

# 1-14: spline suave sobre log(mx) con ancla en 15 para continuidad
predict_child_mu <- function(
    g,
    ages,
    mu15_boundary = NULL,
    eps_q = 1e-12
) {
  if (length(ages) == 0) {
    return(data.table(age_start = integer(), mu_base = numeric()))
  }
  
  child_ref <- copy(g[
    age_start >= 1 & age_start <= 10 & age_interval_open == FALSE,
    .(age_start, age_interval_width, age_interval_open, qx, mx)
  ])
  
  if (nrow(child_ref) == 0) {
    return(data.table(age_start = ages, mu_base = rep(1e-8, length(ages))))
  }
  
  child_ref[, mx_ref := fifelse(!is.na(mx) & mx > 0, mx, qx_to_mu(qx, age_interval_width, eps_q))]
  child_ref[, age_mid := interval_midpoint(age_start, age_interval_width, age_interval_open)]
  child_ref <- child_ref[is.finite(age_mid) & is.finite(mx_ref) & mx_ref > 0]
  
  anchor_x <- child_ref$age_mid
  anchor_y <- log(child_ref$mx_ref)
  
  if (!is.null(mu15_boundary) && is.finite(mu15_boundary) && mu15_boundary > 0) {
    anchor_x <- c(anchor_x, 15)
    anchor_y <- c(anchor_y, log(mu15_boundary))
  }
  
  ord <- order(anchor_x)
  anchor_x <- anchor_x[ord]
  anchor_y <- anchor_y[ord]
  
  if (length(anchor_x) >= 4) {
    fit <- smooth.spline(x = anchor_x, y = anchor_y, spar = 0.65)
    mu_hat <- exp(predict(fit, x = ages)$y)
  } else if (length(anchor_x) >= 2) {
    mu_hat <- exp(approx(anchor_x, anchor_y, xout = ages, rule = 2)$y)
  } else {
    mu_hat <- rep(exp(anchor_y[1]), length(ages))
  }
  
  mu_hat <- pmax(mu_hat, 1e-12)
  
  data.table(age_start = ages, mu_base = mu_hat)
}

# 15-84: Gompertz via lm(log(mx) ~ age)
predict_adult_mu <- function(
    g,
    ages,
    eps_q = 1e-12,
    adult_age_cut = 15L,
    old_age_cut = 85L
) {
  if (length(ages) == 0) {
    return(data.table(age_start = integer(), mu_base = numeric()))
  }
  
  adult_ref <- copy(g[
    age_start >= adult_age_cut & age_start < old_age_cut & age_interval_open == FALSE,
    .(age_start, age_interval_width, age_interval_open, qx, mx)
  ])
  
  if (nrow(adult_ref) == 0) {
    return(data.table(age_start = ages, mu_base = rep(1e-8, length(ages))))
  }
  
  adult_ref[, mx_ref := fifelse(!is.na(mx) & mx > 0, mx, qx_to_mu(qx, age_interval_width, eps_q))]
  adult_ref[, age_mid := interval_midpoint(age_start, age_interval_width, age_interval_open)]
  adult_ref <- adult_ref[is.finite(age_mid) & is.finite(mx_ref) & mx_ref > 0]
  
  if (nrow(adult_ref) >= 2) {
    fit <- lm(log(mx_ref) ~ age_mid, data = adult_ref)
    mu_hat <- exp(predict(fit, newdata = data.frame(age_mid = ages)))
  } else {
    mu_hat <- rep(adult_ref$mx_ref[1], length(ages))
  }
  
  data.table(age_start = ages, mu_base = pmax(mu_hat, 1e-12))
}

# 85+: Kannisto-like sobre qx (logit lineal en edad)
predict_old_qx_kannisto_like <- function(
    g,
    ages,
    eps_q = 1e-12,
    old_age_cut = 85L
) {
  if (length(ages) == 0) {
    return(data.table(age_start = integer(), qx_base = numeric()))
  }
  
  old_ref <- copy(g[
    age_start >= max(80L, old_age_cut - 5L) & age_interval_open == FALSE,
    .(age_start, age_interval_width, age_interval_open, qx)
  ])
  
  old_ref[, age_mid := interval_midpoint(age_start, age_interval_width, age_interval_open)]
  old_ref[, qx_ref := clip_qx(qx, eps = 1e-8)]
  old_ref <- old_ref[is.finite(age_mid) & is.finite(qx_ref) & qx_ref > 0 & qx_ref < 1]
  
  if (nrow(old_ref) >= 2) {
    fit <- lm(qlogis(qx_ref) ~ age_mid, data = old_ref)
    qx_hat <- plogis(predict(fit, newdata = data.frame(age_mid = ages)))
  } else if (nrow(old_ref) == 1) {
    qx_hat <- rep(old_ref$qx_ref[1], length(ages))
  } else {
    qx_hat <- rep(0.5, length(ages))
  }
  
  qx_hat <- pmin(0.995, pmax(1e-8, qx_hat))
  
  data.table(age_start = ages, qx_base = qx_hat)
}

# ------------------------------------------------------------
# Construcción de forma base single-age
# ------------------------------------------------------------
build_single_age_shape <- function(
    g,
    age_last_closed = 109L,
    age_open_start = 110L,
    eps_q = 1e-12,
    adult_age_cut = 15L,
    old_age_cut = 85L
) {
  stopifnot(is.data.table(g))
  g <- copy(g)
  setorder(g, age_start)
  
  ages_closed <- 0:age_last_closed
  
  base <- data.table(
    age_start = ages_closed,
    qx_base = NA_real_,
    mu_base = NA_real_
  )
  
  # edad 0 fija
  row0 <- g[age_start == 0]
  if (nrow(row0) == 0) {
    stop("No se encontró el intervalo de edad 0 en la tabla abridged del grupo.")
  }
  
  q0 <- clip_qx(row0$qx[1], eps = eps_q)
  base[age_start == 0, `:=`(
    qx_base = q0,
    mu_base = qx_to_mu(q0, width = 1, eps = eps_q)
  )]
  
  # tramo adulto primero para obtener frontera en 15
  ages_adult <- ages_closed[ages_closed >= adult_age_cut & ages_closed < old_age_cut]
  adult_dt <- predict_adult_mu(
    g = g,
    ages = ages_adult,
    eps_q = eps_q,
    adult_age_cut = adult_age_cut,
    old_age_cut = old_age_cut
  )
  
  mu15_boundary <- adult_dt[age_start == adult_age_cut, mu_base][1]
  if (is.na(mu15_boundary) || !is.finite(mu15_boundary)) mu15_boundary <- 1e-6
  
  # infancia 1-14 con empalme al adulto
  ages_child <- ages_closed[ages_closed >= 1 & ages_closed < adult_age_cut]
  child_dt <- predict_child_mu(
    g = g,
    ages = ages_child,
    mu15_boundary = mu15_boundary,
    eps_q = eps_q
  )
  
  if (nrow(child_dt) > 0) {
    base[child_dt, on = "age_start", mu_base := i.mu_base]
    base[child_dt, on = "age_start", qx_base := mu_to_qx(i.mu_base, width = 1)]
  }
  
  if (nrow(adult_dt) > 0) {
    base[adult_dt, on = "age_start", mu_base := i.mu_base]
    base[adult_dt, on = "age_start", qx_base := mu_to_qx(i.mu_base, width = 1)]
  }
  
  # cola 85-109 con Kannisto-like
  ages_old_closed <- ages_closed[ages_closed >= old_age_cut]
  old_dt_closed <- predict_old_qx_kannisto_like(
    g = g,
    ages = ages_old_closed,
    eps_q = eps_q,
    old_age_cut = old_age_cut
  )
  
  if (nrow(old_dt_closed) > 0) {
    base[old_dt_closed, on = "age_start", qx_base := i.qx_base]
    base[old_dt_closed, on = "age_start", mu_base := qx_to_mu(i.qx_base, width = 1, eps = eps_q)]
  }
  
  # Relleno final si quedara algún hueco
  if (any(is.na(base$mu_base))) {
    miss_ages <- base[is.na(mu_base), age_start]
    obs <- base[!is.na(mu_base)]
    
    if (nrow(obs) >= 2) {
      mu_fill <- approx(
        x = obs$age_start,
        y = log(obs$mu_base),
        xout = miss_ages,
        rule = 2
      )$y
      base[age_start %in% miss_ages, mu_base := exp(mu_fill)]
      base[age_start %in% miss_ages, qx_base := mu_to_qx(mu_base, width = 1)]
    } else {
      base[is.na(mu_base), mu_base := 1e-8]
      base[is.na(qx_base), qx_base := mu_to_qx(mu_base, width = 1)]
    }
  }
  
  # -------- NUEVO: extrapolar explícitamente 110 con Kannisto-like --------
  old_dt_open <- predict_old_qx_kannisto_like(
    g = g,
    ages = age_open_start,
    eps_q = eps_q,
    old_age_cut = old_age_cut
  )
  
  qx_open <- old_dt_open$qx_base[1]
  if (is.na(qx_open) || !is.finite(qx_open)) {
    # fallback conservador, pero ya no copia literalmente 109 si hay tendencia estimable
    qx109 <- tail(base$qx_base, 1)
    qx108 <- tail(base$qx_base, 2)[1]
    growth <- qx109 / pmax(qx108, 1e-12)
    qx_open <- pmin(0.995, qx109 * growth)
  }
  
  mu_open <- qx_to_mu(qx_open, width = 1, eps = eps_q)
  
  base[, mu_base := pmax(mu_base, 1e-12)]
  base[, qx_base := clip_qx(qx_base, eps = eps_q)]
  
  list(
    base_closed = base[],
    age_open_start = age_open_start,
    age_last_closed = age_last_closed,
    qx_open = qx_open,
    mu_open = mu_open
  )
}

# ------------------------------------------------------------
# Recalibración por intervalos abridged cerrados
# ------------------------------------------------------------
calibrate_single_age_to_abridged_qx <- function(
    g,
    base_shape_obj,
    eps_q = 1e-12
) {
  g <- copy(g)
  setorder(g, age_start)
  
  one <- copy(base_shape_obj$base_closed)
  one[, `:=`(
    qx = NA_real_,
    mx = NA_real_
  )]
  
  g[, age_end := age_start + age_interval_width]
  
  # calibrar solo intervalos cerrados
  g_closed <- g[age_interval_open == FALSE]
  g_open   <- g[age_interval_open == TRUE]
  
  for (i in seq_len(nrow(g_closed))) {
    x0 <- g_closed$age_start[i]
    x1 <- g_closed$age_end[i]
    idx <- one$age_start >= x0 & one$age_start < x1
    
    if (!any(idx)) {
      stop("No se encontraron edades simples para el intervalo abridged [", x0, ", ", x1, ").")
    }
    
    nqx_target <- clip_qx(g_closed$qx[i], eps = eps_q)
    H_target <- -log(1 - pmin(1 - eps_q, nqx_target))
    
    mu_base <- pmax(one$mu_base[idx], 1e-12)
    H_base <- sum(mu_base)
    
    if (!is.finite(H_base) || H_base <= 0) {
      mu_cal <- rep(H_target / sum(idx), sum(idx))
    } else {
      scale_fac <- H_target / H_base
      mu_cal <- mu_base * scale_fac
    }
    
    qx_cal <- mu_to_qx(mu_cal, width = 1)
    
    one[idx, `:=`(
      mx = mu_cal,
      qx = qx_cal
    )]
  }
  
  # si existe intervalo abierto abridged, mantener forma modelada desde su inicio
  if (nrow(g_open) > 0) {
    x0_open <- g_open$age_start[1]
    idx_open_tail <- one$age_start >= x0_open
    
    if (any(idx_open_tail)) {
      one[idx_open_tail, mx := pmax(mu_base, 1e-12)]
      one[idx_open_tail, qx := mu_to_qx(mx, width = 1)]
    }
  }
  
  # agregar fila final 110+ como intervalo abierto real
  open_row <- data.table(
    age_start = base_shape_obj$age_open_start,
    qx = clip_qx(base_shape_obj$qx_open, eps = eps_q),
    mx = pmax(base_shape_obj$mu_open, 1e-12)
  )
  
  out <- rbind(
    one[, .(age_start, qx, mx)],
    open_row,
    fill = TRUE
  )
  
  out[, age_interval_width := 1L]
  out[, age_interval_open := FALSE]
  out[age_start == base_shape_obj$age_open_start, age_interval_open := TRUE]
  
  setorder(out, age_start)
  out[]
}

# ------------------------------------------------------------
# Reconstrucción de tabla de vida single-age
# ------------------------------------------------------------
reconstruct_single_age_life_table <- function(one, radix = 100000, eps_q = 1e-12) {
  one <- copy(one)
  setorder(one, age_start)
  
  one[, qx := clip_qx(qx, eps = eps_q)]
  one[, mx := pmax(safe_num(mx), 0)]
  
  one[, lx := NA_real_]
  one[1, lx := radix]
  
  if (nrow(one) > 1) {
    for (i in 2:nrow(one)) {
      one[i, lx := one[i - 1, lx] * (1 - one[i - 1, qx])]
    }
  }
  
  one[, dx := lx * qx]
  one[, Lx := fifelse(mx > 0, dx / mx, lx)]
  
  idx_open <- which(one$age_interval_open %in% TRUE)
  if (length(idx_open) == 1) {
    one[idx_open, Lx := lx / pmax(mx, 1e-12)]
  }
  
  one[, Tx := rev(cumsum(rev(Lx)))]
  one[, ex := Tx / pmax(lx, 1e-12)]
  one[, px_5y := 1 - qx]
  
  one[]
}

# ------------------------------------------------------------
# Expandir abridged -> single-age
# ------------------------------------------------------------
expand_life_table_single_age <- function(
    dt_abridged,
    age_open_start = 110L,
    radix = 100000,
    eps_q = 1e-12,
    adult_age_cut = 15L,
    old_age_cut = 85L
) {
  stopifnot(is.data.table(dt_abridged))
  
  needed <- c(
    "year_id",
    "location_id", "location_name",
    "sex_id", "sex_source_value",
    "life_table_label",
    "age_start", "age_interval_width", "age_interval_open",
    "qx",
    "source_year_left", "source_year_right", "interp_weight", "interp_scale"
  )
  
  miss <- setdiff(needed, names(dt_abridged))
  if (length(miss) > 0) {
    stop("Faltan columnas requeridas en dt_abridged: ", paste(miss, collapse = ", "))
  }
  
  dt <- copy(dt_abridged)
  dt[, age_start := as.integer(age_start)]
  dt[, age_interval_width := as.integer(age_interval_width)]
  dt[, qx := clip_qx(qx, eps = eps_q)]
  dt[, age_interval_open := as.logical(age_interval_open)]
  if ("mx" %in% names(dt)) dt[, mx := safe_num(mx)]
  
  group_vars <- c(
    "year_id",
    "location_id", "location_name",
    "sex_id", "sex_source_value",
    "life_table_label",
    "source_year_left", "source_year_right",
    "interp_weight", "interp_scale"
  )
  
  age_last_closed <- age_open_start - 1L
  setorderv(dt, c(group_vars, "age_start"))
  
  out <- dt[, {
    g <- copy(.SD)
    setorder(g, age_start)
    
    shape_obj <- build_single_age_shape(
      g = g,
      age_last_closed = age_last_closed,
      age_open_start = age_open_start,
      eps_q = eps_q,
      adult_age_cut = adult_age_cut,
      old_age_cut = old_age_cut
    )
    
    one <- calibrate_single_age_to_abridged_qx(
      g = g,
      base_shape_obj = shape_obj,
      eps_q = eps_q
    )
    
    one <- reconstruct_single_age_life_table(
      one = one,
      radix = radix,
      eps_q = eps_q
    )
    
    one[]
  }, by = group_vars]
  
  setcolorder(
    out,
    c(
      group_vars,
      "age_start", "age_interval_width", "age_interval_open",
      "mx", "qx", "px_5y", "lx", "dx", "Lx", "Tx", "ex"
    )
  )
  
  out[]
}

# ------------------------------------------------------------
# Reagregar single-age -> abridged
# ------------------------------------------------------------
reaggregate_to_abridged <- function(
    dt_single,
    dt_abridged_ref,
    eps = 1e-12
) {
  stopifnot(is.data.table(dt_single), is.data.table(dt_abridged_ref))
  
  group_vars <- c(
    "year_id",
    "location_id", "location_name",
    "sex_id", "sex_source_value",
    "life_table_label",
    "source_year_left", "source_year_right",
    "interp_weight", "interp_scale"
  )
  
  ref <- copy(dt_abridged_ref)
  ref[, age_start := as.integer(age_start)]
  ref[, age_interval_width := as.integer(age_interval_width)]
  ref[, age_interval_open := as.logical(age_interval_open)]
  
  sin <- copy(dt_single)
  sin[, age_start := as.integer(age_start)]
  
  setorderv(ref, c(group_vars, "age_start"))
  setorderv(sin, c(group_vars, "age_start"))
  
  ref[, age_end := age_start + age_interval_width]
  
  max_single_by_group <- sin[
    ,
    .(max_single_age = max(age_start, na.rm = TRUE)),
    by = group_vars
  ]
  
  ref <- merge(
    ref,
    max_single_by_group,
    by = group_vars,
    all.x = TRUE,
    sort = FALSE
  )
  
  ref[age_interval_open == TRUE, age_end := max_single_age + 1L]
  
  out <- ref[, {
    g <- copy(.SD)
    ss_all <- sin[.BY, on = group_vars, nomatch = 0]
    
    ans <- rbindlist(lapply(seq_len(nrow(g)), function(i) {
      x0 <- g$age_start[i]
      x1 <- g$age_end[i]
      
      ss <- ss_all[age_start >= x0 & age_start < x1]
      
      qx_hat <- 1 - prod(pmax(eps, 1 - ss$qx))
      mx_hat <- sum(ss$dx) / pmax(sum(ss$Lx), eps)
      ex_hat <- ss_all[age_start == x0, ex][1]
      
      data.table(
        age_start = x0,
        age_interval_width = g$age_interval_width[i],
        age_interval_open = g$age_interval_open[i],
        qx_hat = qx_hat,
        mx_hat = mx_hat,
        ex_hat = ex_hat
      )
    }))
    
    ans
  }, by = group_vars]
  
  cmp <- merge(
    ref[, c(group_vars, "age_start", "age_interval_width", "age_interval_open", "qx", "mx", "ex"), with = FALSE],
    out,
    by = c(group_vars, "age_start", "age_interval_width", "age_interval_open"),
    all.x = TRUE,
    sort = FALSE
  )
  
  cmp[, `:=`(
    err_qx = qx_hat - qx,
    err_mx = mx_hat - mx,
    err_ex = ex_hat - ex,
    rel_err_qx = (qx_hat - qx) / pmax(abs(qx), 1e-6),
    rel_err_mx = (mx_hat - mx) / pmax(abs(mx), 1e-6),
    rel_err_ex = (ex_hat - ex) / pmax(abs(ex), 1e-6)
  )]
  
  cmp[, is_open_ref := age_interval_open %in% TRUE]
  cmp[is_open_ref == TRUE, `:=`(
    err_qx = NA_real_,
    rel_err_qx = NA_real_
  )]
  
  cmp[]
}

# ------------------------------------------------------------
# QC específico del intervalo abierto
# ------------------------------------------------------------
qc_open_interval <- function(dt_single) {
  stopifnot(is.data.table(dt_single))
  
  group_vars <- c(
    "year_id",
    "location_id", "location_name",
    "sex_id", "sex_source_value",
    "life_table_label",
    "source_year_left", "source_year_right",
    "interp_weight", "interp_scale"
  )
  
  dt <- copy(dt_single)
  setorderv(dt, c(group_vars, "age_start"))
  
  open_qc <- dt[, {
    last_closed_age <- max(age_start[age_interval_open == FALSE], na.rm = TRUE)
    open_age <- age_start[age_interval_open == TRUE][1]
    
    row_last <- .SD[age_start == last_closed_age]
    row_open <- .SD[age_start == open_age]
    
    mx_last <- row_last$mx[1]
    qx_last <- row_last$qx[1]
    ex_last <- row_last$ex[1]
    
    mx_open <- row_open$mx[1]
    qx_open <- row_open$qx[1]
    ex_open <- row_open$ex[1]
    
    data.table(
      last_closed_age = last_closed_age,
      open_age = open_age,
      mx_last = mx_last,
      mx_open = mx_open,
      qx_last = qx_last,
      qx_open = qx_open,
      ex_last = ex_last,
      ex_open = ex_open,
      mx_jump_ratio = mx_open / pmax(mx_last, 1e-12),
      qx_jump_ratio = qx_open / pmax(qx_last, 1e-12),
      ex_drop = ex_last - ex_open,
      chk_mx_open_pos = is.finite(mx_open) && mx_open > 0,
      chk_qx_open_range = is.finite(qx_open) && qx_open >= 0 && qx_open <= 1,
      chk_ex_open_pos = is.finite(ex_open) && ex_open > 0,
      chk_tail_not_implausible = is.finite(mx_open / pmax(mx_last, 1e-12)) &&
        (mx_open / pmax(mx_last, 1e-12)) >= 0.8 &&
        (mx_open / pmax(mx_last, 1e-12)) <= 5
    )
  }, by = group_vars]
  
  open_qc[, any_open_fail := !(
    chk_mx_open_pos &
      chk_qx_open_range &
      chk_ex_open_pos &
      chk_tail_not_implausible
  )]
  
  open_qc[]
}

# ------------------------------------------------------------
# Validación de la tabla single-age
# ------------------------------------------------------------
validate_life_table_single_age <- function(
    dt_single,
    dt_abridged_ref = NULL,
    qx_gold_tol = 1e-10
) {
  stopifnot(is.data.table(dt_single))
  
  group_vars <- c(
    "year_id",
    "location_id", "location_name",
    "sex_id", "sex_source_value",
    "life_table_label",
    "source_year_left", "source_year_right",
    "interp_weight", "interp_scale"
  )
  
  dt <- copy(dt_single)
  setorderv(dt, c(group_vars, "age_start"))
  
  dt[, chk_qx_range := (qx >= -1e-12 & qx <= 1 + 1e-12)]
  dt[, chk_mx_range := (mx >= -1e-12)]
  dt[, chk_lx_nonneg := (lx >= -1e-8)]
  dt[, chk_ex_nonneg := (ex >= -1e-8)]
  
  dt[, chk_lx_mono := (lx <= shift(lx, type = "lag", fill = Inf) + 1e-8), by = group_vars]
  dt[, chk_Tx_mono := (Tx <= shift(Tx, type = "lag", fill = Inf) + 1e-8), by = group_vars]
  
  dt[, err_dx := abs(dx - lx * qx)]
  dt[, chk_dx := err_dx <= 1e-6]
  
  dt[, lx_next := shift(lx, type = "lead"), by = group_vars]
  dt[, err_lx := abs((lx - dx) - lx_next)]
  dt[is.na(lx_next), err_lx := 0]
  dt[, chk_lx := err_lx <= 1e-6]
  dt[, lx_next := NULL]
  
  dt[, err_ex := abs(ex - Tx / pmax(lx, 1e-12))]
  dt[, chk_ex := err_ex <= 1e-6]
  
  # closed intervals only
  dt[, err_mx_cons := NA_real_]
  dt[, chk_mx_cons := TRUE]
  dt[age_interval_open == FALSE, err_mx_cons := abs(mx - dx / pmax(Lx, 1e-12))]
  dt[age_interval_open == FALSE, chk_mx_cons := err_mx_cons <= 1e-8]
  
  summary_dt <- dt[, .(
    n_rows = .N,
    age_min = min(age_start),
    age_max = max(age_start),
    e0 = ex[age_start == 0][1],
    e60 = ex[age_start == 60][1],
    e80 = ex[age_start == 80][1],
    e109 = ex[age_start == 109][1],
    e110_open = ex[age_interval_open == TRUE][1],
    any_fail = any(
      !chk_qx_range | !chk_mx_range | !chk_lx_nonneg |
        !chk_lx_mono | !chk_Tx_mono | !chk_ex_nonneg |
        !chk_dx | !chk_lx | !chk_ex | !chk_mx_cons
    ),
    n_fail_rows = sum(
      !chk_qx_range | !chk_mx_range | !chk_lx_nonneg |
        !chk_lx_mono | !chk_Tx_mono | !chk_ex_nonneg |
        !chk_dx | !chk_lx | !chk_ex | !chk_mx_cons
    )
  ), by = group_vars]
  
  cmp <- NULL
  if (!is.null(dt_abridged_ref)) {
    cmp <- reaggregate_to_abridged(dt_single = dt, dt_abridged_ref = dt_abridged_ref)
    
    cons <- cmp[is_open_ref == FALSE, .(
      mean_abs_rel_err_qx = mean(abs(rel_err_qx), na.rm = TRUE),
      max_abs_rel_err_qx  = max(abs(rel_err_qx), na.rm = TRUE),
      mean_abs_rel_err_mx = mean(abs(rel_err_mx), na.rm = TRUE),
      max_abs_rel_err_mx  = max(abs(rel_err_mx), na.rm = TRUE),
      mean_abs_rel_err_ex = mean(abs(rel_err_ex), na.rm = TRUE),
      max_abs_rel_err_ex  = max(abs(rel_err_ex), na.rm = TRUE)
    ), by = group_vars]
    
    summary_dt <- merge(summary_dt, cons, by = group_vars, all.x = TRUE)
    summary_dt[, cons_fail := max_abs_rel_err_qx > qx_gold_tol]
  } else {
    summary_dt[, cons_fail := FALSE]
  }
  
  open_qc <- qc_open_interval(dt)
  summary_dt <- merge(
    summary_dt,
    open_qc[, c(group_vars, "mx_jump_ratio", "qx_jump_ratio", "ex_drop", "any_open_fail"), with = FALSE],
    by = group_vars,
    all.x = TRUE
  )
  
  fail_rows <- dt[!(
    chk_qx_range & chk_mx_range & chk_lx_nonneg &
      chk_lx_mono & chk_Tx_mono & chk_ex_nonneg &
      chk_dx & chk_lx & chk_ex & chk_mx_cons
  )]
  
  list(
    summary = summary_dt,
    failures = fail_rows,
    reagg_cmp = cmp,
    open_qc = open_qc
  )
}

# ------------------------------------------------------------
# 1) Leer input abridged anual
# ------------------------------------------------------------
cat("==> Leyendo annual abridged: ", CFG$in_file, "\n", sep = "")
dt_abr <- fread(CFG$in_file)

if ("age_start" %in% names(dt_abr)) dt_abr[, age_start := as.integer(age_start)]
if ("age_interval_width" %in% names(dt_abr)) dt_abr[, age_interval_width := as.integer(age_interval_width)]
if ("year_id" %in% names(dt_abr)) dt_abr[, year_id := as.integer(year_id)]
if ("sex_id" %in% names(dt_abr)) dt_abr[, sex_id := as.integer(sex_id)]
if ("location_id" %in% names(dt_abr)) dt_abr[, location_id := as.integer(location_id)]
if ("qx" %in% names(dt_abr)) dt_abr[, qx := safe_num(qx)]
if ("mx" %in% names(dt_abr)) dt_abr[, mx := safe_num(mx)]
if ("age_interval_open" %in% names(dt_abr)) dt_abr[, age_interval_open := as.logical(age_interval_open)]

# ------------------------------------------------------------
# 2) Expandir a edad simple
# ------------------------------------------------------------
cat("==> Expandiendo a edad simple cerrada 0-", CFG$age_last_closed, " y abierto ", CFG$age_open_start, "+ ...\n", sep = "")

dt_single <- expand_life_table_single_age(
  dt_abridged = as.data.table(dt_abr),
  age_open_start = CFG$age_open_start,
  radix = CFG$radix,
  eps_q = CFG$eps_q,
  adult_age_cut = CFG$adult_age_cut,
  old_age_cut = CFG$old_age_cut
)

# ------------------------------------------------------------
# 3) QC y validación
# ------------------------------------------------------------
cat("==> Validando coherencia interna y conservación de qx...\n")

qc <- validate_life_table_single_age(
  dt_single = dt_single,
  dt_abridged_ref = as.data.table(dt_abr),
  qx_gold_tol = CFG$qx_gold_tol
)

qc_sum  <- qc$summary
qc_fail <- qc$failures
qc_reag <- qc$reagg_cmp
qc_open <- qc$open_qc

fwrite(qc_sum,  file.path(CFG$qc_dir, "qc_life_table_single_age_summary.csv"))
fwrite(qc_fail, file.path(CFG$qc_dir, "qc_life_table_single_age_failures.csv"))
if (!is.null(qc_reag)) {
  fwrite(qc_reag, file.path(CFG$qc_dir, "qc_life_table_single_age_reaggregate_compare.csv"))
}
if (!is.null(qc_open)) {
  fwrite(qc_open, file.path(CFG$qc_dir, "qc_life_table_single_age_open_interval.csv"))
}

cat("==> Resumen QC por grupo:\n")
qc_print(qc_sum)

n_fail_groups <- qc_sum[any_fail == TRUE | cons_fail == TRUE | any_open_fail == TRUE, .N]
if (n_fail_groups > 0) {
  cat("==> Ejemplos de grupos FAIL:\n")
  print(qc_sum[any_fail == TRUE | cons_fail == TRUE | any_open_fail == TRUE][1:min(20, .N)])
  stop("Se detectaron FAIL en script 03. Revisa los qc_* antes de propagar.")
}

# ------------------------------------------------------------
# 4) Validación básica reutilizando helper del proyecto
# ------------------------------------------------------------
qc_basic <- qc_life_table_basic(as.data.table(dt_single))
write_qc_list(qc_basic, file.path(CFG$qc_dir, "qc_basic"))

# ------------------------------------------------------------
# 5) Validación contra spec single-age
# ------------------------------------------------------------
cat("==> Validando contra spec single-age...\n")

spec <- read_spec(CFG$spec_file)
validate_by_spec(as.data.table(dt_single), spec)

# ------------------------------------------------------------
# 6) Diccionario extendido
# ------------------------------------------------------------
cat("==> Construyendo diccionario extendido...\n")

dict <- dict_from_spec(
  spec = spec,
  dataset_version = "v1.4.0",
  run_id = NA_character_,
  config_dir = P$CONFIG_DIR
)

dict_ext <- enrich_dict_with_stats(dict, as.data.table(dt_single))

# ------------------------------------------------------------
# 7) Export
# ------------------------------------------------------------
out_csv  <- file.path(CFG$out_dir, CFG$out_file)
out_dict <- file.path(CFG$out_dir, CFG$out_dict)

fwrite(as.data.table(dt_single), out_csv)
fwrite(as.data.table(dict_ext), out_dict)

cat("==> OK export dataset: ", out_csv, "\n", sep = "")
cat("==> OK export dictionary: ", out_dict, "\n", sep = "")
cat("==> OK QC dir: ", CFG$qc_dir, "\n", sep = "")
cat("==> Script 03 finalizado correctamente.\n")