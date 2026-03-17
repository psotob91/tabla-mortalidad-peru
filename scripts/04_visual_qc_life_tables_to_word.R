# scripts/04_visual_qc_qx_only_life_tables_to_pdf.R
# ------------------------------------------------------------
# Script 04
# QC gráfico SOLO de qx
# - Script 02: annual abridged vs anclas empíricas
# - Script 03:
#     A) QC de forma: single-age 1qx vs annual-equivalent 1qx derivado del abridged
#     B) QC de conservación: nqx abridged objetivo vs nqx reagregado desde single-age
# - Salida: PDFs multipágina
# - Vistas:
#     * escala log10
#     * escala lineal
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(grDevices)
  library(grid)
})

# ------------------------------------------------------------
# Sources del proyecto
# ------------------------------------------------------------
source(here("R/io_utils.R"))
source(here("config/parametros-proyecto.R"))

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------
P <- paths_mort()

CFG <- list(
  in_empirical_abridged = file.path(P$FINAL_DIR, "ref_life_table_mortality.csv"),
  in_annual_abridged    = file.path(P$FINAL_DIR, "annual_abridged", "ref_life_table_mortality_annual_abridged.csv"),
  in_single_age         = file.path(P$FINAL_DIR, "single_age", "ref_life_table_mortality_single_age.csv"),
  
  out_dir               = here("reports", "visual_qc_life_table_mortality"),
  metric                = "qx",
  years_target          = if (exists("TARGET_YEARS_LIFE_TABLE")) TARGET_YEARS_LIFE_TABLE else 2018:2024,
  location_ids          = 0:25,
  
  pdf_width             = 14,
  pdf_height            = 8.5,
  base_size             = 11,
  eps_q                 = 1e-12
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Validaciones de existencia
# ------------------------------------------------------------
req_files <- c(
  CFG$in_empirical_abridged,
  CFG$in_annual_abridged,
  CFG$in_single_age
)

miss_files <- req_files[!file.exists(req_files)]
if (length(miss_files) > 0) {
  stop("Faltan archivos de entrada:\n- ", paste(miss_files, collapse = "\n- "))
}

# ------------------------------------------------------------
# Lectura
# ------------------------------------------------------------
emp <- fread(CFG$in_empirical_abridged)
abr <- fread(CFG$in_annual_abridged)
sng <- fread(CFG$in_single_age)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
safe_num <- function(x) suppressWarnings(as.numeric(x))

sex_label_from_id <- function(x) {
  fifelse(
    x == 8507L, "Hombres",
    fifelse(x == 8532L, "Mujeres", paste0("sex_id=", x))
  )
}

metric_label <- function(v) {
  switch(
    v,
    qx = "Probabilidad de morir (qx)",
    v
  )
}

scale_label <- function(scale_mode) {
  if (scale_mode == "log10") "escala log10" else "escala lineal"
}

clip_qx <- function(x, eps = 1e-12) {
  pmin(1 - eps, pmax(0, safe_num(x)))
}

clean_dt <- function(dt) {
  x <- copy(as.data.table(dt))
  
  int_cols <- intersect(
    c(
      "year_id", "period_start_year", "period_end_year", "age_start",
      "age_interval_width", "location_id", "sex_id",
      "source_year_left", "source_year_right"
    ),
    names(x)
  )
  for (cc in int_cols) x[, (cc) := as.integer(get(cc))]
  
  num_cols <- intersect(c("qx", "mx", "ex", "lx", "interp_weight", "dx", "Lx"), names(x))
  for (cc in num_cols) x[, (cc) := safe_num(get(cc))]
  
  if ("location_name" %in% names(x)) x[, location_name := as.character(location_name)]
  if ("sex_source_value" %in% names(x)) x[, sex_source_value := as.character(sex_source_value)]
  if ("life_table_label" %in% names(x)) x[, life_table_label := as.character(life_table_label)]
  if ("interp_scale" %in% names(x)) x[, interp_scale := as.character(interp_scale)]
  if ("age_interval_open" %in% names(x)) x[, age_interval_open := as.logical(age_interval_open)]
  
  x[, sex_label := sex_label_from_id(sex_id)]
  x[, sex_label := factor(sex_label, levels = c("Hombres", "Mujeres"))]
  x[]
}

emp <- clean_dt(emp)
abr <- clean_dt(abr)
sng <- clean_dt(sng)

loc_map <- unique(abr[, .(location_id, location_name)])
setorder(loc_map, location_id)
loc_map <- loc_map[location_id %in% CFG$location_ids]

prep_qx_for_plot <- function(dt, qx_col = "qx", scale_mode = c("log10", "linear")) {
  scale_mode <- match.arg(scale_mode)
  x <- copy(dt)
  
  if (!(qx_col %in% names(x))) return(x)
  
  if (scale_mode == "log10") {
    x[get(qx_col) <= 0 | is.na(get(qx_col)), (qx_col) := NA_real_]
  } else {
    x[, (qx_col) := safe_num(get(qx_col))]
  }
  
  x
}

count_bad_qx_for_log <- function(dt, qx_col = "qx") {
  if (!(qx_col %in% names(dt))) return(0L)
  sum(!is.na(dt[[qx_col]]) & dt[[qx_col]] <= 0)
}

apply_qx_scale <- function(p, scale_mode = c("log10", "linear")) {
  scale_mode <- match.arg(scale_mode)
  
  if (scale_mode == "log10") {
    p + scale_y_log10(labels = label_number(accuracy = 0.001))
  } else {
    p + scale_y_continuous(labels = label_number(accuracy = 0.001))
  }
}

qc_theme <- function() {
  theme_minimal(base_size = CFG$base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9.2),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 8.8),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.2),
      panel.grid.major.y = element_line(linewidth = 0.2),
      plot.margin = margin(8, 18, 8, 8)
    )
}

subtitle_scale_note <- function(scale_mode, n_bad_total) {
  if (scale_mode == "log10" && n_bad_total > 0) {
    paste0(" Vista en escala log10. Se omitieron ", n_bad_total, " valores <= 0 para permitir la transformación.")
  } else if (scale_mode == "linear") {
    " Vista en escala lineal sin transformar para detectar valores negativos, qx > 1, picos o discontinuidades raras."
  } else {
    ""
  }
}

get_scale_modes <- function() c("log10", "linear")

# ------------------------------------------------------------
# Datos de contraste para Script 02
# annual abridged vs anclas empíricas
# ------------------------------------------------------------
get_empirical_overlay_for_location <- function(loc_id, sex_ids = c(8507L, 8532L)) {
  da <- abr[location_id == loc_id & sex_id %in% sex_ids]
  if (nrow(da) == 0) return(data.table())
  
  map_yr <- unique(
    da[, .(
      year_id,
      source_year_left,
      source_year_right,
      sex_id,
      sex_label,
      life_table_label
    )]
  )
  
  emp_sub <- emp[
    location_id == loc_id & sex_id %in% sex_ids,
    .(
      location_id, location_name,
      period_start_year,
      sex_id, sex_label,
      life_table_label,
      age_start, age_interval_width, age_interval_open,
      qx
    )
  ]
  
  left_dt <- merge(
    map_yr[, .(
      year_id,
      source_year_anchor = source_year_left,
      anchor_side = "Empírico izquierdo",
      sex_id, sex_label, life_table_label
    )],
    emp_sub,
    by.x = c("source_year_anchor", "sex_id", "sex_label", "life_table_label"),
    by.y = c("period_start_year", "sex_id", "sex_label", "life_table_label"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  right_dt <- merge(
    map_yr[, .(
      year_id,
      source_year_anchor = source_year_right,
      anchor_side = "Empírico derecho",
      sex_id, sex_label, life_table_label
    )],
    emp_sub,
    by.x = c("source_year_anchor", "sex_id", "sex_label", "life_table_label"),
    by.y = c("period_start_year", "sex_id", "sex_label", "life_table_label"),
    all.x = TRUE,
    allow.cartesian = TRUE
  )
  
  out <- rbind(left_dt, right_dt, fill = TRUE)
  out <- unique(out)
  out[, anchor_side := factor(anchor_side, levels = c("Empírico izquierdo", "Empírico derecho"))]
  out[]
}

# ------------------------------------------------------------
# Helpers específicos Script 03
# ------------------------------------------------------------

# A) referencia anual equivalente desde abridged
get_abridged_annual_equivalent_for_location <- function(loc_id) {
  db <- copy(abr[location_id == loc_id])
  if (nrow(db) == 0) return(data.table())
  
  db[, qx_eq_1y := fifelse(
    age_interval_open %in% TRUE,
    qx,
    1 - (1 - clip_qx(qx, eps = CFG$eps_q))^(1 / pmax(age_interval_width, 1))
  )]
  
  db[]
}

# B) reagregación single-age -> abridged para comparar nqx
reaggregate_single_to_abridged_for_location <- function(loc_id) {
  ds <- copy(sng[location_id == loc_id])
  db <- copy(abr[location_id == loc_id])
  
  if (nrow(ds) == 0 || nrow(db) == 0) return(data.table())
  
  group_vars <- c(
    "year_id",
    "location_id", "location_name",
    "sex_id", "sex_label",
    "sex_source_value",
    "life_table_label",
    "source_year_left", "source_year_right",
    "interp_weight", "interp_scale"
  )
  
  setorderv(ds, c(group_vars, "age_start"))
  setorderv(db, c(group_vars, "age_start"))
  
  max_single <- ds[, .(max_single_age = max(age_start, na.rm = TRUE)), by = group_vars]
  db <- merge(db, max_single, by = group_vars, all.x = TRUE, sort = FALSE)
  
  db[, age_end := age_start + age_interval_width]
  db[age_interval_open %in% TRUE, age_end := max_single_age + 1L]
  
  out <- db[, {
    g <- copy(.SD)
    ss_all <- ds[.BY, on = group_vars, nomatch = 0]
    
    ans <- rbindlist(lapply(seq_len(nrow(g)), function(i) {
      x0 <- g$age_start[i]
      x1 <- g$age_end[i]
      
      ss <- ss_all[age_start >= x0 & age_start < x1]
      qx_hat <- 1 - prod(pmax(CFG$eps_q, 1 - ss$qx))
      
      data.table(
        age_start = x0,
        age_interval_width = g$age_interval_width[i],
        age_interval_open = g$age_interval_open[i],
        qx_hat_abr = qx_hat,
        qx_target_abr = g$qx[i]
      )
    }))
    
    ans
  }, by = group_vars]
  
  out[]
}

# ------------------------------------------------------------
# Script 02
# annual abridged vs anclas empíricas
# ------------------------------------------------------------
plot_script02_qx <- function(loc_id, scale_mode = c("log10", "linear")) {
  scale_mode <- match.arg(scale_mode)
  
  da0 <- abr[location_id == loc_id]
  ov0 <- get_empirical_overlay_for_location(loc_id)
  
  if (nrow(da0) == 0) return(NULL)
  
  loc_name <- da0$location_name[1]
  
  da <- prep_qx_for_plot(da0, qx_col = "qx", scale_mode = scale_mode)
  ov <- prep_qx_for_plot(ov0, qx_col = "qx", scale_mode = scale_mode)
  
  da[, year_f := factor(year_id, levels = CFG$years_target)]
  ov[, year_f := factor(year_id, levels = CFG$years_target)]
  
  n_bad_total <- count_bad_qx_for_log(da0, "qx") + count_bad_qx_for_log(ov0, "qx")
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, color = "grey50") +
    geom_line(
      data = da,
      aes(x = age_start, y = qx, group = 1),
      linewidth = 0.75,
      color = "#1F4E79",
      alpha = 0.95,
      na.rm = TRUE
    ) +
    geom_point(
      data = da,
      aes(x = age_start, y = qx),
      size = 0.9,
      color = "#1F4E79",
      alpha = 0.8,
      na.rm = TRUE
    ) +
    geom_point(
      data = ov,
      aes(x = age_start, y = qx, shape = anchor_side),
      size = 1.5,
      stroke = 0.7,
      color = "black",
      fill = "white",
      na.rm = TRUE
    ) +
    facet_grid(sex_label ~ year_f, scales = "free_y") +
    scale_shape_manual(values = c(1, 17), drop = FALSE) +
    labs(
      title = paste0(loc_name, " | Script 02 | ", metric_label(CFG$metric), " | ", scale_label(scale_mode)),
      subtitle = paste0(
        "Línea azul: tabla annual abridged interpolada. ",
        "Símbolos negros: anclas empíricas izquierda y derecha usadas para la interpolación.",
        subtitle_scale_note(scale_mode, n_bad_total)
      ),
      x = "Edad de inicio del intervalo",
      y = metric_label(CFG$metric),
      shape = "Serie"
    ) +
    qc_theme()
  
  apply_qx_scale(p, scale_mode = scale_mode)
}

# ------------------------------------------------------------
# Script 03A
# QC de forma:
# single-age 1qx vs annual-equivalent 1qx derivado del abridged
# ------------------------------------------------------------
plot_script03_qx_shape <- function(loc_id, scale_mode = c("log10", "linear")) {
  scale_mode <- match.arg(scale_mode)
  
  ds0 <- sng[location_id == loc_id]
  db0 <- get_abridged_annual_equivalent_for_location(loc_id)
  
  if (nrow(ds0) == 0 || nrow(db0) == 0) return(NULL)
  
  loc_name <- ds0$location_name[1]
  
  ds <- prep_qx_for_plot(ds0, qx_col = "qx", scale_mode = scale_mode)
  db <- prep_qx_for_plot(db0, qx_col = "qx_eq_1y", scale_mode = scale_mode)
  
  ds[, year_f := factor(year_id, levels = CFG$years_target)]
  db[, year_f := factor(year_id, levels = CFG$years_target)]
  
  n_bad_total <- count_bad_qx_for_log(ds0, "qx") + count_bad_qx_for_log(db0, "qx_eq_1y")
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, color = "grey50") +
    geom_step(
      data = ds,
      aes(x = age_start, y = qx, group = 1),
      linewidth = 0.72,
      color = "#8B1E3F",
      alpha = 0.95,
      direction = "hv",
      na.rm = TRUE
    ) +
    geom_point(
      data = db,
      aes(x = age_start, y = qx_eq_1y),
      size = 1.5,
      shape = 21,
      stroke = 0.6,
      color = "black",
      fill = "white",
      na.rm = TRUE
    ) +
    facet_grid(sex_label ~ year_f, scales = "free_y") +
    labs(
      title = paste0(loc_name, " | Script 03A | ", metric_label(CFG$metric), " | ", scale_label(scale_mode)),
      subtitle = paste0(
        "Línea vino: 1qx de la tabla a edad simple. ",
        "Puntos negros: 1qx equivalente anual derivado de la tabla abridged.",
        subtitle_scale_note(scale_mode, n_bad_total)
      ),
      x = "Edad exacta / edad de inicio",
      y = metric_label(CFG$metric)
    ) +
    qc_theme()
  
  apply_qx_scale(p, scale_mode = scale_mode)
}

# ------------------------------------------------------------
# Script 03B
# QC de conservación:
# nqx abridged objetivo vs nqx reagregado desde single-age
# ------------------------------------------------------------
plot_script03_qx_reaggregated <- function(loc_id, scale_mode = c("log10", "linear")) {
  scale_mode <- match.arg(scale_mode)
  
  rg0 <- reaggregate_single_to_abridged_for_location(loc_id)
  if (nrow(rg0) == 0) return(NULL)
  
  loc_name <- rg0$location_name[1]
  
  rg <- copy(rg0)
  rg <- prep_qx_for_plot(rg, qx_col = "qx_hat_abr", scale_mode = scale_mode)
  rg <- prep_qx_for_plot(rg, qx_col = "qx_target_abr", scale_mode = scale_mode)
  
  rg[, year_f := factor(year_id, levels = CFG$years_target)]
  
  n_bad_total <- count_bad_qx_for_log(rg0, "qx_hat_abr") + count_bad_qx_for_log(rg0, "qx_target_abr")
  
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25, color = "grey50") +
    geom_line(
      data = rg,
      aes(x = age_start, y = qx_hat_abr, group = 1),
      linewidth = 0.75,
      color = "#8B1E3F",
      alpha = 0.95,
      na.rm = TRUE
    ) +
    geom_point(
      data = rg,
      aes(x = age_start, y = qx_hat_abr),
      size = 1.0,
      color = "#8B1E3F",
      alpha = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      data = rg,
      aes(x = age_start, y = qx_target_abr),
      size = 1.8,
      shape = 21,
      stroke = 0.7,
      color = "black",
      fill = "white",
      na.rm = TRUE
    ) +
    facet_grid(sex_label ~ year_f, scales = "free_y") +
    labs(
      title = paste0(loc_name, " | Script 03B | ", metric_label(CFG$metric), " | ", scale_label(scale_mode)),
      subtitle = paste0(
        "Línea/puntos vino: nqx reagregado desde la tabla single-age. ",
        "Puntos negros: nqx objetivo de la tabla annual abridged.",
        subtitle_scale_note(scale_mode, n_bad_total)
      ),
      x = "Edad de inicio del intervalo abridged",
      y = metric_label(CFG$metric)
    ) +
    qc_theme()
  
  apply_qx_scale(p, scale_mode = scale_mode)
}

# ------------------------------------------------------------
# PDF helpers
# ------------------------------------------------------------
open_pdf_device <- function(file_pdf) {
  dir.create(dirname(file_pdf), recursive = TRUE, showWarnings = FALSE)
  cairo_pdf(
    filename = file_pdf,
    width = CFG$pdf_width,
    height = CFG$pdf_height,
    onefile = TRUE
  )
}

close_pdf_device <- function() {
  dev.off()
}

plot_pdf_cover <- function(title_main, subtitle_main) {
  grid.newpage()
  pushViewport(viewport())
  grid.text(
    title_main,
    x = 0.5, y = 0.72,
    gp = gpar(fontsize = 22, fontface = "bold")
  )
  grid.text(
    subtitle_main,
    x = 0.5, y = 0.62,
    gp = gpar(fontsize = 12)
  )
  grid.text(
    paste0("Fecha de generación: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    x = 0.5, y = 0.53,
    gp = gpar(fontsize = 10)
  )
  grid.text(
    paste0(
      "Cobertura: location_id 0:25 | años ",
      min(CFG$years_target), "-", max(CFG$years_target)
    ),
    x = 0.5, y = 0.47,
    gp = gpar(fontsize = 10)
  )
  popViewport()
}

plot_interpretation_page_script02 <- function() {
  txt_title <- "Guía breve de interpretación visual | Script 02"
  
  txt_body <- paste(
    "1) La vista log10 sirve para comparar el patrón relativo de qx entre edades, sexos y años.",
    "2) La vista lineal sirve para detectar qx negativa, qx > 1, picos aislados o discontinuidades bruscas.",
    "3) qx debe cumplir 0 <= qx <= 1.",
    "4) En población humana, qx suele ser relativamente alta en edad 0, cae en infancia, alcanza un mínimo en niñez/adolescencia temprana y luego aumenta con la edad.",
    "5) La línea interpolada debe ser coherente con las anclas empíricas izquierda y derecha.",
    sep = "\n\n"
  )
  
  grid.newpage()
  pushViewport(viewport())
  grid.text(
    txt_title,
    x = 0.5, y = 0.92,
    gp = gpar(fontsize = 18, fontface = "bold")
  )
  grid.text(
    txt_body,
    x = 0.06, y = 0.82,
    just = c("left", "top"),
    gp = gpar(fontsize = 11)
  )
  popViewport()
}

plot_interpretation_page_script03 <- function() {
  txt_title <- "Guía breve de interpretación visual | Script 03"
  
  txt_body <- paste(
    "1) Script 03A evalúa la FORMA de 1qx a edad simple comparándola contra un 1qx equivalente anual derivado del abridged. Aquí ambas series sí son comparables en magnitud.",
    "2) Script 03B evalúa la CONSERVACIÓN de qx: compara el nqx abridged objetivo contra el nqx reagregado desde la tabla single-age.",
    "3) Lo importante en 03A es la plausibilidad demográfica de la forma por edad: caída infantil, mínimo en edades tempranas y aumento progresivo posterior.",
    "4) Lo importante en 03B es que la línea vino y los puntos negros casi se superpongan; eso indica que la expansión single-age respeta el qx abridged objetivo.",
    "5) qx debe cumplir 0 <= qx <= 1 en todas las vistas.",
    sep = "\n\n"
  )
  
  grid.newpage()
  pushViewport(viewport())
  grid.text(
    txt_title,
    x = 0.5, y = 0.92,
    gp = gpar(fontsize = 18, fontface = "bold")
  )
  grid.text(
    txt_body,
    x = 0.06, y = 0.82,
    just = c("left", "top"),
    gp = gpar(fontsize = 11)
  )
  popViewport()
}

# ------------------------------------------------------------
# PDF Script 02
# ------------------------------------------------------------
build_pdf_script02 <- function() {
  out_pdf <- file.path(CFG$out_dir, "visual_qc_script02_qx_annual_abridged_vs_empirical.pdf")
  
  idx <- data.table()
  
  open_pdf_device(out_pdf)
  on.exit(close_pdf_device(), add = TRUE)
  
  plot_pdf_cover(
    title_main = "Inspección visual | Script 02 | qx",
    subtitle_main = "Annual abridged vs anclas empíricas"
  )
  plot_interpretation_page_script02()
  
  for (loc_id in loc_map$location_id) {
    loc_name <- loc_map[location_id == loc_id, location_name][1]
    
    for (scale_mode in get_scale_modes()) {
      p <- plot_script02_qx(loc_id = loc_id, scale_mode = scale_mode)
      
      if (!is.null(p)) {
        print(p)
        
        idx <- rbind(
          idx,
          data.table(
            section = "script02",
            location_id = loc_id,
            location_name = loc_name,
            metric = "qx",
            scale_mode = scale_mode,
            pdf_file = normalizePath(out_pdf, winslash = "/", mustWork = FALSE)
          ),
          fill = TRUE
        )
      }
    }
  }
  
  idx
}

# ------------------------------------------------------------
# PDF Script 03 - FORMA
# ------------------------------------------------------------
build_pdf_script03_shape <- function() {
  out_pdf <- file.path(CFG$out_dir, "visual_qc_script03A_qx_shape_single_age_vs_annual_equivalent.pdf")
  
  idx <- data.table()
  
  open_pdf_device(out_pdf)
  on.exit(close_pdf_device(), add = TRUE)
  
  plot_pdf_cover(
    title_main = "Inspección visual | Script 03A | qx",
    subtitle_main = "Forma: single-age 1qx vs annual-equivalent 1qx"
  )
  plot_interpretation_page_script03()
  
  for (loc_id in loc_map$location_id) {
    loc_name <- loc_map[location_id == loc_id, location_name][1]
    
    for (scale_mode in get_scale_modes()) {
      p <- plot_script03_qx_shape(loc_id = loc_id, scale_mode = scale_mode)
      
      if (!is.null(p)) {
        print(p)
        
        idx <- rbind(
          idx,
          data.table(
            section = "script03A_shape",
            location_id = loc_id,
            location_name = loc_name,
            metric = "qx",
            scale_mode = scale_mode,
            pdf_file = normalizePath(out_pdf, winslash = "/", mustWork = FALSE)
          ),
          fill = TRUE
        )
      }
    }
  }
  
  idx
}

# ------------------------------------------------------------
# PDF Script 03 - CONSERVACIÓN
# ------------------------------------------------------------
build_pdf_script03_reagg <- function() {
  out_pdf <- file.path(CFG$out_dir, "visual_qc_script03B_qx_reaggregated_vs_abridged.pdf")
  
  idx <- data.table()
  
  open_pdf_device(out_pdf)
  on.exit(close_pdf_device(), add = TRUE)
  
  plot_pdf_cover(
    title_main = "Inspección visual | Script 03B | qx",
    subtitle_main = "Conservación: nqx reagregado vs nqx abridged objetivo"
  )
  plot_interpretation_page_script03()
  
  for (loc_id in loc_map$location_id) {
    loc_name <- loc_map[location_id == loc_id, location_name][1]
    
    for (scale_mode in get_scale_modes()) {
      p <- plot_script03_qx_reaggregated(loc_id = loc_id, scale_mode = scale_mode)
      
      if (!is.null(p)) {
        print(p)
        
        idx <- rbind(
          idx,
          data.table(
            section = "script03B_reagg",
            location_id = loc_id,
            location_name = loc_name,
            metric = "qx",
            scale_mode = scale_mode,
            pdf_file = normalizePath(out_pdf, winslash = "/", mustWork = FALSE)
          ),
          fill = TRUE
        )
      }
    }
  }
  
  idx
}

# ------------------------------------------------------------
# Ejecutar
# ------------------------------------------------------------
cat("==> Generando PDFs de visual QC SOLO para qx...\n")

idx02   <- build_pdf_script02()
idx03A  <- build_pdf_script03_shape()
idx03B  <- build_pdf_script03_reagg()

pdf_index <- data.table(
  section = c("script02", "script03A_shape", "script03B_reagg"),
  pdf_file = normalizePath(
    c(
      file.path(CFG$out_dir, "visual_qc_script02_qx_annual_abridged_vs_empirical.pdf"),
      file.path(CFG$out_dir, "visual_qc_script03A_qx_shape_single_age_vs_annual_equivalent.pdf"),
      file.path(CFG$out_dir, "visual_qc_script03B_qx_reaggregated_vs_abridged.pdf")
    ),
    winslash = "/",
    mustWork = FALSE
  )
)

plot_index <- rbindlist(list(idx02, idx03A, idx03B), fill = TRUE)

fwrite(pdf_index, file.path(CFG$out_dir, "visual_qc_qx_pdf_documents_index.csv"))
fwrite(plot_index, file.path(CFG$out_dir, "visual_qc_qx_pdf_plots_index.csv"))

cat("==> OK PDF script 02: ",
    file.path(CFG$out_dir, "visual_qc_script02_qx_annual_abridged_vs_empirical.pdf"),
    "\n", sep = "")
cat("==> OK PDF script 03A: ",
    file.path(CFG$out_dir, "visual_qc_script03A_qx_shape_single_age_vs_annual_equivalent.pdf"),
    "\n", sep = "")
cat("==> OK PDF script 03B: ",
    file.path(CFG$out_dir, "visual_qc_script03B_qx_reaggregated_vs_abridged.pdf"),
    "\n", sep = "")
cat("==> OK carpeta de salida: ", CFG$out_dir, "\n", sep = "")
cat("==> Script 04 finalizado correctamente.\n")