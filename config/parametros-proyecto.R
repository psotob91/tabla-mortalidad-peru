# config/parametros-proyecto.R
# ---------------------------------------------
# Parámetros modificables del proyecto
# ---------------------------------------------

# Años de análisis principales
# (Modificar aquí si cambia la longitud de seguimiento)
YEARS_ANALISIS <- 2018:2024

# Chequeo suave: recomendar mínimo 5 años consecutivos
if (length(YEARS_ANALISIS) < 5) {
  warning(
    "Se recomienda usar al menos 5 años consecutivos en YEARS_ANALISIS ",
    "para poder modelar adecuadamente la tendencia."
  )
}

TARGET_YEARS_LIFE_TABLE <- 2018:2024
LIFE_TABLE_INTERP_SCALE <- "cloglog"
LIFE_TABLE_ANCHOR_YEAR_VAR <- "period_start_year"
