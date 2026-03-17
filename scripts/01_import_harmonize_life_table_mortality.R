# scripts/01_import_harmonize_life_table_mortality.R

library(here)
library(readxl)
library(dplyr)
library(stringr)
library(janitor)
library(data.table)
library(stringi)

source(here("R/io_utils.R"))
source(here("R/spec_utils.R"))
source(here("R/dictionary_utils.R"))
source(here("R/qc_life_table_utils.R"))

DROP_LIMA_METROPOLITANA <- TRUE

P <- paths_mort()

in_xlsx <- file.path(P$RAW_DIR, "mortalidad_combinada_peru.xlsx")
if (!file.exists(in_xlsx)) stop("No existe: ", in_xlsx)

spec_path <- file.path(P$CONFIG_DIR, "spec_life_table_mortality.yml")
if (!file.exists(spec_path)) stop("No existe: ", spec_path)

maestro_loc_path <- file.path(P$CONFIG_DIR, "maestro_location_dept.csv")
if (!file.exists(maestro_loc_path)) stop("No existe: ", maestro_loc_path)

maestro_loc <- data.table::fread(maestro_loc_path)

norm_key <- function(x) {
  x |>
    as.character() |>
    str_trim() |>
    stringi::stri_trans_general("Latin-ASCII") |>
    str_to_upper()
}

recode_dept <- function(x) {
  y <- norm_key(x)
  
  y <- y |>
    str_replace("^DEPARTAMENTO\\s+DE\\s+", "") |>
    str_replace("^DEPARTAMENTO\\s+DEL\\s+", "") |>
    str_replace("^DPTO\\.?\\s+DE\\s+", "") |>
    str_replace("^REGION\\s+DE\\s+", "") |>
    str_replace("^GOBIERNO\\s+REGIONAL\\s+DE\\s+", "") |>
    str_squish()
  
  y <- dplyr::case_when(
    y %in% c("CALLAO", "PROV. CONST. DEL CALLAO") ~ "PROVINCIA CONSTITUCIONAL DEL CALLAO",
    TRUE ~ y
  )
  
  y
}

# ----------------------------
# 1) Import
# ----------------------------
df_raw <- readxl::read_excel(in_xlsx, sheet = 1) |>
  janitor::clean_names()

# ----------------------------
# 2) Detectar columna P(x,x+5)
# ----------------------------
px_col <- intersect(
  c("p_x_x_5", "px_x_5", "p_x_x_5y", "p_x_x_5_years"),
  names(df_raw)
)[1]

if (is.na(px_col)) {
  stop(
    "No encuentro columna de P(x,x+5) tras clean_names().\n",
    "Columnas disponibles:\n- ",
    paste(names(df_raw), collapse = "\n- ")
  )
}

# ----------------------------
# 3) Excluir Lima Metropolitana
# ----------------------------
df_stage <- df_raw |>
  mutate(dept_src_norm = norm_key(departamento))

if (DROP_LIMA_METROPOLITANA) {
  dropped <- df_stage |>
    filter(dept_src_norm == "LIMA METROPOLITANA")
  
  if (nrow(dropped) > 0) {
    dir.create(P$QC_DIR, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      as.data.table(dropped),
      file.path(P$QC_DIR, "qc_dropped_lima_metropolitana_rows.csv")
    )
    message(
      "⚠️ Se excluyeron filas de LIMA METROPOLITANA: ", nrow(dropped),
      " (ver QC: qc_dropped_lima_metropolitana_rows.csv)"
    )
  }
  
  df_stage <- df_stage |>
    filter(dept_src_norm != "LIMA METROPOLITANA")
}

# ----------------------------
# 4) Harmonización base
# ----------------------------
df_clean <- df_stage |>
  mutate(
    n = suppressWarnings(as.integer(as.numeric(n))),
    px_5y = .data[[px_col]] |>
      as.character() |>
      str_replace(",", ".") |>
      str_extract("[0-9.]+") |>
      as.numeric(),
    period_label = as.character(periodo),
    period_start_year = as.integer(str_sub(period_label, 1, 4)),
    period_end_year   = as.integer(str_sub(period_label, 6, 9)),
    location_name = as.character(departamento),
    sex_source_value = as.character(sexo),
    life_table_label = as.character(tabla)
  ) |>
  rename(
    age_start          = edad_x,
    age_interval_width = n,
    mx                 = m_x_n,
    qx                 = q_x_n,
    lx                 = l_x,
    dx                 = d_x_n,
    Lx                 = l_x_n,
    Tx                 = t_x,
    ex                 = e_x
  ) |>
  mutate(
    across(c(mx, qx, lx, dx, Lx, Tx, ex), as.numeric),
    age_start = as.integer(age_start),
    age_interval_width = as.integer(age_interval_width)
  ) |>
  mutate(sex_key = norm_key(sex_source_value)) |>
  mutate(
    sex_id = case_when(
      sex_key %in% c("HOMBRE", "MASCULINO", "MALE", "M") ~ 8507L,
      sex_key %in% c("MUJER", "FEMENINO", "FEMALE", "F") ~ 8532L,
      TRUE ~ NA_integer_
    )
  ) |>
  mutate(dept_key = recode_dept(location_name)) |>
  left_join(
    maestro_loc |>
      mutate(dept_key = norm_key(location_name)) |>
      select(location_id, dept_key),
    by = "dept_key"
  ) |>
  select(-sex_key, -dept_key, -dept_src_norm) |>
  select(
    period_label, period_start_year, period_end_year,
    location_id, location_name,
    sex_id, sex_source_value,
    life_table_label,
    age_start, age_interval_width,
    mx, qx, px_5y, lx, dx, Lx, Tx, ex
  )

# ----------------------------
# 5) Inferir age_interval_width cuando falte
#    - Si falta en una edad intermedia: usar diff con la siguiente edad
#    - Si falta en la última edad: dejar NA (intervalo abierto)
# ----------------------------
df_clean <- df_clean |>
  arrange(
    period_start_year, period_end_year,
    location_id, sex_id, life_table_label, age_start
  ) |>
  group_by(
    period_start_year, period_end_year,
    location_id, sex_id, life_table_label
  ) |>
  mutate(
    next_age_start = lead(age_start),
    age_interval_width = if_else(
      is.na(age_interval_width) & !is.na(next_age_start),
      as.integer(next_age_start - age_start),
      age_interval_width
    ),
    age_interval_open = is.na(age_interval_width)
  ) |>
  ungroup() |>
  select(-next_age_start)

# ----------------------------
# 6) Checks de mapeo
# ----------------------------
if (anyNA(df_clean$sex_id)) {
  bad <- df_clean |>
    filter(is.na(sex_id)) |>
    distinct(sex_source_value) |>
    arrange(sex_source_value)
  stop(
    "Hay sex_source_value sin mapeo a sex_id. Valores:\n",
    paste(bad$sex_source_value, collapse = ", ")
  )
}

if (anyNA(df_clean$location_id)) {
  bad <- df_clean |>
    filter(is.na(location_id)) |>
    distinct(location_name) |>
    arrange(location_name)
  stop(
    "Hay location_name sin mapeo a location_id. Valores:\n",
    paste(bad$location_name, collapse = ", ")
  )
}

# ----------------------------
# 7) QC duro: age_interval_width NA solo puede estar en la última edad del estrato
# ----------------------------
qc_age_open_bad <- df_clean |>
  group_by(period_start_year, period_end_year, location_id, sex_id, life_table_label) |>
  mutate(max_age_start = max(age_start, na.rm = TRUE)) |>
  ungroup() |>
  filter(is.na(age_interval_width) & age_start != max_age_start)

if (nrow(qc_age_open_bad) > 0) {
  dir.create(P$QC_DIR, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    as.data.table(qc_age_open_bad),
    file.path(P$QC_DIR, "qc_bad_open_age_interval_rows.csv")
  )
  stop(
    "Hay age_interval_width = NA en edades que NO son la última del estrato. ",
    "Revisa QC: qc_bad_open_age_interval_rows.csv"
  )
}

# ----------------------------
# 8) Validación por spec
# ----------------------------
spec <- read_spec(spec_path)
dt <- as.data.table(df_clean)

validate_by_spec(dt, spec)

# ----------------------------
# 9) Diccionario
# ----------------------------
dict <- dict_from_spec(
  spec,
  dataset_version = "v1.0.0",
  run_id = NA_character_,
  config_dir = P$CONFIG_DIR
)
dict_ext <- enrich_dict_with_stats(dict, dt)

# ----------------------------
# 10) QC específico
# ----------------------------
qc_list <- qc_life_table_basic(dt)
write_qc_list(qc_list, P$QC_DIR)

# QC adicional: conteo de intervalos abiertos
qc_open_summary <- df_clean |>
  count(
    period_start_year, period_end_year,
    location_id, sex_id, life_table_label,
    age_interval_open,
    name = "n_rows"
  )

data.table::fwrite(
  as.data.table(qc_open_summary),
  file.path(P$QC_DIR, "qc_open_age_interval_summary.csv")
)

# ----------------------------
# 11) Export
# ----------------------------
dir.create(P$FINAL_DIR, recursive = TRUE, showWarnings = FALSE)

out_csv  <- file.path(P$FINAL_DIR, "ref_life_table_mortality.csv")
out_dict <- file.path(P$FINAL_DIR, "ref_life_table_mortality_dictionary_ext.csv")

data.table::fwrite(as.data.table(df_clean), out_csv)
data.table::fwrite(as.data.table(dict_ext), out_dict)

message("✅ Export final: ", out_csv)
message("✅ Diccionario ext: ", out_dict)
message("✅ QC en: ", P$QC_DIR)
message("✅ Columna P(x,x+5) usada: ", px_col)
message("✅ DROP_LIMA_METROPOLITANA = ", DROP_LIMA_METROPOLITANA)
