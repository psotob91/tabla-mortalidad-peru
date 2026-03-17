# R/io_utils.R
library(here)

paths_mort <- function() {
  list(
    RAW_DIR     = here("data", "raw", "life_table_mortality"),
    QC_DIR      = here("data", "derived", "qc", "life_table_mortality"),
    FINAL_DIR   = here("data", "final", "life_table_mortality"),
    CONFIG_DIR  = here("config")
  )
}