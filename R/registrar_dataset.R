#' Registrar dataset y diccionario en disco y catálogo maestro
#'
#' @param data          data.frame / tibble con los datos limpios.
#' @param diccionario   data.frame con metadatos de variables (p.ej. diccionario_inei_ext).
#' @param dataset_id    Identificador corto, único y estable (p.ej. "inei_poblacion_1995_2030").
#' @param dataset_name  Nombre descriptivo del dataset.
#' @param version       Versión del dataset (p.ej. "v1.0").
#' @param module        Nombre del módulo/carpeta bajo datos/limpios (p.ej. "poblacion-proyectada-inei").
#' @param base_dir      Carpeta base para datos (por defecto "datos").
#' @param source        Descripción breve de la fuente (INEI, EsSalud, etc.).
#' @param file_stem     Prefijo de los archivos (por defecto = dataset_id).
#' @param write_parquet Guardar en parquet (TRUE/FALSE).
#' @param write_csv     Guardar en CSV (TRUE/FALSE).
#' @param write_rds     Guardar en RDS (TRUE/FALSE).
#' @param update_catalog Actualizar/crear catalogo_datasets.csv (TRUE/FALSE).
#' @param catalog_path  Ruta al catálogo (por defecto "datos/catalogo_datasets.csv").
#' @param register_metadata_db  Registrar diccionario en BD (TRUE/FALSE).
#' @param con           Conexión DBI (si register_metadata_db = TRUE).
#' @param results_schema Esquema de resultados (por defecto "results").
#' @param metadata_table Nombre de la tabla de metadatos (por defecto "metadata_variable").
#' @param table_name    Nombre de la tabla lógica en la BD a la que pertenecen estas variables (p.ej. "population_result").
#'
#' @return Lista con rutas de archivos y, si aplica, catálogo actualizado.
#' @export
registrar_dataset <- function(
    data,
    diccionario,
    dataset_id,
    dataset_name,
    version       = "v1.0",
    module        = "sin-modulo",
    base_dir      = "datos",
    source        = "",
    file_stem     = dataset_id,
    write_parquet = TRUE,
    write_csv     = TRUE,
    write_rds     = TRUE,
    update_catalog      = TRUE,
    catalog_path        = file.path(base_dir, "catalogo_datasets.csv"),
    register_metadata_db = FALSE,
    con              = NULL,
    results_schema   = "results",
    metadata_table   = "metadata_variable",
    table_name       = NULL
) {
  # Paquetes requeridos
  requireNamespace("readr", quietly = TRUE)
  requireNamespace("arrow", quietly = TRUE)
  requireNamespace("writexl", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  
  # 1. Crear carpeta destino para este módulo
  data_dir <- file.path(base_dir, "limpios", module)
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 2. Rutas de archivos de datos
  path_parquet <- file.path(data_dir, paste0(file_stem, ".parquet"))
  path_csv     <- file.path(data_dir, paste0(file_stem, ".csv"))
  path_rds     <- file.path(data_dir, paste0(file_stem, ".rds"))
  
  # 3. Rutas de archivos de diccionario
  dict_stem    <- paste0(file_stem, "_diccionario_ext")
  dict_rds     <- file.path(data_dir, paste0(dict_stem, ".rds"))
  dict_xlsx    <- file.path(data_dir, paste0(dict_stem, ".xlsx"))
  dict_csv     <- file.path(data_dir, paste0(dict_stem, ".csv"))
  
  # 4. Guardar datos
  if (write_parquet) {
    arrow::write_parquet(data, path_parquet)
  }
  if (write_csv) {
    readr::write_csv(data, path_csv)
  }
  if (write_rds) {
    saveRDS(data, path_rds)
  }
  
  # 5. Guardar diccionario extendido
  saveRDS(diccionario, dict_rds)
  writexl::write_xlsx(diccionario, dict_xlsx)
  readr::write_csv(diccionario, dict_csv)
  
  # 6. Actualizar catálogo maestro de datasets
  catalogo_actualizado <- NULL
  if (update_catalog) {
    if (file.exists(catalog_path)) {
      catalogo_existente <- readr::read_csv(catalog_path, show_col_types = FALSE)
    } else {
      catalogo_existente <- tibble::tibble(
        dataset_id      = character(),
        dataset_name    = character(),
        version         = character(),
        module          = character(),
        data_path_parquet = character(),
        data_path_csv   = character(),
        data_path_rds   = character(),
        dictionary_path = character(),
        source          = character(),
        created_at      = as.POSIXct(character())
      )
    }
    
    nuevo_registro <- tibble::tibble(
      dataset_id        = dataset_id,
      dataset_name      = dataset_name,
      version           = version,
      module            = module,
      data_path_parquet = ifelse(write_parquet, path_parquet, NA_character_),
      data_path_csv     = ifelse(write_csv, path_csv, NA_character_),
      data_path_rds     = ifelse(write_rds, path_rds, NA_character_),
      dictionary_path   = dict_xlsx,
      source            = source,
      created_at        = Sys.time()
    )
    
    catalogo_actualizado <- catalogo_existente %>%
      dplyr::filter(!(dataset_id == nuevo_registro$dataset_id & version == nuevo_registro$version)) %>%
      dplyr::bind_rows(nuevo_registro)
    
    readr::write_csv(catalogo_actualizado, catalog_path)
  }
  
  # 7. Registrar diccionario en BD (results.metadata_variable)
  if (register_metadata_db) {
    if (is.null(con)) {
      stop("Para register_metadata_db = TRUE debes pasar una conexión DBI válida en `con`.")
    }
    if (is.null(table_name)) {
      stop("Para register_metadata_db = TRUE debes indicar `table_name` (p.ej. 'population_result').")
    }
    
    # Crear tabla de metadatos si no existe
    ddl_meta <- sprintf("
      CREATE TABLE IF NOT EXISTS %s.%s (
        table_name        VARCHAR(100) NOT NULL,
        column_name       VARCHAR(100) NOT NULL,
        label             TEXT,
        data_type         VARCHAR(50),
        section           VARCHAR(50),
        is_categorical    BOOLEAN,
        omop_domain_id    VARCHAR(50),
        is_omop_concept_id BOOLEAN,
        example_values    TEXT,
        notes             TEXT,
        PRIMARY KEY (table_name, column_name)
      );
    ", results_schema, metadata_table)
    
    DBI::dbExecute(con, ddl_meta)
    
    # Asegurar columnas mínimas en diccionario
    dicc_db <- diccionario
    
    if (!"variable" %in% names(dicc_db)) {
      stop("El diccionario debe tener al menos una columna llamada 'variable' con los nombres de las columnas.")
    }
    if (!"label" %in% names(dicc_db)) {
      dicc_db$label <- NA_character_
    }
    if (!"tipo" %in% names(dicc_db)) {
      dicc_db$tipo <- NA_character_
    }
    if (!"section" %in% names(dicc_db)) {
      dicc_db$section <- NA_character_
    }
    if (!"is_categorical" %in% names(dicc_db)) {
      dicc_db$is_categorical <- NA
    }
    if (!"omop_domain_id" %in% names(dicc_db)) {
      dicc_db$omop_domain_id <- NA_character_
    }
    if (!"is_omop_concept_id" %in% names(dicc_db)) {
      dicc_db$is_omop_concept_id <- NA
    }
    if (!"example_values" %in% names(dicc_db)) {
      dicc_db$example_values <- NA_character_
    }
    if (!"notes" %in% names(dicc_db)) {
      dicc_db$notes <- NA_character_
    }
    
    dicc_db <- dicc_db %>%
      dplyr::mutate(
        table_name  = table_name,
        column_name = .data$variable
      ) %>%
      dplyr::select(
        table_name,
        column_name,
        label,
        data_type      = .data$tipo,
        section,
        is_categorical,
        omop_domain_id,
        is_omop_concept_id,
        example_values,
        notes
      )
    
    # Borrar registros previos de esa tabla si ya existen
    delete_sql <- sprintf(
      "DELETE FROM %s.%s WHERE table_name = $1;",
      results_schema, metadata_table
    )
    DBI::dbExecute(con, delete_sql, params = list(table_name))
    
    # Insertar nuevos metadatos
    DBI::dbWriteTable(
      con,
      DBI::Id(schema = results_schema, table = metadata_table),
      dicc_db,
      append = TRUE,
      row.names = FALSE
    )
  }
  
  # 8. Devolver lista con info útil
  invisible(list(
    data_dir        = data_dir,
    data_paths      = list(
      parquet = if (write_parquet) path_parquet else NA_character_,
      csv     = if (write_csv) path_csv else NA_character_,
      rds     = if (write_rds) path_rds else NA_character_
    ),
    dict_paths      = list(
      rds  = dict_rds,
      xlsx = dict_xlsx,
      csv  = dict_csv
    ),
    catalog_path    = if (update_catalog) catalog_path else NULL,
    catalog_updated = if (update_catalog) catalogo_actualizado else NULL
  ))
}
