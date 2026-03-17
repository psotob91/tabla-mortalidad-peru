# scripts/00_print_tree.R

library(here)
library(fs)
library(stringr)
library(purrr)

project_root <- here::here()

# -----------------------------
# Config
# -----------------------------
max_depth <- 6  # ajusta si quieres más/menos profundidad

exclude_dirs <- c(
  ".git", ".Rproj.user", "renv", ".quarto", ".github",
  "node_modules", "__pycache__", ".idea",
  "data/derived", "data/staging", "data/tmp", "data/cache",
  "datos/derivados", "datos/staging", "datos/tmp", "datos/cache"
)

exclude_files_regex <- paste0(
  "(",
  paste(
    c(
      "\\.DS_Store$", "\\.Rhistory$", "\\.RData$",
      "\\.log$", "\\.tmp$", "\\.bak$",
      "\\.zip$", "\\.7z$", "\\.tar$", "\\.gz$"
    ),
    collapse = "|"
  ),
  ")"
)

# -----------------------------
# Helpers
# -----------------------------
norm_sep <- function(p) gsub("\\\\", "/", p)

is_excluded_dir <- function(p) {
  rp <- norm_sep(path_rel(p, project_root))
  parts <- unlist(str_split(rp, "/"))
  any(parts %in% exclude_dirs)
}

is_excluded_file <- function(p) {
  grepl(exclude_files_regex, path_file(p), ignore.case = TRUE)
}

depth_of <- function(p) {
  rp <- norm_sep(path_rel(p, project_root))
  if (rp == "." || rp == "") return(0L)
  length(unlist(str_split(rp, "/")))
}

# Construye árbol tipo "├─" "└─"
make_tree_lines <- function(paths, root) {
  rp <- norm_sep(path_rel(paths, root))
  rp <- rp[rp != "."]
  rp <- sort(unique(rp))
  
  # Solo hasta max_depth
  rp <- rp[map_int(rp, \(x) length(str_split(x, "/")[[1]])) <= max_depth]
  
  # Convertir a “nodos” por nivel
  nodes <- str_split(rp, "/")
  
  # Armar set de prefixes (para saber si un nodo tiene hijos)
  all_prefixes <- unique(unlist(map(nodes, \(x) {
    map_chr(seq_along(x), \(i) paste(x[1:i], collapse = "/"))
  })))
  
  # Para cada prefix, decidir si es folder o file (file si coincide exacto con rp y no tiene hijos)
  has_children <- function(prefix) {
    any(startsWith(rp, paste0(prefix, "/")))
  }
  
  # Listar nodes finales (folders + files)
  # Mantener folders necesarios para estructura aunque no estén en rp (por max_depth ya filtrado)
  items <- unique(c(all_prefixes, rp))
  items <- sort(items)
  
  lines <- map_chr(items, function(item) {
    parts <- str_split(item, "/")[[1]]
    d <- length(parts)
    name <- tail(parts, 1)
    
    indent <- if (d > 1) paste(rep("  ", d - 1), collapse = "") else ""
    
    # Marcar carpetas con "/" (si tiene hijos)
    if (has_children(item) && !endsWith(name, "/")) {
      paste0(indent, "└─ ", name, "/")
    } else {
      paste0(indent, "└─ ", name)
    }
  })
  
  lines
}

# -----------------------------
# Recolectar rutas
# -----------------------------
all_dirs <- dir_ls(project_root, type = "directory", recurse = TRUE)
all_files <- dir_ls(project_root, type = "file", recurse = TRUE)

dirs <- all_dirs[!map_lgl(all_dirs, is_excluded_dir)]
files <- all_files[!map_lgl(all_files, is_excluded_dir)]
files <- files[!map_lgl(files, is_excluded_file)]

# Para árbol: incluir dirs + files, pero limitar profundidad visual
paths_for_tree <- c(dirs, files)

tree_lines <- c(
  paste0("PROJECT ROOT: ", norm_sep(project_root)),
  "",
  make_tree_lines(paths_for_tree, project_root)
)

cat(paste0(tree_lines, collapse = "\n"))

# -----------------------------
# Guardar a txt para compartir
# -----------------------------
out_txt <- file.path(project_root, "tree_project.txt")
writeLines(tree_lines, out_txt)
message("\n\n✅ Árbol guardado en: ", out_txt)