library(dplyr)

prepare_groups <- function(df, group_col = "Group") {
  df %>%
    mutate(
      !!group_col := factor(
        .data[[group_col]],
        levels = group_levels
      )
    )
}

# -----------------------------------------------------------------------------
# Cell-line display labels
# -----------------------------------------------------------------------------
# The underlying Cell_Line values stay as the plate codes "1.6" / "1.6.24"
# (stable keys for joins and stats). These helpers provide the rendered
# strings used in figures and tables.

# Figures: italic MYH7 via plotmath expressions. Use with a facet labeller,
# e.g. facet_wrap(~ Cell_Line, labeller = cell_line_fig_labeller).
cell_line_fig_labels <- c(
  "1.6"    = "italic('MYH7')~'R403L'",
  "1.6.24" = "italic('MYH7')~'ISO'"
)
cell_line_fig_labeller <-
  ggplot2::as_labeller(cell_line_fig_labels, ggplot2::label_parsed)

# Tables: plain string with plate code retained in parentheses.
cell_line_table_labels <- c(
  "1.6"    = "MYH7 R403L (1.6)",
  "1.6.24" = "MYH7 ISO (1.6.24)"
)

# Rewrite the Cell_Line column of a data.frame to the table-display strings.
# Call right before writing a CSV / xlsx sheet.
relabel_cell_line_table <- function(df, col = "Cell_Line") {
  if (col %in% names(df)) {
    vals   <- as.character(df[[col]])
    mapped <- cell_line_table_labels[vals]
    mapped[is.na(mapped)] <- vals[is.na(mapped)]   # keep unknown codes as-is
    df[[col]] <- unname(mapped)
  }
  df
}

read_cyteseer_singlecell <- function(root_dir, plate_id, skip = 71) {
  if (!dir.exists(root_dir)) {
    warning("Directory does not exist: ", root_dir)
    return(NULL)
  }
  all_csv <- list.files(path = root_dir, pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
  files <- all_csv[grepl("Primary_DataTable\\.csv$", basename(all_csv),
                         ignore.case = TRUE)]
  if (length(files) == 0) {
    warning("No Primary_DataTable.csv files found in: ", root_dir)
    return(NULL)
  }
  purrr::map_dfr(files, function(f) {
    well <- basename(dirname(f))
    df <- read.csv(f, skip = skip, header = TRUE, fill = TRUE,
                   check.names = TRUE)
    df$Plate       <- plate_id
    df$Well        <- well
    df$source_file <- f
    df
  })
}

load_plate_maps <- function() {
  plate_map_16 <- read.delim("../1_data/Plate_Map_1.6.txt",
                             stringsAsFactors = FALSE) %>%
    mutate(Well      = as.character(Well),
           Condition = as.character(Condition),
           Cell_Line = as.character(Cell_Line),
           Plate     = "1.6")

  plate_map_1624 <- read.delim("../1_data/Plate_Map_1.6.24.txt",
                               stringsAsFactors = FALSE) %>%
    mutate(Well      = as.character(Well),
           Condition = as.character(Condition),
           Cell_Line = as.character(Cell_Line),
           Plate     = "1.6.24")

  bind_rows(
    plate_map_16   %>% select(Plate, Well, Condition, Cell_Line),
    plate_map_1624 %>% select(Plate, Well, Condition, Cell_Line)
  )
}

load_and_merge_cyteseer <- function(cyteseer_dirs, plate_map_all, skip = 71) {
  singlecell_raw <- purrr::pmap_dfr(
    list(cyteseer_dirs$root_dir, cyteseer_dirs$plate_id),
    ~ read_cyteseer_singlecell(..1, ..2, skip = skip)
  )

  cat("  Raw rows loaded:", nrow(singlecell_raw), "\n")

  singlecell <- singlecell_raw %>%
    left_join(plate_map_all, by = c("Plate", "Well")) %>%
    mutate(
      Cell_Line = Plate,
      well_uid  = paste(Plate, Well, sep = "_")
    ) %>%
    group_by(well_uid) %>%
    mutate(
      cell_id  = as.character(row_number()),
      cell_uid = paste(well_uid, cell_id, sep = "_")
    ) %>%
    ungroup()

  if (exists("group_levels", envir = .GlobalEnv)) {
    singlecell$Condition <- factor(singlecell$Condition,
                                   levels = get("group_levels", envir = .GlobalEnv))
  }
  singlecell$Cell_Line <- factor(singlecell$Cell_Line, levels = c("1.6", "1.6.24"))

  cat("  Cells after merge:", nrow(singlecell), "\n")
  singlecell
}
