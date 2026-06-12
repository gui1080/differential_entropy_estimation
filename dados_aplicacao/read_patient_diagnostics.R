# ============================================================
#  PTB-XL: Build merged dataframe
#  Joins patient/record metadata with SCP diagnostic descriptions
# ============================================================
#
# Required packages (install once):
# install.packages(c("readr", "dplyr", "jsonlite", "tidyr"))

library(readr)
library(dplyr)
library(jsonlite)
library(tidyr)

# ── 0. Configuration ─────────────────────────────────────────
BASE_PATH <- "C:/Users/gbp19/Documents/differential_entropy_estimation-master/dados_ecg"


# ── 1. Load source files ─────────────────────────────────────
db <- read_csv(
  file.path(BASE_PATH, "ptbxl_database.csv"),
  show_col_types = FALSE
)

scp_ref <- read_csv(
  file.path(BASE_PATH, "scp_statements.csv"),
  show_col_types = FALSE
) %>%
  rename(scp_code = 1)   # unnamed first column → rename to "scp_code"


# ── 2. Select only the needed columns from the database ──────
db_slim <- db %>%
  select(ecg_id, patient_id, age, sex, nurse, device, report, scp_codes)


# ── 3. Expand scp_codes: one row per (record × SCP code) ─────
#
# scp_codes is stored as a Python dict string, e.g.:
#   {'NORM': 100.0, 'SR': 0.0}
# Steps:
#   a) convert single → double quotes  (Python → JSON)
#   b) parse JSON to a named list
#   c) convert to a tidy tibble with scp_code + confidence columns

parse_scp <- function(scp_string) {
  tryCatch({
    json    <- gsub("'", '"', scp_string)
    parsed  <- fromJSON(json)
    tibble(
      scp_code   = names(parsed),
      confidence = unlist(parsed)
    )
  }, error = function(e) {
    tibble(scp_code = NA_character_, confidence = NA_real_)
  })
}

db_expanded <- db_slim %>%
  mutate(scp_parsed = lapply(scp_codes, parse_scp)) %>%
  select(-scp_codes) %>%
  unnest(scp_parsed)


# ── 4. Join with scp_statements to get descriptions ──────────
scp_ref_slim <- scp_ref %>%
  select(
    scp_code,
    description,
    statement_category        = `Statement Category`,
    scp_ecg_statement_desc    = `SCP-ECG Statement Description`
  )

final_df <- db_expanded %>%
  left_join(scp_ref_slim, by = "scp_code")


# ── 5. Inspect result ─────────────────────────────────────────
cat("=== Final dataframe: dimensions ===\n")
cat("Rows:", nrow(final_df), " | Cols:", ncol(final_df), "\n\n")

cat("=== Column names ===\n")
print(colnames(final_df))
cat("\n")

cat("=== First 10 rows ===\n")
print(head(final_df, 10))

output_dir <- file.path(dirname(BASE_PATH), "resultados_ecg")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_path <- file.path(output_dir, "patient_diagnostics.csv")
write_csv(final_df, output_path)

cat("\nCSV saved to:", output_path, "\n")
