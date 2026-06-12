# ============================================================
# PTB-XL: Read all WFDB ECG records and save one RDS file
# Output format: list(patient_id = data.frame(time, lead, mV))
# Works with records100 (100 Hz, _lr files)
# ============================================================
#
# Required packages (install once):
# install.packages(c("readr", "dplyr", "tidyr"))

library(readr)
library(dplyr)
library(tidyr)

# ---- 0. Configuration ---------------------------------------

BASE_PATH <- normalizePath(
  "C:/Users/gbp19/Documents/differential_entropy_estimation-master/dados_ecg",
  winslash = "/",
  mustWork = TRUE
)

OUTPUT_RDS <- file.path(dirname(BASE_PATH), "ecgs_por_paciente.rds")


# ---- 1. Parse the .hea header file --------------------------

read_hea <- function(hea_path) {
  lines <- readLines(hea_path)

  top <- strsplit(trimws(lines[1]), "\\s+")[[1]]
  n_leads <- as.integer(top[2])
  fs <- as.numeric(top[3])
  n_samp <- as.integer(top[4])

  leads <- lapply(lines[2:(n_leads + 1)], function(line) {
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    gain_str <- parts[3]
    gain <- as.numeric(sub("\\(.*", "", gain_str))
    base <- as.numeric(gsub(".*\\(|\\).*", "", gain_str))

    list(
      file = parts[1],
      fmt = as.integer(parts[2]),
      gain = gain,
      baseline = base,
      lead_name = parts[9]
    )
  })

  list(
    n_leads = n_leads,
    fs = fs,
    n_samp = n_samp,
    leads = leads
  )
}


# ---- 2. Read the binary .dat file ---------------------------
# Format 16 = 16-bit signed integers, little-endian, interleaved by sample.

read_dat <- function(dat_path, hea) {
  raw_con <- file(dat_path, "rb")
  on.exit(close(raw_con), add = TRUE)

  raw_int <- readBin(
    raw_con,
    integer(),
    n = hea$n_samp * hea$n_leads,
    size = 2,
    signed = TRUE,
    endian = "little"
  )

  mat <- matrix(raw_int, nrow = hea$n_samp, ncol = hea$n_leads, byrow = TRUE)

  for (i in seq_len(hea$n_leads)) {
    gain <- hea$leads[[i]]$gain
    baseline <- hea$leads[[i]]$baseline
    mat[, i] <- (mat[, i] - baseline) / gain
  }

  colnames(mat) <- vapply(hea$leads, `[[`, character(1), "lead_name")
  mat
}


# ---- 3. Read one ECG record as a long data frame -------------

read_ecg_record <- function(record_row) {
  record_rel <- record_row$filename_lr
  hea_path <- file.path(BASE_PATH, paste0(record_rel, ".hea"))
  dat_path <- file.path(BASE_PATH, paste0(record_rel, ".dat"))

  if (!file.exists(hea_path)) {
    stop("Header file not found: ", hea_path)
  }

  if (!file.exists(dat_path)) {
    stop("Data file not found: ", dat_path)
  }

  hea <- read_hea(hea_path)
  signal <- read_dat(dat_path, hea)
  time_sec <- seq(0, by = 1 / hea$fs, length.out = nrow(signal))
  lead_levels <- colnames(signal)

  as.data.frame(signal) %>%
    mutate(time = time_sec) %>%
    pivot_longer(-time, names_to = "lead", values_to = "mV") %>%
    mutate(lead = factor(lead, levels = lead_levels)) %>%
    select(time, lead, mV)
}

patient_key <- function(record_row) {
  patient_id <- suppressWarnings(as.integer(record_row$patient_id))

  if (is.na(patient_id)) {
    return(paste0("unknown_ecg_", record_row$ecg_id))
  }

  as.character(patient_id)
}


# ---- 4. Read all records and save one RDS --------------------

db <- read_csv(file.path(BASE_PATH, "ptbxl_database.csv"), show_col_types = FALSE)

records_to_read <- db %>%
  filter(!is.na(filename_lr), nzchar(filename_lr))

total_records <- nrow(records_to_read)
ecgs_por_paciente <- list()
read_errors <- list()

cat("Reading", total_records, "ECG records from", BASE_PATH, "\n")

for (i in seq_len(total_records)) {
  record_row <- records_to_read[i, , drop = FALSE]
  key <- patient_key(record_row)

  ecg_df <- tryCatch(
    read_ecg_record(record_row),
    error = function(e) {
      read_errors[[length(read_errors) + 1]] <<- tibble(
        ecg_id = record_row$ecg_id,
        patient_id = record_row$patient_id,
        filename_lr = record_row$filename_lr,
        error = conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(ecg_df)) {
    if (is.null(ecgs_por_paciente[[key]])) {
      ecgs_por_paciente[[key]] <- ecg_df
    } else {
      ecgs_por_paciente[[key]] <- bind_rows(ecgs_por_paciente[[key]], ecg_df)
    }
  }

  if (i %% 250 == 0 || i == total_records) {
    cat(sprintf("Read %d/%d records\n", i, total_records))
  }
}

errors_df <- if (length(read_errors) == 0) {
  tibble()
} else {
  bind_rows(read_errors)
}

attr(ecgs_por_paciente, "base_path") <- BASE_PATH
attr(ecgs_por_paciente, "generated_at") <- Sys.time()
attr(ecgs_por_paciente, "n_patients") <- length(ecgs_por_paciente)
attr(ecgs_por_paciente, "n_records") <- total_records - nrow(errors_df)
attr(ecgs_por_paciente, "read_errors") <- errors_df

saveRDS(ecgs_por_paciente, OUTPUT_RDS)

cat("\nRDS saved to:", OUTPUT_RDS, "\n")
cat("Patients:", length(ecgs_por_paciente), "\n")
cat("Records :", total_records - nrow(errors_df), "\n")

if (nrow(errors_df) > 0) {
  cat("Records with errors:", nrow(errors_df), "\n")
  print(errors_df)
}

cat("\nExample after loading:\n")
cat("dados <- readRDS(\"ecgs_por_paciente.rds\")\n")
cat("dados[[\"19005\"]]\n")
