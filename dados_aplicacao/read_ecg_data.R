# ============================================================
#  PTB-XL: Read a WFDB ECG record + look up its diagnostic
#  Works with records100 (100 Hz, _lr files)
# ============================================================
#
# Required packages (install once):
# install.packages(c("readr", "dplyr", "jsonlite", "ggplot2", "tidyr"))

library(readr)
library(dplyr)
library(jsonlite)
library(ggplot2)
library(tidyr)

# ── 0. Configuration ────────────────────────────────────────
BASE_PATH  <- "c:/Users/cesar/Downloads/ptb-xl-a-large-publicly-available-electrocardiography-dataset-1.0.3/ptb-xl-a-large-publicly-available-electrocardiography-dataset-1.0.3/"
RECORD_REL <- "records100/21000/21837_lr"   # change this to load a different record


# ── 1. Parse the .hea header file ───────────────────────────
read_hea <- function(hea_path) {
  lines <- readLines(hea_path)
  
  # Line 1: record-level info
  top     <- strsplit(trimws(lines[1]), "\\s+")[[1]]
  n_leads <- as.integer(top[2])
  fs      <- as.numeric(top[3])
  n_samp  <- as.integer(top[4])
  
  # Lines 2+: one per lead
  leads <- lapply(lines[2:(n_leads + 1)], function(l) {
    parts <- strsplit(trimws(l), "\\s+")[[1]]
    gain_str <- parts[3]                           # e.g. "1000.0(0)/mV"
    gain  <- as.numeric(sub("\\(.*", "", gain_str))  # 1000.0
    base  <- as.numeric(gsub(".*\\(|\\).*", "", gain_str))  # 0
    list(
      file      = parts[1],
      fmt       = as.integer(parts[2]),
      gain      = gain,
      baseline  = base,
      lead_name = parts[9]
    )
  })
  
  list(
    n_leads = n_leads,
    fs      = fs,
    n_samp  = n_samp,
    leads   = leads
  )
}

hea_path <- file.path(BASE_PATH, paste0(RECORD_REL, ".hea"))
hea      <- read_hea(hea_path)

cat("Record     :", RECORD_REL, "\n")
cat("Leads      :", hea$n_leads, "\n")
cat("Sample rate:", hea$fs, "Hz\n")
cat("Samples    :", hea$n_samp, "( =", hea$n_samp / hea$fs, "seconds )\n")
cat("Lead names :", sapply(hea$leads, `[[`, "lead_name"), "\n\n")


# ── 2. Read the binary .dat file ────────────────────────────
# Format 16 = 16-bit signed integers, little-endian, interleaved by sample
read_dat <- function(dat_path, hea) {
  raw_con <- file(dat_path, "rb")
  raw_int <- readBin(raw_con, integer(), n = hea$n_samp * hea$n_leads,
                     size = 2, signed = TRUE, endian = "little")
  close(raw_con)
  
  # Reshape: rows = samples, cols = leads
  mat <- matrix(raw_int, nrow = hea$n_samp, ncol = hea$n_leads, byrow = TRUE)
  
  # Apply gain + baseline to convert ADC units → mV
  for (i in seq_len(hea$n_leads)) {
    gain     <- hea$leads[[i]]$gain
    baseline <- hea$leads[[i]]$baseline
    mat[, i] <- (mat[, i] - baseline) / gain
  }
  
  colnames(mat) <- sapply(hea$leads, `[[`, "lead_name")
  mat
}

dat_path <- file.path(BASE_PATH, paste0(RECORD_REL, ".dat"))
signal   <- read_dat(dat_path, hea)

cat("Signal matrix shape:", nrow(signal), "samples x", ncol(signal), "leads\n")
cat("Value range (mV)   : [", round(min(signal), 3), ",", round(max(signal), 3), "]\n\n")


# ── 3. Look up the diagnostic from ptbxl_database.csv ───────
db <- read_csv(file.path(BASE_PATH, "ptbxl_database.csv"),
               show_col_types = FALSE)

scp_ref <- read_csv(file.path(BASE_PATH, "scp_statements.csv"),
                    show_col_types = FALSE) %>%
  rename(code = 1)   # the first column has no header → rename it to "code"

# Extract ecg_id from the record name (the number before _lr)
ecg_id <- as.integer(sub(".*/(\\d+)_lr$", "\\1", RECORD_REL))

row <- db %>% filter(ecg_id == !!ecg_id)

if (nrow(row) == 0) stop("Record not found in ptbxl_database.csv")

cat("=== Patient Info ===\n")
cat("ECG ID     :", row$ecg_id, "\n")
cat("Patient ID :", row$patient_id, "\n")
cat("Age        :", row$age, "\n")
cat("Sex        :", ifelse(row$sex == 0, "Female", "Male"), "\n")
cat("Recorded   :", row$recording_date, "\n")
cat("Report     :", row$report, "\n")
cat("Strat fold :", row$strat_fold, "\n\n")

# Parse scp_codes (stored as Python dict string → convert to JSON)
scp_raw  <- row$scp_codes
scp_json <- gsub("'", '"', scp_raw)       # single → double quotes
scp_list <- fromJSON(scp_json)             # named numeric vector

cat("=== SCP Codes (raw) ===\n")
print(scp_list)
cat("\n")

# Join with scp_statements to get human-readable descriptions
scp_df <- tibble(
  code       = names(scp_list),
  confidence = unlist(scp_list)
) %>%
  left_join(scp_ref %>% select(code, description, diagnostic_class, diagnostic_subclass),
            by = "code")

cat("=== Diagnostic Summary ===\n")
print(scp_df %>% select(code, confidence, diagnostic_class, description))


# ── 4. Plot the 12-lead ECG ──────────────────────────────────
time_sec <- seq(0, by = 1 / hea$fs, length.out = nrow(signal))

ecg_long <- as.data.frame(signal) %>%
  mutate(time = time_sec) %>%
  pivot_longer(-time, names_to = "lead", values_to = "mV") %>%
  mutate(lead = factor(lead, levels = colnames(signal)))  # preserve lead order

p <- ggplot(ecg_long, aes(x = time, y = mV)) +
  geom_line(color = "#e63946", linewidth = 0.35) +
  facet_wrap(~lead, ncol = 2, scales = "free_y") +
  labs(
    title    = paste0("PTB-XL ECG — Record ", ecg_id),
    subtitle = paste0("Diagnosis: ", paste(scp_df$code, collapse = ", "),
                      " | Age: ", row$age,
                      " | Sex: ", ifelse(row$sex == 0, "Female", "Male")),
    x = "Time (s)",
    y = "Amplitude (mV)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  )

print(p)
# ggsave("ecg_plot.png", p, width = 10, height = 8, dpi = 150)  # uncomment to save
