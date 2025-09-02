library(readxl)
library(dplyr)
library(stringr)
library(readr)
library(fs)

# === INPUT ===
info_file <- "studyInformation/studyInformation_publication.xlsx"
sheet <- 1
tmap_folder <- "tMaps_resliced/Amygdala"
output_file <- "sdm_table.txt"

# === LOAD STUDY INFO ===
info <- read_excel(info_file, sheet = sheet) %>%
  mutate(
    stem_clean = str_trim(`image name stem`),
    sample_size = `sample size amygdala`
  )

# === GET ALL NII FILES ===
tmap_files <- dir_ls(tmap_folder, regexp = "\\.nii$", recurse = FALSE)

# === EXTRACT FILE NAMES WITHOUT EXTENSION ===
studies <- path_file(tmap_files) %>%
  str_remove("\\.nii$")

# === BUILD DATA FRAME ===
sdm_table <- tibble(
  study = studies,
  matched_stem = NA_character_,
  n1 = NA_real_,
  t_thr = 0
)

# === MATCH STUDIES TO STEMS ===
for (i in seq_len(nrow(sdm_table))) {
  s <- sdm_table$study[i]
  
  # Exact match logic:
  match_idx <- which(
    info$stem_clean == s |
      str_starts(s, paste0(info$stem_clean, "_")) |
      (info$stem_clean == "Benzait2023" & str_starts(s, "BenzaitUnpublished"))
  )
  
  if (length(match_idx) == 1) {
    sdm_table$matched_stem[i] <- info$stem_clean[match_idx]
    sdm_table$n1[i] <- info$sample_size[match_idx]
  }
}

# === WARN IF MISSING VALUES ===
unmatched <- filter(sdm_table, is.na(n1))
if (nrow(unmatched) > 0) {
  warning("⚠️ Some image files could not be matched to sample sizes:\n",
          paste(unmatched$study, collapse = "\n"))
}

# === FINALIZE TABLE ===
final_table <- sdm_table %>%
  select(study, n1, t_thr) %>%
  filter(!is.na(n1))

# === WRITE TAB-DELIMITED FILE ===
write_tsv(final_table, output_file)

cat("✔️ Created", output_file, "with", nrow(final_table), "entries.\n")


library(R.utils)
library(fs)

# Path to your SDM-ready folder
root_dir <- "SDMtest/SDMdata"

# Recursively find all .nii files
nii_files <- dir_ls(root_dir, regexp = "\\.nii$", recurse = TRUE)

# Compress each file to .nii.gz (leaving original untouched or deleted if needed)
for (f in nii_files) {
  gz_file <- paste0(f, ".gz")
  if (!file_exists(gz_file)) {
    gzip(f, destname = gz_file, overwrite = TRUE)
    # Optionally: remove the original .nii file
    # file_delete(f)
  }
}

cat("✔️ All .nii files gzipped (retaining originals unless you delete them).\n")
