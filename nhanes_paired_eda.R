# ==============================================================================
# NHANES PAIRED EDA
# Analyzes chemical co-exposure for 10 dataset pairs across cycles 2003-2017
# Outputs: eda_results list (used by nhanes_paired_eda.Rmd)
#
# USAGE:
#   Step 1 — Run this script once in R to download all CDC data:
#              source("nhanes_paired_eda.R")
#            This writes nhanes_eda_cache.rds next to the script (~30 min).
#
#   Step 2 — Knit the report:
#              rmarkdown::render("nhanes_paired_eda.Rmd")
#            The Rmd loads the cache instantly; no downloads happen.
#
#   To force a re-download, delete nhanes_eda_cache.rds and re-run Step 1.
# ==============================================================================

library(nhanesA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(scales)
library(corrplot)
library(gridExtra)

options(width = 180)

CACHE_FILE <- "nhanes_eda_cache.rds"

# If cache exists, load it and skip all downloads ─────────────────────────────
if (file.exists(CACHE_FILE)) {
  message("\n── Cache found: loading ", CACHE_FILE, " ──────────────────────────")
  cache         <- readRDS(CACHE_FILE)
  merged_pairs  <- cache$merged_pairs
  lod_summaries <- cache$lod_summaries
  lab_store     <- cache$lab_store
  demo_store    <- cache$demo_store
  needed        <- cache$needed
  message("   Loaded ", length(merged_pairs), " pairs from cache. Skipping downloads.")
  SKIP_DOWNLOADS <- TRUE
} else {
  SKIP_DOWNLOADS <- FALSE
}

# ==============================================================================
# SECTION 1 — CONFIGURATION
# ==============================================================================

PAIRS <- list(
  list(a = "PBCD_EXT", b = "PFC"),
  list(a = "IHGEM",    b = "PFC"),
  list(a = "PFC",      b = "PHTHTE"),
  list(a = "PFAS",     b = "PHTHTE"),
  list(a = "SSPFAS",   b = "PHTHTE"),
  list(a = "PBCD_EXT", b = "PFAS"),
  list(a = "IHGEM",    b = "PFAS"),
  list(a = "PBCD_EXT", b = "SSPFAS"),
  list(a = "IHGEM",    b = "SSPFAS"),
  list(a = "PFAS",     b = "UPHOPM")
)

# CDC file codes for every dataset × cycle
FILE_CODES <- list(
  PBCD_EXT = c("2011" = "PbCd_G",   "2013" = "PBCD_H",   "2015" = "PBCD_I",   "2017" = "PBCD_J"),
  IHGEM    = c("2007" = "IHg_E",    "2009" = "IHg_F",    "2011" = "IHgEM_G",  "2013" = "IHGEM_H",
               "2015" = "IHGEM_I",  "2017" = "IHGEM_J"),
  PFC      = c("2003" = "L24PFC_C", "2005" = "PFC_D",    "2007" = "PFC_E",    "2009" = "PFC_F",
               "2011" = "PFC_G"),
  PFAS     = c("2013" = "PFAS_H",   "2015" = "PFAS_I",   "2017" = "PFAS_J"),
  SSPFAS   = c("2013" = "SSPFAS_H", "2017" = "SSPFAS_J"),
  PHTHTE   = c("2003" = "L24PH_C",  "2005" = "PHTHTE_D", "2007" = "PHTHTE_E", "2009" = "PHTHTE_F",
               "2011" = "PHTHTE_G", "2013" = "PHTHTE_H", "2015" = "PHTHTE_I", "2017" = "PHTHTE_J"),
  UPHOPM   = c("2011" = "UPHOPM_G", "2013" = "UPHOPM_H", "2015" = "UPHOPM_I")
)

DEMO_CODES <- c(
  "2003" = "DEMO_C", "2005" = "DEMO_D", "2007" = "DEMO_E", "2009" = "DEMO_F",
  "2011" = "DEMO_G", "2013" = "DEMO_H", "2015" = "DEMO_I", "2017" = "DEMO_J"
)

# Demographic variables to retain
DEMO_VARS <- c("SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDEXPRG",
               "INDFMPIR", "SDMVSTRA", "SDMVPSU", "WTMEC2YR")

# ==============================================================================
# SECTION 2 — HELPER FUNCTIONS
# ==============================================================================

# Chemical result columns: LBX* or URX*, excluding LOD comment columns (*LC)
chem_cols <- function(df) {
  nms <- names(df)
  nms[grepl("^(LBX|URX)", nms) & !grepl("LC$", nms)]
}

# LOD comment/flag columns -- NHANES uses several patterns:
#   LBD*LC (most blood/serum)  URD*LC (most urine)
#   LBX*LC (some older files)  URX*LC (some older urine)
lod_flag_cols <- function(df) {
  nms <- names(df)
  nms[grepl("^(LBD|URD|LBX|URX).+LC$", nms)]
}

# Paired LOD flag column -- tries standard swap first (LBX->LBD, URX->URD),
# then falls back to same-prefix+LC used in some older cycles
paired_lod_flag <- function(col, df_names = NULL) {
  standard <- paste0(sub("^LBX", "LBD", sub("^URX", "URD", col)), "LC")
  fallback  <- paste0(col, "LC")
  if (!is.null(df_names)) {
    if (standard %in% df_names) return(standard)
    if (fallback  %in% df_names) return(fallback)
  }
  standard
}

# Safe download — returns NULL with a message instead of stopping
safe_nhanes <- function(file_code, label = file_code) {
  Sys.sleep(0.4)
  df <- tryCatch(nhanes(file_code), error = function(e) {
    message("    [SKIP] ", label, " — ", conditionMessage(e))
    NULL
  })
  if (is.null(df)) {
    message("    [SKIP] ", label, " — CDC returned NULL (file not available)")
  }
  df
}

# Geometric mean (NA/zero-safe)
geo_mean <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  exp(mean(log(x)))
}

# LOD summary table for one dataset (one or many cycles stacked)
build_lod_summary <- function(lab_list, dataset_key) {
  purrr::imap_dfr(lab_list, function(df, yr) {
    if (is.null(df)) return(NULL)
    chems <- chem_cols(df)
    purrr::map_dfr(chems, function(col) {
      flag <- paired_lod_flag(col, names(df))
      n_total     <- sum(!is.na(df[[col]]))
      n_below_lod <- if (flag %in% names(df)) sum(df[[flag]] == 1, na.rm = TRUE) else NA_integer_
      data.frame(
        dataset     = dataset_key,
        year        = as.integer(yr),
        chemical    = col,
        n_total     = n_total,
        n_below_lod = n_below_lod,
        pct_below   = round(n_below_lod / n_total * 100, 1),
        stringsAsFactors = FALSE
      )
    })
  })
}

# ==============================================================================
# SECTION 3 — DOWNLOAD DEMOGRAPHICS
# ==============================================================================

if (!SKIP_DOWNLOADS) {
  
  message("\n── Downloading demographics ──────────────────────────────────────────────")
  demo_store <- list()
  
  for (yr in names(DEMO_CODES)) {
    code <- DEMO_CODES[yr]
    message("  ", code, "  (", yr, ")")
    df <- safe_nhanes(code, paste0("DEMO / ", yr))
    if (!is.null(df)) {
      keep <- intersect(DEMO_VARS, names(df))
      demo_store[[yr]] <- df[, keep, drop = FALSE] %>% mutate(year = as.integer(yr))
    }
  }
  
  message("  Loaded: ", sum(!sapply(demo_store, is.null)), " cycles")
  
  # ==============================================================================
  # SECTION 4 — DOWNLOAD LAB DATASETS
  # ==============================================================================
  
  message("\n── Downloading lab datasets ──────────────────────────────────────────────")
  
  # Only download datasets that appear in at least one pair
  needed <- unique(c(sapply(PAIRS, `[[`, "a"), sapply(PAIRS, `[[`, "b")))
  lab_store <- setNames(vector("list", length(needed)), needed)
  
  for (ds in needed) {
    lab_store[[ds]] <- list()
    for (yr in names(FILE_CODES[[ds]])) {
      code <- FILE_CODES[[ds]][[yr]]
      message("  ", code, "  (", ds, " / ", yr, ")")
      df <- safe_nhanes(code, paste0(ds, " / ", yr))
      if (!is.null(df)) {
        df$year <- as.integer(yr)
        lab_store[[ds]][[yr]] <- df
      }
    }
  }
  
  # ==============================================================================
  # SECTION 5 — MERGE PAIRS WITH DEMOGRAPHICS
  # ==============================================================================
  
  message("\n── Building merged pair datasets ─────────────────────────────────────────")
  
  merged_pairs <- list()
  
  for (pair in PAIRS) {
    a_key      <- pair$a
    b_key      <- pair$b
    pair_label <- paste0(a_key, "_vs_", b_key)
    
    yrs_a      <- names(lab_store[[a_key]])[!sapply(lab_store[[a_key]], is.null)]
    yrs_b      <- names(lab_store[[b_key]])[!sapply(lab_store[[b_key]], is.null)]
    common_yrs <- intersect(yrs_a, yrs_b)
    
    if (length(common_yrs) == 0) {
      message("  [SKIP] ", pair_label, " — no overlapping cycles")
      next
    }
    
    message("  ", pair_label, "  (cycles: ", paste(common_yrs, collapse = ", "), ")")
    
    yr_dfs <- list()
    for (yr in common_yrs) {
      df_a <- lab_store[[a_key]][[yr]]
      df_b <- lab_store[[b_key]][[yr]]
      demo <- demo_store[[yr]]
      if (is.null(df_a) || is.null(df_b) || is.null(demo)) next
      
      # Columns to keep from each lab file: SEQN + result + LOD flag cols
      keep_a <- unique(c("SEQN", chem_cols(df_a), lod_flag_cols(df_a)))
      keep_b <- unique(c("SEQN", chem_cols(df_b), lod_flag_cols(df_b)))
      
      sub_a <- df_a[, intersect(keep_a, names(df_a)), drop = FALSE]
      sub_b <- df_b[, intersect(keep_b, names(df_b)), drop = FALSE]
      
      # Suffix non-SEQN columns to avoid clashes in the merged frame
      rename_nonsqn <- function(df, suffix) {
        idx <- names(df) != "SEQN"
        names(df)[idx] <- paste0(names(df)[idx], ".", suffix)
        df
      }
      sub_a <- rename_nonsqn(sub_a, a_key)
      sub_b <- rename_nonsqn(sub_b, b_key)
      
      merged <- sub_a %>%
        inner_join(sub_b, by = "SEQN") %>%
        inner_join(demo,  by = "SEQN") %>%
        mutate(
          year      = as.integer(yr),
          pair      = pair_label,
          gender    = factor(as.integer(RIAGENDR), levels = 1:2, labels = c("Male", "Female")),
          race_eth  = factor(as.integer(RIDRETH1), levels = 1:5,
                             labels = c("Mexican American", "Other Hispanic",
                                        "Non-Hispanic White", "Non-Hispanic Black", "Other")),
          age_group = cut(RIDAGEYR,
                          breaks = c(0, 18, 45, 65, Inf),
                          labels = c("<18", "18–44", "45–64", "65+"),
                          right  = FALSE)
        )
      
      yr_dfs[[yr]] <- merged
    }
    
    if (length(yr_dfs) > 0) {
      merged_pairs[[pair_label]] <- bind_rows(yr_dfs)
      message("    → n = ", nrow(merged_pairs[[pair_label]]),
              " across ", length(yr_dfs), " cycle(s)")
    }
  }
  
  # ==============================================================================
  # SECTION 6 — LOD SUMMARIES
  # ==============================================================================
  
  message("\n── Building LOD summaries ────────────────────────────────────────────────")
  
  lod_summaries <- setNames(
    lapply(needed, function(ds) build_lod_summary(lab_store[[ds]], ds)),
    needed
  )
  
  # Save cache so the Rmd doesn't need to re-download
  message("\n── Saving cache to ", CACHE_FILE, " ──────────────────────────────────")
  saveRDS(list(
    merged_pairs  = merged_pairs,
    lod_summaries = lod_summaries,
    lab_store     = lab_store,
    demo_store    = demo_store,
    needed        = needed
  ), CACHE_FILE)
  message("   Cache saved (", round(file.size(CACHE_FILE) / 1e6, 1), " MB)")
  
} # end if (!SKIP_DOWNLOADS)

# ==============================================================================
# SECTION 7 — EDA PLOT FUNCTIONS
# ==============================================================================

# Colour palette
PAL_YEAR   <- scale_color_brewer(palette = "Set1", name = "Cycle")
PAL_GENDER <- scale_fill_manual(values  = c("Male" = "#4472C4", "Female" = "#ED7D31"), name = "Sex")
PAL_RACE   <- scale_fill_brewer(palette = "Set2", name = "Race/Ethnicity")

theme_nhanes <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, colour = "grey40"),
          strip.text = element_text(face = "bold"))
}

# 1. Raw + log10 histograms side-by-side for one chemical
plot_dist <- function(df, col, title_prefix = "") {
  vals <- df[[col]][!is.na(df[[col]]) & df[[col]] > 0]
  if (length(vals) < 5) return(NULL)
  
  d <- data.frame(x = vals)
  p1 <- ggplot(d, aes(x = x)) +
    geom_histogram(bins = 60, fill = "#4472C4", alpha = 0.85) +
    labs(title = paste0(title_prefix, col), subtitle = "Raw scale", x = col, y = "Count") +
    theme_nhanes()
  
  p2 <- ggplot(d, aes(x = x)) +
    geom_histogram(bins = 60, fill = "#ED7D31", alpha = 0.85) +
    scale_x_log10(labels = label_comma()) +
    labs(subtitle = "Log10 scale", x = paste0("log10(", col, ")"), y = "Count") +
    theme_nhanes()
  
  grid.arrange(p1, p2, ncol = 2)
}

# 2. Geometric means by cycle for one chemical
plot_gm_by_cycle <- function(df, col, dataset_key) {
  gm_df <- df %>%
    filter(!is.na(.data[[col]]), .data[[col]] > 0) %>%
    group_by(year) %>%
    summarise(gm = geo_mean(.data[[col]]), n = n(), .groups = "drop")
  
  if (nrow(gm_df) == 0) return(NULL)
  
  ggplot(gm_df, aes(x = factor(year), y = gm, group = 1)) +
    geom_line(colour = "#4472C4", linewidth = 1) +
    geom_point(aes(size = n), colour = "#4472C4") +
    geom_text(aes(label = comma(n)), vjust = -1.2, size = 3, colour = "grey40") +
    labs(title = paste0("Geometric mean by cycle: ", col),
         subtitle = paste0("Dataset: ", dataset_key),
         x = "Cycle start year", y = "Geometric mean", size = "N") +
    theme_nhanes()
}

# 3. Geometric means by sex × age group
plot_gm_by_demo <- function(df, col, dataset_key) {
  gm_df <- df %>%
    filter(!is.na(.data[[col]]), .data[[col]] > 0,
           !is.na(gender), !is.na(age_group)) %>%
    group_by(gender, age_group) %>%
    summarise(gm = geo_mean(.data[[col]]), n = n(), .groups = "drop")
  
  if (nrow(gm_df) == 0) return(NULL)
  
  ggplot(gm_df, aes(x = age_group, y = gm, fill = gender)) +
    geom_col(position = "dodge", alpha = 0.9) +
    geom_text(aes(label = round(gm, 2)),
              position = position_dodge(width = 0.9), vjust = -0.4, size = 3) +
    PAL_GENDER +
    labs(title = paste0("Geometric mean by sex & age: ", col),
         subtitle = dataset_key, x = "Age group", y = "Geometric mean") +
    theme_nhanes()
}

# 4. Scatter A vs B (log-log) coloured by cycle
plot_scatter <- function(df, col_a, col_b) {
  d <- df %>%
    filter(!is.na(.data[[col_a]]), !is.na(.data[[col_b]]),
           .data[[col_a]] > 0, .data[[col_b]] > 0) %>%
    mutate(year_f = factor(year))
  
  if (nrow(d) < 5) return(NULL)
  
  ggplot(d, aes(x = .data[[col_a]], y = .data[[col_b]], colour = year_f)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                colour = "black", linewidth = 0.8, linetype = "dashed") +
    scale_x_log10(labels = label_comma()) +
    scale_y_log10(labels = label_comma()) +
    PAL_YEAR +
    labs(title = paste0(col_a, " vs ", col_b),
         x = paste0("log10(", col_a, ")"),
         y = paste0("log10(", col_b, ")")) +
    theme_nhanes()
}

# 5. Correlation heatmap between all A and B chemicals
plot_correlation <- function(df, cols_a, cols_b, pair_label) {
  all_cols <- c(cols_a, cols_b)
  all_cols <- all_cols[all_cols %in% names(df)]
  if (length(all_cols) < 2) return(NULL)
  
  mat <- df %>%
    select(all_of(all_cols)) %>%
    mutate(across(everything(), ~ log1p(pmax(., 0)))) %>%
    cor(use = "pairwise.complete.obs")
  
  # Shorten col names for display
  rownames(mat) <- colnames(mat) <- sub("\\.(PBCD_EXT|IHGEM|PFC|PFAS|SSPFAS|PHTHTE|UPHOPM)$", "", colnames(mat))
  
  corrplot(mat, method = "color", type = "upper", tl.cex = 0.65,
           addCoef.col = "black", number.cex = 0.55,
           title = pair_label, mar = c(0, 0, 2, 0))
}

# 6. LOD summary bar chart
plot_lod <- function(lod_df, dataset_key) {
  d <- lod_df %>% filter(!is.na(pct_below))
  if (nrow(d) == 0) return(NULL)
  
  ggplot(d, aes(x = reorder(chemical, pct_below), y = pct_below, fill = factor(year))) +
    geom_col(position = "dodge", alpha = 0.9) +
    coord_flip() +
    scale_fill_brewer(palette = "Set1", name = "Cycle") +
    labs(title = paste0("% Below LOD: ", dataset_key),
         x = NULL, y = "% below limit of detection") +
    theme_nhanes()
}

# ==============================================================================
# SECTION 8 — ASSEMBLE EDA RESULTS
# ==============================================================================

message("\n── Assembling EDA results ────────────────────────────────────────────────")

eda_results <- list()

for (pair_label in names(merged_pairs)) {
  parts <- str_split(pair_label, "_vs_")[[1]]
  a_key <- parts[1]
  b_key <- parts[2]
  df    <- merged_pairs[[pair_label]]
  
  # Find suffixed chemical columns in the merged frame
  # Result cols have .DATASET suffix; LOD flag cols contain "LC."
  cols_a <- grep(paste0("^LBX|^URX"), sub(paste0("\\.", a_key, "$"), "", names(df)), value = FALSE)
  # Simpler: grep by suffix
  cols_a_all <- grep(paste0("\\.", a_key, "$"), names(df), value = TRUE)
  cols_b_all <- grep(paste0("\\.", b_key, "$"), names(df), value = TRUE)
  
  # Split into result vs LOD flag columns
  cols_a_res  <- cols_a_all[grepl("^LBX|^URX", sub(paste0("\\.", a_key, "$"), "", cols_a_all))]
  cols_a_lod  <- cols_a_all[grepl("LC\\.", cols_a_all) | grepl("^LBD|^URD", sub(paste0("\\.", a_key, "$"), "", cols_a_all))]
  cols_b_res  <- cols_b_all[grepl("^LBX|^URX", sub(paste0("\\.", b_key, "$"), "", cols_b_all))]
  cols_b_lod  <- cols_b_all[grepl("LC\\.", cols_b_all) | grepl("^LBD|^URD", sub(paste0("\\.", b_key, "$"), "", cols_b_all))]
  
  message("  ", pair_label, ": ", length(cols_a_res), " A-chemicals, ",
          length(cols_b_res), " B-chemicals, n=", nrow(df))
  
  eda_results[[pair_label]] <- list(
    pair_label = pair_label,
    a_key      = a_key,
    b_key      = b_key,
    df         = df,
    cols_a     = cols_a_res,
    cols_b     = cols_b_res,
    lod_a      = lod_summaries[[a_key]],
    lod_b      = lod_summaries[[b_key]]
  )
}

message("\nDone. ", length(eda_results), " pairs ready.")
message("Access via: eda_results[[\"PFAS_vs_PHTHTE\"]]$df  etc.")
message("\nAvailable pairs:")
for (p in names(eda_results)) {
  r <- eda_results[[p]]
  message(sprintf("  %-25s  n = %5d  A-chems: %d  B-chems: %d",
                  p, nrow(r$df), length(r$cols_a), length(r$cols_b)))
}


# ==============================================================================
# SECTION 9 — DASHBOARD EXPORT (Power BI / Tableau)
# ==============================================================================
# Run this section after eda_results is built.
# Produces 5 flat CSV files ready to connect directly to Power BI or Tableau.
# All files are written to the working directory.
# ==============================================================================

message("\n── Exporting dashboard CSVs ──────────────────────────────────────────────")

DATASET_LABELS <- c(
  PBCD_EXT   = "Cadmium, Lead, Hg, Se & Mn (Blood)",
  IHGEM      = "Mercury: Inorganic, Ethyl & Methyl (Blood)",
  PFC        = "Polyfluoroalkyl Chemicals (Serum)",
  PFAS       = "PFAS (Serum)",
  SSPFAS     = "PFAS Surplus (Serum)",
  PHTHTE     = "Phthalates & Plasticizers (Urine)",
  UPHOPM     = "Pyrethroids, Herbicides & OP Metabolites (Urine)"
)

# ── 1. Long-format chemical values with demographics ──────────────────────────
# One row per person × chemical × pair × cycle
message("  Building dash_chemicals_long.csv ...")

chem_long_rows <- list()

for (pair_label in names(eda_results)) {
  r     <- eda_results[[pair_label]]
  df    <- r$df
  a_key <- r$a_key
  b_key <- r$b_key
  
  emit_rows <- function(cols, ds_key) {
    for (col in cols) {
      # Strip the .DATASET suffix to get the raw CDC variable name
      raw_col <- sub(paste0("\\.", ds_key, "$"), "", col)
      # LOD flag col (LBD*LC or URD*LC)
      flag_col      <- paired_lod_flag(raw_col, names(df))
      flag_col_full <- paste0(flag_col, ".", ds_key)
      
      vals <- df[[col]]
      flags <- if (flag_col_full %in% names(df)) df[[flag_col_full]] else NA_integer_
      
      chem_long_rows[[length(chem_long_rows) + 1]] <<- data.frame(
        pair          = pair_label,
        dataset_a     = a_key,
        dataset_b     = b_key,
        dataset       = ds_key,
        dataset_label = unname(DATASET_LABELS[ds_key]),
        chemical      = raw_col,
        year          = df$year,
        SEQN          = df$SEQN,
        value         = vals,
        value_log10   = ifelse(!is.na(vals) & vals > 0, log10(vals), NA_real_),
        below_lod     = as.integer(flags == 1),
        gender        = as.character(df$gender),
        age           = df$RIDAGEYR,
        age_group     = as.character(df$age_group),
        race_eth      = as.character(df$race_eth),
        pir           = df$INDFMPIR,  # poverty-income ratio
        stringsAsFactors = FALSE
      )
    }
  }
  
  emit_rows(r$cols_a, a_key)
  emit_rows(r$cols_b, b_key)
}

dash_chemicals_long <- bind_rows(chem_long_rows) %>%
  filter(!is.na(value)) %>%
  arrange(pair, dataset, chemical, year)

write.csv(dash_chemicals_long, "dash_chemicals_long.csv", row.names = FALSE)
message("    Saved: dash_chemicals_long.csv  (",
        format(nrow(dash_chemicals_long), big.mark=","), " rows)")

# ── 2. Geometric means by chemical × cycle × sex × age group ─────────────────
message("  Building dash_geomeans.csv ...")

dash_geomeans <- dash_chemicals_long %>%
  filter(!is.na(value), value > 0, is.na(below_lod) | below_lod == 0) %>%
  group_by(pair, dataset, dataset_label, chemical, year, gender, age_group) %>%
  summarise(
    n         = n(),
    geo_mean  = exp(mean(log(value))),
    median    = median(value),
    p25       = quantile(value, 0.25),
    p75       = quantile(value, 0.75),
    .groups   = "drop"
  )

write.csv(dash_geomeans, "dash_geomeans.csv", row.names = FALSE)
message("    Saved: dash_geomeans.csv  (",
        format(nrow(dash_geomeans), big.mark=","), " rows)")

# ── 3. Pairwise correlations per pair × cycle ─────────────────────────────────
message("  Building dash_correlations.csv ...")

cor_rows <- list()

for (pair_label in names(eda_results)) {
  r     <- eda_results[[pair_label]]
  df    <- r$df
  a_key <- r$a_key
  b_key <- r$b_key
  
  for (yr in sort(unique(df$year))) {
    df_yr <- df %>% filter(year == yr)
    for (ca in r$cols_a) {
      for (cb in r$cols_b) {
        va <- df_yr[[ca]]; vb <- df_yr[[cb]]
        ok <- !is.na(va) & !is.na(vb) & va > 0 & vb > 0
        if (sum(ok) < 5) next
        r_val <- tryCatch(
          cor(log10(va[ok]), log10(vb[ok])),
          error = function(e) NA_real_
        )
        cor_rows[[length(cor_rows) + 1]] <- data.frame(
          pair       = pair_label,
          dataset_a  = a_key,
          dataset_b  = b_key,
          chemical_a = sub(paste0("\\.", a_key, "$"), "", ca),
          chemical_b = sub(paste0("\\.", b_key, "$"), "", cb),
          year       = yr,
          n          = sum(ok),
          r          = round(r_val, 4),
          r_squared  = round(r_val^2, 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

dash_correlations <- bind_rows(cor_rows) %>% arrange(pair, year, chemical_a, chemical_b)
write.csv(dash_correlations, "dash_correlations.csv", row.names = FALSE)
message("    Saved: dash_correlations.csv  (",
        format(nrow(dash_correlations), big.mark=","), " rows)")

# ── 4. LOD summary flat table ─────────────────────────────────────────────────
message("  Building dash_lod_summary.csv ...")

dash_lod <- bind_rows(lod_summaries, .id = "dataset") %>%
  mutate(dataset_label = unname(DATASET_LABELS[dataset]),
         lod_flag = dplyr::case_when(
           pct_below > 50 ~ "Majority below LOD (caution)",
           pct_below > 20 ~ "Some below LOD",
           TRUE           ~ "Mostly detected (reliable)"
         )) %>%
  arrange(dataset, chemical, year)

write.csv(dash_lod, "dash_lod_summary.csv", row.names = FALSE)
message("    Saved: dash_lod_summary.csv  (",
        format(nrow(dash_lod), big.mark=","), " rows)")

# ── 5. Pair overlap (N shared participants per pair × cycle) ──────────────────
message("  Building dash_pair_overlap.csv ...")

dash_overlap <- purrr::imap_dfr(eda_results, function(r, lbl) {
  r$df %>%
    group_by(year) %>%
    summarise(n_shared = n(), .groups = "drop") %>%
    mutate(
      pair          = lbl,
      dataset_a     = r$a_key,
      dataset_b     = r$b_key,
      label_a       = unname(DATASET_LABELS[r$a_key]),
      label_b       = unname(DATASET_LABELS[r$b_key])
    )
}) %>% arrange(pair, year)

write.csv(dash_overlap, "dash_pair_overlap.csv", row.names = FALSE)
message("    Saved: dash_pair_overlap.csv  (",
        format(nrow(dash_overlap), big.mark=","), " rows)")

message("\n── Dashboard export complete ─────────────────────────────────────────────")
message("  Files written to: ", getwd())
message("  Connect all 5 CSVs in Power BI / Tableau as separate tables.")
message("  Join on: pair + dataset + chemical + year + SEQN as appropriate.")
message("")
message("  Suggested dashboard pages:")
message("    1. Overview        — dash_pair_overlap.csv (cards + cycle filter)")
message("    2. Trends          — dash_geomeans.csv    (line chart, slicer: chemical/sex/age)")
message("    3. Distributions   — dash_chemicals_long.csv (histogram, box plot)")
message("    4. A vs B          — dash_correlations.csv (heatmap, scatter)")
message("    5. Detection rates — dash_lod_summary.csv  (bar chart, colour by lod_flag)")