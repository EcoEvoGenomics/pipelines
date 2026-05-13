suppressMessages(library(tidyverse))
suppressMessages(library(forcats))
suppressMessages(library(gtools))

# -------- Parse command line arguments ----------------------------------------
args <- commandArgs(trailing = TRUE)
if (length(args) != 8) {
  stop("USAGE: Rscript plot_xpehh.R
  <path to list of xpehh/ihs.csv file paths>
  <path to list of xpehh/ihs.cand.csv file paths>
  <path to annotations .gff>
  <path to chromosome name conversion .tsv>
  <base p-value used to call candidate regions (e.g. 0.01)>
  <main plot width in mm>
  <candidate region plots width in mm>
  <plot height per scan in mm>"
  )
}
input_scans <- args[1]
input_cands <- args[2]
input_gff <- args[3]
chr_conversion_table <- args[4]
base_pval <- as.double(args[5])
width_mm <- as.integer(args[6])
cand_mm <- as.integer(args[7])
height_mm <- as.integer(args[8])

scan_files <- readLines(input_scans)
cand_files <- readLines(input_cands)

n_scans <- length(scan_files)
n_cands <- length(cand_files)
stopifnot(n_cands == n_scans)

# -------- Parse chromosomes rename tsv ----------------------------------------
renamed_chrs <- read.table(chr_conversion_table)$V1
names(renamed_chrs) <- read.table(chr_conversion_table)$V2

# -------- Parse input xp-EHH and IHS .CSV files -------------------------------
print("Parsing input files ...")

scan_list <- vector("list", length(scan_files))
cand_list <- vector("list", length(cand_files))

for (scan_index in seq_along(scan_files)) {

  file_path <- scan_files[scan_index]

  file_is_ihs <- "IHS" %in% colnames(read.csv(file_path, nrows = 1))
  file_nrows <- as.integer(system(paste("wc -l <", file_path), intern = TRUE))
  scan_nsnps <- file_nrows - 1

  tmp <- paste(tempfile(), ".csv.tmp", sep = "")
  awk_fields <- ifelse(file_is_ihs, "'{print $1,$2,$4}'", "'{print $1,$2,$8}'")
  paste("cat", file_path, "| awk -F,", awk_fields, ">", tmp) |> system()

  parsed_scanfile <- read.table(tmp, header = TRUE) |>
    as_tibble() |>
    drop_na(LOGPVALUE) |>
    mutate(SCAN = file_path) |>
    mutate(SCAN_BONFERRONI_THRESHOLD = -log10(base_pval / scan_nsnps))

  scan_list[[scan_index]] <- parsed_scanfile
  rm(parsed_scanfile)
  invisible(gc())

}

for (cand_index in seq_along(cand_files)) {

  file_path <- cand_files[cand_index]
  scan_index <- cand_index # Must be ordered identically in input
  matched_scan_path <- scan_files[scan_index]

  parsed_candfile <- read.csv(file_path, header = TRUE) |>
    as_tibble() |>
    select(CHR, START, END, MAX_MRK) |>
    mutate(SCAN = matched_scan_path)

  cand_list[[cand_index]] <- parsed_candfile
  rm(parsed_candfile)
  invisible(gc())

}

parsed_scans <- bind_rows(scan_list)
parsed_cands <- bind_rows(cand_list)

format_scans_cands <- function(scans_or_cands, renamed_chrs) {
  scans_or_cands <- scans_or_cands |>
    mutate(SCAN = factor(SCAN, levels = unique(SCAN))) |>
    mutate(CHR = factor(CHR, levels = gtools::mixedsort(unique(CHR)))) |>
    mutate(across(CHR, \(x) fct_recode(x, !!!renamed_chrs)))
}

parsed_scans <- format_scans_cands(parsed_scans, renamed_chrs)
parsed_cands <- format_scans_cands(parsed_cands, renamed_chrs)

print("Input files parsed")

# -------- Parse annotations from GFF file -------------------------------------
print("Parsing annotations from GFF ...")

annots <- input_gff |>
  read_tsv(
    comment = "#",
    show_col_types = FALSE,
    col_names = c(
      "CHR", "SOURCE", "FEATURE",
      "START", "END", "SCORE", "STRAND", "FRAME",
      "ATTRIBUTE"
    )
  ) |>
  filter(FEATURE == "gene") |>
  mutate(CHR = factor(CHR)) |>
  mutate(across(CHR, \(x) fct_recode(x, !!!renamed_chrs)))

print("GFF parsed")

# -------- Plotting constants --------------------------------------------------
textsize <- 6
constant_y_max <- max(parsed_scans$LOGPVALUE)
constant_y_min <- -(constant_y_max / 6)
y_gap <- constant_y_min / 3
cand_padding <- 1e5

scan_names <- unique(parsed_scans$SCAN)
scan_labels <- tools::file_path_sans_ext((basename(scan_files)))
scan_labels <- tools::file_path_sans_ext(scan_labels)
scan_labels_for_paths <- scan_labels
scan_labels <- str_replace_all(scan_labels, "_", " - ")
names(scan_labels_for_paths) <- scan_names
names(scan_labels) <- scan_names
chr_names <- levels(parsed_scans$CHR)
chr_labels <- chr_names
for (chr_index in seq_along(chr_labels)) {
  chr <- chr_labels[chr_index]
  chr_size <- max(parsed_scans$POSITION[parsed_scans$CHR == chr])
  if (nchar(chr) > (((chr_size / 1e9) * width_mm / 1.75))) {
    chr_labels[chr_index] <- ""
  }
}
names(chr_labels) <- chr_names
facet_labels <- c(chr_labels, scan_labels)

axislabels_common <- list(
  ylab(expression(bold(bolditalic("-log")[10] ~ "P")))
)

theme_common <- theme(
  axis.line = element_line(colour = "black", linewidth = 0.1),
  axis.text = element_text(colour = "black", size = textsize),
  axis.ticks.length = unit(0.5, "mm"),
  axis.ticks = element_line(colour = "black", linewidth = 0.1),
  axis.title.y = element_text(
    angle = 0,
    hjust = 0,
    vjust = 1,
    colour = "black",
    size = textsize,
    face = "bold"
  ),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.spacing.x = unit(0, "mm"),
  panel.spacing.y = unit(5, "mm"),
  strip.background = element_blank(),
  strip.clip = "off",
  strip.text.x = element_text(
    colour = "black",
    size = textsize,
    margin = margin(0, 0, 0, 0, unit = "mm")
  ),
  strip.text.y.left = element_text(
    angle = 0,
    colour = "black",
    hjust = 0,
    vjust = 1,
    margin = margin(
      l = 0.75,
      r = -max(nchar(scan_labels)) / (6 / textsize),
      unit = "mm"
    ),
    size = textsize
  )
)

# -------- Plot main Manhattan plot --------------------------------------------
print("Plotting main plot ...")

manhattan <- ggplot(parsed_scans, aes(x = POSITION, y = LOGPVALUE)) +
  coord_cartesian(clip = "off") +
  facet_grid(
    cols = vars(CHR),
    rows = vars(SCAN),
    labeller = as_labeller(facet_labels),
    scales = "free_x",
    space = "free_x",
    switch = "both"
  ) +
  geom_point(
    aes(colour = as.integer(CHR) %% 2 == 0),
    size = 0.25, stroke = 0, show.legend = FALSE
  ) +
  geom_point(
    aes(alpha = LOGPVALUE > SCAN_BONFERRONI_THRESHOLD),
    colour = "red", size = 0.25, stroke = 0, show.legend = FALSE
  ) +
  geom_rect(
    data = parsed_cands,
    aes(
      xmin = START,
      xmax = END,
      ymax = y_gap,
      ymin = y_gap + (constant_y_min - y_gap)
    ),
    inherit.aes = FALSE,
    fill = "red"
  ) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(
    guide = guide_axis(cap = TRUE),
    breaks = seq(
      from = 0,
      to = ceiling(constant_y_max),
      length.out = ceiling(height_mm / 7.5)
    ),
    labels = seq(
      from = 0,
      to = ceiling(constant_y_max),
      length.out = ceiling(height_mm / 7.5)
    ) |> round(digits = 0),
    limits = c(constant_y_min, ceiling(constant_y_max))
  ) +
  scale_alpha_manual(values = c(0, 1)) +
  scale_colour_manual(values = c("black", "grey60")) +
  axislabels_common +
  theme_common +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

outname <- paste(tools::file_path_sans_ext(input_scans), "_main.png", sep = "")
ggsave(
  manhattan,
  filename = outname,
  units = "mm",
  dpi = 1600,
  width = width_mm,
  height = height_mm * length(unique(parsed_scans$SCAN))
)

print(paste("Main plot saved to ", outname, sep = ""))

# -------- Plot candidate region Manhattan plots -------------------------------
print("Plotting candidate region plots and writing candidate genes to file ...")

cand_outputdir <- tools::file_path_sans_ext(input_scans)
dir.create(cand_outputdir)

# First merge overlapping candidate intervals into combined plotting windows
merged_regions <- parsed_cands |>
  arrange(CHR, START, END) |>
  group_by(CHR) |>
  group_modify(\(cands, key) {

    merged <- list()
    current_start <- cands$START[1]
    current_end   <- cands$END[1]

    for (i in seq_len(nrow(cands))) {
      start <- cands$START[i]
      end <- cands$END[i]

      if (start <= current_end) {
        current_end <- max(current_end, end)
      } else {
        merged[[length(merged) + 1]] <- tibble(
          START = current_start,
          END = current_end
        )
        current_start <- start
        current_end <- end
      }
    }

    merged[[length(merged) + 1]] <- tibble(
      START = current_start,
      END = current_end
    )
    bind_rows(merged)

  }) |>
  ungroup() |>
  arrange(CHR, START, END)

for (region_index in seq_len(nrow(merged_regions))) {

  chr <- as.character(merged_regions$CHR[region_index])
  chr_min <- 0
  chr_max <- max(parsed_scans$POSITION[parsed_scans$CHR == chr])

  region_start <- merged_regions$START[region_index]
  region_end <- merged_regions$END[region_index]
  clamped_window_start <- max(chr_min, region_start - cand_padding)
  clamped_window_end <- min(chr_max, region_end + cand_padding)

  chr_outputdir <- paste(cand_outputdir, "/", chr, sep = "")
  dir.create(chr_outputdir)

  region_scans <- parsed_scans |>
    filter(CHR == chr) |>
    filter(POSITION >= clamped_window_start) |>
    filter(POSITION <= clamped_window_end) |>
    mutate(SCAN = factor(SCAN)) |>
    droplevels()

  region_cands <- parsed_cands |>
    filter(CHR == chr) |>
    filter(START >= clamped_window_start) |>
    mutate(END = ifelse(END > chr_max, chr_max, END)) |>
    filter(END <= clamped_window_end) |>
    mutate(SCAN = factor(SCAN, levels = levels(region_scans$SCAN))) |>
    droplevels()

  region_annots <- annots |>
    filter(CHR == chr) |>
    filter(START >= clamped_window_start) |>
    filter(END <= clamped_window_end)

  # Identify candidate genes and output them only if there are annotations
  if (nrow(region_annots) != 0) {

    region_annots <- region_annots |>
      crossing(SCAN = levels(region_scans$SCAN)) |>
      mutate(SCAN = factor(SCAN, levels = levels(region_scans$SCAN))) |>
      rowwise() |>
      mutate(
        IS_CANDGENE = {
          this_scan <- SCAN
          this_start <- START
          this_end <- END
          scan_cands <- region_cands |>
            filter(as.character(SCAN) == as.character(this_scan))
          no_scan_cands <- nrow(scan_cands) == 0
          if (no_scan_cands) {
            FALSE
          } else {
            any(
              this_start <= scan_cands$END &
                this_end >= scan_cands$START
            )
          }
        }
      ) |>
      ungroup()

    # Output gff file with candidate genes before modifying annotations for plot
    for (scan in levels(region_scans$SCAN)) {

      candgenes_in_scan <- region_annots |>
        filter(SCAN == scan, IS_CANDGENE) |>
        select(-SCAN, -IS_CANDGENE)

      no_candgenes <- nrow(candgenes_in_scan) == 0
      if (no_candgenes) {
        next
      }

      candgenes_outname <- paste(
        chr_outputdir, "/", chr,
        "_", format(region_start, scientific = FALSE),
        "_", format(region_end, scientific = FALSE),
        "_", scan_labels_for_paths[which(names(scan_labels_for_paths) == scan)],
        ".gff",
        sep = ""
      )

      candgenes_in_scan |> write_tsv(file = candgenes_outname)
    }

    # Reverse start and end for genes on minus strand for plotting arrows
    for (index in seq_len(nrow(region_annots))) {
      start <- region_annots$START[index]
      end <- region_annots$END[index]
      if (region_annots$STRAND[index] == "-") {
        region_annots$END[index] <- start
        region_annots$START[index] <- end
      }
    }

  }

  region_ymax <- ceiling(max(region_scans$LOGPVALUE))
  region_ymax_padded <- region_ymax * 1.15
  region_yratio <- region_ymax / constant_y_max

  candhattan <- ggplot(region_scans, aes(x = POSITION, y = LOGPVALUE)) +
    coord_cartesian(clip = "off") +
    facet_grid(
      rows = vars(SCAN),
      labeller = as_labeller(facet_labels),
      switch = "y"
    ) +
    geom_rect(
      data = region_cands,
      aes(
        xmin = START,
        xmax = END,
        ymax = region_ymax,
        ymin = 0
      ),
      inherit.aes = FALSE,
      fill = "red", colour = "red",
      linewidth = 0.1, lty = 2,
      alpha = 0.05
    ) +
    geom_point(colour = "black", size = 0.5, stroke = 0, show.legend = FALSE) +
    geom_point(
      aes(alpha = LOGPVALUE > SCAN_BONFERRONI_THRESHOLD),
      colour = "red", size = 0.5, stroke = 0, show.legend = FALSE
    ) +
    xlab(chr) +
    scale_alpha_manual(values = c(0, 1)) +
    scale_x_continuous(
      breaks = seq(
        from = clamped_window_start,
        to = clamped_window_end,
        length.out = max(2, ceiling(width_mm / 40))
      ),
      limits = c(
        clamped_window_start,
        clamped_window_end
      ),
      expand = expansion(mult = c(0, 0.05)),
      guide = guide_axis(cap = "both")
    ) +
    scale_y_continuous(
      guide = guide_axis(cap = TRUE),
      breaks = seq(
        from = 0,
        to = region_ymax_padded,
        length.out = ceiling(height_mm / 7.5)
      ),
      labels = seq(
        from = 0,
        to = region_ymax_padded,
        length.out = ceiling(height_mm / 7.5)
      ) |> round(digits = 0),
      limits = c(
        constant_y_min * region_yratio,
        region_ymax_padded
      )
    ) +
    axislabels_common +
    theme_common +
    theme(
      axis.title.x = element_text(
        colour = "black",
        face = "bold",
        size = textsize,
        hjust = 0.48
      )
    )

  # Add annotation arrows if they exist in the window
  if (nrow(region_annots) != 0) {
    candhattan <- candhattan +
      geom_segment(
        data = region_annots,
        aes(
          x = START,
          xend = END,
          y = (constant_y_min - (constant_y_min - y_gap) / 2) * region_yratio,
          yend = (constant_y_min - (constant_y_min - y_gap) / 2) * region_yratio,
          colour = IS_CANDGENE
        ),
        arrow = arrow(
          angle = 30,
          length = unit(0.75, "mm"),
          type = "closed"
        ),
        position = position_jitter(
          width = 0,
          height = abs((constant_y_min - y_gap) / 3) * region_yratio
        ),
        linewidth = 0.25,
        lineend = "square",
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      scale_colour_manual(values = c("#bbbbff", "blue"))
  }

  candplot_outname <- paste(
    chr_outputdir, "/",
    chr, "_", format(region_start, scientific = FALSE),
    "_", format(region_end, scientific = FALSE), ".png",
    sep = ""
  )
  ggsave(
    candhattan,
    filename = candplot_outname,
    dpi = 1600,
    units = "mm",
    width = cand_mm,
    height = height_mm * length(unique(parsed_scans$SCAN))
  )
  print(
    paste(
      "Candidate region plot ",
      region_index, "/", nrow(merged_regions),
      " saved to ", candplot_outname,
      " and candidate genes extracted from your .gff",
      sep = ""
    )
  )

}

print("Done")
