suppressMessages(library(tidyverse))
suppressMessages(library(gtools))

# -------- Parse command line arguments ----------------------------------------
args <- commandArgs(trailing = TRUE)
if (length(args) != 7) {
  stop("USAGE: Rscript plot_xpehh.R
  <path to list of xpehh/ihs.csv file paths>
  <path to list of xpehh/ihs.cand.csv file paths>
  <path to annotations .gff>
  <base p-value used to call candidate regions (e.g. 0.01)>
  <main plot width in mm>
  <candidate region plots width in mm>
  <plot height per scan in mm>"
  )
}
input_scans <- args[1]
input_cands <- args[2]
input_gff <- args[3]
base_pval <- as.double(args[4])
width_mm <- as.integer(args[5])
cand_mm <- as.integer(args[6])
height_mm <- as.integer(args[7])

n_scans <- length(readLines(input_scans))
n_cands <- length(readLines(input_cands))
stopifnot(n_cands == n_scans)

# -------- Parse input xp-EHH and IHS .CSV files -------------------------------
print("Parsing input files ...")

parsed_scans <- data.frame()
parsed_cands <- data.frame()

for (file_path in readLines(input_scans)) {

  file_is_ihs <- "IHS" %in% colnames(read.csv(file_path, nrows = 1))
  file_nrows <- as.integer(system(paste("wc -l <", file_path), intern = TRUE))

  tmp <- paste(tempfile(), ".csv.tmp", sep = "")
  awk_fields <- ifelse(file_is_ihs, "'{print $1,$2,$4}'", "'{print $1,$2,$8}'")
  paste("cat", file_path, "| awk -F,", awk_fields, ">", tmp) |> system()

  parsed_scanfile <- read.table(tmp, header = TRUE) |>
    as_tibble() |>
    drop_na(LOGPVALUE) |>
    mutate(SCAN = file_path) |>
    mutate(SCAN_BONFERRONI_THRESHOLD = -log10(base_pval / file_nrows))

  parsed_scans <- rbind(parsed_scans, parsed_scanfile) # Memory inefficient
  rm(parsed_scanfile)
  invisible(gc())

}

for (cand_index in seq_len(n_cands)) {

  file_path <- readLines(input_cands)[cand_index]
  scan_number <- rep(seq_len(n_scans), n_cands / n_scans)[cand_index]
  cand_number <- ceiling(cand_index / n_scans)
  matched_scan_path <- readLines(input_scans)[scan_number]

  parsed_candfile <- read.csv(file_path, header = TRUE) |>
    as_tibble() |>
    select(CHR, START, END, MAX_MRK) |>
    mutate(SCAN = matched_scan_path) |>
    mutate(CAND = cand_number)

  parsed_cands <- rbind(parsed_cands, parsed_candfile) # Memory inefficient
  rm(parsed_candfile)
  invisible(gc())

}

format_scans_cands <- function(scans_or_cands) {
  scans_or_cands |>
    mutate(SCAN = factor(SCAN, levels = unique(SCAN))) |>
    mutate(CHR = factor(CHR, levels = gtools::mixedsort(unique(CHR))))
}

parsed_scans <- format_scans_cands(parsed_scans)
parsed_cands <- format_scans_cands(parsed_cands)

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
  filter(FEATURE == "gene")

print("GFF parsed")

# -------- Plotting constants --------------------------------------------------
textsize <- 6
constant_y_max <- max(parsed_scans$LOGPVALUE)
constant_y_min <- -(constant_y_max / 6)
y_gap <- constant_y_min / 3
cand_padding <- 1e5

scan_names <- unique(parsed_scans$SCAN)
scan_labels <- tools::file_path_sans_ext((basename(readLines(input_scans))))
scan_labels <- tools::file_path_sans_ext(scan_labels)
scan_labels <- str_replace_all(scan_labels, "_", " - ")
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

scales_common <- list(
  scale_alpha_manual(values = c(0, 1)),
  scale_colour_manual(values = c("black", "grey60"))
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
      group = CAND,
      xmin = START,
      xmax = END,
      ymax = y_gap + (((constant_y_min - y_gap) / CAND) * (CAND - 1)),
      ymin = y_gap + (((constant_y_min - y_gap) / CAND) * CAND)
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
  axislabels_common +
  scales_common +
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

for (cand_index in seq_len(nrow(parsed_cands))) {

  chr <- as.character(parsed_cands$CHR[cand_index])
  focal_scan <- parsed_cands$SCAN[cand_index]
  region_start <- parsed_cands$START[cand_index]
  region_end <- parsed_cands$END[cand_index]

  chr_outputdir <- paste(cand_outputdir, "/", chr, sep = "")
  dir.create(chr_outputdir)

  region_scans <- parsed_scans |>
    filter(CHR == chr) |>
    filter(POSITION > region_start - cand_padding) |>
    filter(POSITION < region_end + cand_padding)

  region_cands <- parsed_cands |>
    filter(CHR == chr) |>
    filter(START > region_start - cand_padding) |>
    filter(END < region_end + cand_padding)

  region_annots <- annots |>
    filter(CHR == chr) |>
    filter(START > region_start - cand_padding) |>
    filter(END < region_end + cand_padding)

  # Output gff file with candidate genes before modifying annotations for plot
  candgenes_outname <- paste(
    chr_outputdir, "/",
    chr, "_", format(region_start, scientific = FALSE),
    "_", format(region_end, scientific = FALSE), ".gff",
    sep = ""
  )
  candgenes <- filter(
    region_annots,
    (START > region_start & START < region_end) |
      (END > region_start & END < region_end)
  )
  candgenes |> write_tsv(file = candgenes_outname)

  # Reverse start and end for genes on minus strand for plotting arrows
  for (index in seq_len(nrow(region_annots))) {
    start <- region_annots$START[index]
    end <- region_annots$END[index]
    if (region_annots$STRAND[index] == "-") {
      region_annots$END[index] <- start
      region_annots$START[index] <- end
    }
  }

  region_ymax <- ceiling(max(region_scans$LOGPVALUE))
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
        ymin = 0,
        ymax = max(region_scans$LOGPVALUE)
      ),
      inherit.aes = FALSE,
      fill = "red", colour = "red",
      linewidth = 0.1, lty = 2,
      alpha = 0.1
    ) +
    geom_point(colour = "black", size = 0.5, stroke = 0, show.legend = FALSE) +
    geom_point(
      aes(alpha = LOGPVALUE > SCAN_BONFERRONI_THRESHOLD),
      colour = "red", size = 0.5, stroke = 0, show.legend = FALSE
    ) +
    geom_segment(
      data = region_annots,
      aes(
        x = START,
        xend = END,
        y = (constant_y_min - ((constant_y_min - y_gap)) / 2) * region_yratio,
        yend = (constant_y_min - ((constant_y_min - y_gap)) / 2) * region_yratio
      ),
      arrow = arrow(
        angle = 30,
        length = unit(1, "mm"),
        type = "closed"
      ),
      position = position_jitter(
        width = 0,
        height = abs((constant_y_min - y_gap) / 3) * region_yratio
      ),
      colour = "blue",
      linewidth = 0.4,
      lineend = "square",
      inherit.aes = FALSE
    ) +
    xlab(chr) +
    scale_x_continuous(
      breaks = seq(
        from = region_start - cand_padding,
        to = region_end + cand_padding,
        length.out = max(2, ceiling(width_mm / 40))
      ),
      expand = expansion(mult = c(0.01, 0.05)),
      guide = guide_axis(cap = "both")
    ) +
    scale_y_continuous(
      guide = guide_axis(cap = TRUE),
      breaks = seq(
        from = 0,
        to = ceiling(max(region_scans$LOGPVALUE)),
        length.out = ceiling(height_mm / 7.5)
      ),
      labels = seq(
        from = 0,
        to = ceiling(max(region_scans$LOGPVALUE)),
        length.out = ceiling(height_mm / 7.5)
      ) |> round(digits = 0),
      limits = c(
        constant_y_min * region_yratio,
        ceiling(max(region_scans$LOGPVALUE))
      )
    ) +
    axislabels_common +
    scales_common +
    theme_common +
    theme(
      axis.title.x = element_text(
        colour = "black",
        face = "bold",
        size = textsize,
        hjust = 0.48
      )
    )

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
      cand_index, "/", nrow(parsed_cands),
      " saved to ", candplot_outname,
      " and candidate genes extracted from your .gff",
      sep = ""
    )
  )

}

print("Done")
