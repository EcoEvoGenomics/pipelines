suppressMessages(library(tidyverse))
suppressMessages(library(gtable))
suppressMessages(library(grid))
suppressMessages(library(ggrepel))

# -------- Parse command line arguments ----------------------------------------
args <- commandArgs(trailing = TRUE)
if (length(args) != 6) {
  stop("USAGE: Rscript plot_xpehh.R
  <path to list of xpehh/ihs.csv file paths>
  <path to list of xpehh/ihs.cand.csv file paths>
  <path to annotations .gff>
  <base p-value used to call candidate regions (e.g. 0.01)>
  <plot width in mm>
  <plot height per scan in mm>"
  )
}
input_scans <- args[1]
input_cands <- args[2]
input_gff <- args[3]
base_pval <- as.double(args[4])
width_mm <- as.integer(args[5])
height_mm <- as.integer(args[6])

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

for (i in seq_len(n_cands)) {

  file_path <- readLines(input_cands)[i]

  scan_number <- rep(seq_len(n_scans), n_cands / n_scans)[i]
  cand_number <- ceiling(i / n_scans)

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

format_and_filter <- function(scans_or_cands) {
  scans_or_cands |>
    mutate(SCAN = factor(SCAN, levels = unique(SCAN))) |>
    filter(CHR != "chrLGE22") |>
    filter(CHR != "chrZ") |>
    filter(CHR != "mtDNA") |>
    mutate(
      CHR = factor(
        CHR,
        levels = paste("chr", c(1, "1A", seq(2, 28)[-15]), sep = "")
      )
    )
}

parsed_scans <- format_and_filter(parsed_scans)
parsed_cands <- format_and_filter(parsed_cands)

print("Input files parsed")

# -------- Parse annotations from GFF file (courtesy of Jack Harper) -----------
print("Parsing annotations from GFF ...")

annots <- read_tsv(
  input_gff,
  col_names = c(
    "CHR", "SOURCE", "FEATURE",
    "START", "END",
    "SCORE", "STRAND", "FRAME",
    "ATTRIBUTE"
  ),
  show_col_types = FALSE
) |>
  filter(FEATURE == "gene") |>
  mutate(HOMOLOG = str_extract(ATTRIBUTE, "(?<=Note=)[^;]+")) |>
  mutate(SYMBOL = str_extract(HOMOLOG, "(?<=Similar to )[^:]+")) |>
  mutate(SYMBOL = toupper(SYMBOL)) |>
  mutate(
    DIRECTED_SYMBOL = ifelse(
      STRAND == "+",
      paste(SYMBOL, ">"), # For an actual arrow, use "\u2192"
      paste("<", SYMBOL)  # For an actual arrow, use "\u2190"
    )
  )

print("GFF parsed")

# -------- Plotting constants --------------------------------------------------
constant_y_max <- max(parsed_scans$LOGPVALUE)
n_cand_searches <- max(parsed_cands$CAND)
constant_y_min <- -(constant_y_max / 10)
y_gap <- constant_y_min / 2
cand_padding <- 5e5

chr_names <- paste("chr", c(1, "1A", seq(2, 28)[-15]), sep = "")
chr_labels <- c(1, "1A", seq(2, 28)[-15])
chr_labels[21:28] <- ""
names(chr_labels) <- chr_names

scan_names <- unique(parsed_scans$SCAN)
scan_labels <- basename(readLines(input_scans))
names(scan_labels) <- scan_names

facet_labels <- c(chr_labels, scan_labels)

theme_common <- theme(
  axis.line.y = element_line(colour = "grey45", linewidth = 0.1),
  axis.line.x = element_blank(),
  axis.text.y = element_text(colour = "grey45", size = 4),
  axis.text.x = element_blank(),
  axis.ticks.y = element_line(colour = "grey45", linewidth = 0.1),
  axis.ticks.length.y = unit(0.5, "mm"),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(colour = "grey45", size = 4, face = "bold"),
  axis.title.x = element_text(
    colour = "grey45", size = 4, face = "bold",
    margin = margin(0, 0, 0, 0, unit = "mm")
  ),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.spacing = unit(0, "mm"),
  strip.background = element_blank(),
  strip.text.x = element_text(
    colour = "grey45",
    size = 4,
    margin = margin(0, 0, 0, 0, unit = "mm")
  ),
  strip.text.y = element_blank()
)

# -------- Plot main Manhattan plot --------------------------------------------

print("Plotting main plot ...")

# Draw Manhattan plots
manhattan <- ggplot(parsed_scans, aes(x = POSITION, y = LOGPVALUE)) +
  coord_cartesian(clip = "on") +
  facet_grid(
    cols = vars(CHR),
    rows = vars(SCAN),
    scales = "free_x",
    space = "free_x",
    switch = "x",
    labeller = as_labeller(facet_labels)
  ) +
  geom_point(alpha = 0.1, size = 0.25, stroke = NA, show.legend = FALSE) +
  geom_point(
    aes(alpha = LOGPVALUE > SCAN_BONFERRONI_THRESHOLD),
    colour = "red", size = 0.25, stroke = NA, show.legend = FALSE
  ) +
  labs(
    x = "Chromosome",
    y = expression(bold(bolditalic("-log")[10] ~ "P"))
  ) +
  scale_alpha_manual(values = c(0, 1)) +
  scale_y_continuous(
    guide = guide_axis(cap = TRUE),
    breaks = seq(
      from = 0,
      to = ceiling(constant_y_max * 0.8),
      length.out = ceiling(height_mm / 7.5)
    ),
    limits = c(constant_y_min, ceiling(constant_y_max) * 1.1)
  ) +
  scale_x_continuous(expand = expansion(add = 3e6)) +
  theme_common

# Draw candidate regions from candidate region scan
manhattan <- manhattan +
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
  )

# Draw boxes around candidate region scans
manhattan <- manhattan +
  geom_rect(
    data = parsed_scans |>
      group_by(CHR) |>
      summarise(START = min(POSITION), END = max(POSITION)),
    aes(xmin = START, xmax = END),
    ymax = y_gap,
    ymin = constant_y_min,
    inherit.aes = FALSE,
    colour = "grey45", fill = NA, linewidth = 0.1
  )

# Add horizontal labels for each scan above the y axis
add_scan_labels <- function(plot, labels) {

  g <- ggplotGrob(plot)

  panels <- g$layout[grepl("^panel", g$layout$name), , drop = FALSE]
  panels <- panels[order(panels$t, panels$l), ]
  panel_row_t <- unique(panels$t)
  row_ranges <- lapply(panel_row_t, function(tval) {
    rp <- panels[panels$t == tval, ]
    list(l = min(rp$l), r = max(rp$r))
  })

  for (i in rev(seq_along(panel_row_t))) {
    tval <- panel_row_t[i]
    lcol <- row_ranges[[i]]$l
    rcol <- row_ranges[[i]]$r

    g <- gtable_add_rows(g, heights = unit(height_mm / 5, "mm"), pos = tval - 1)

    lab_grob <- textGrob(
      label = labels[i],
      x = unit(0, "npc"),
      gp = gpar(fontsize = 4, fontface = "bold", col = "grey45"),
      just = "left", vjust = unit(height_mm / 5, "mm")
    )

    g <- gtable_add_grob(
      g,
      grobs = lab_grob,
      t = tval, l = lcol, b = tval, r = rcol,
      name = paste0("row-label-", i),
      clip = "off",
      z = Inf
    )
  }

  return(g)
}

manhattan <- add_scan_labels(manhattan, scan_labels)

# -------- Output main plot ----------------------------------------------------
outname <- paste(tools::file_path_sans_ext(input_scans), ".png", sep = "")
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

print("Plotting candidate region plots ...")

for (i in seq_len(nrow(parsed_cands))) {

  chr <- parsed_cands$CHR[i]
  focal_scan <- parsed_cands$SCAN[i]
  region_start <- parsed_cands$START[i]
  region_end <- parsed_cands$END[i]

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

  region_labelled_annots <- region_annots |>
    filter(START > region_start) |>
    filter(END < region_end) |>
    mutate(SCAN = focal_scan)

  candhattan <- ggplot(region_scans, aes(x = POSITION, y = LOGPVALUE)) +
    coord_cartesian(clip = "on") +
    facet_grid(rows = vars(SCAN), labeller = as_labeller(facet_labels)) +
    geom_rect(
      data = region_cands,
      aes(xmin = START, xmax = END, ymin = constant_y_min, ymax = y_gap),
      inherit.aes = FALSE,
      fill = "red"
    ) +
    geom_rect(
      data = region_annots,
      aes(xmin = START, xmax = END, ymin = constant_y_min, ymax = y_gap),
      inherit.aes = FALSE,
      fill = "grey15"
    ) +
    annotate(
      "rect",
      xmin = region_start - cand_padding,
      xmax = region_end + cand_padding,
      ymin = constant_y_min,
      ymax = y_gap,
      colour = "grey45",
      fill = NA,
      linewidth = 0.1
    ) +
    geom_text_repel(
      data = region_labelled_annots,
      aes(x = END - ((END - START) / 2), y = y_gap, label = DIRECTED_SYMBOL),
      ylim = c(constant_y_max * 0.4, constant_y_max * 0.8),
      box.padding = 0.1,
      direction = "y",
      min.segment.length = 0,
      segment.size = 0.05,
      segment.color = "grey85",
      colour = "grey70",
      size = 1,
      inherit.aes = FALSE
    ) +
    geom_point(
      alpha = 0.1, fill = "black", stroke = NA, pch = 21, size = 0.6,
      show.legend = FALSE
    ) +
    geom_point(
      aes(alpha = LOGPVALUE > SCAN_BONFERRONI_THRESHOLD),
      fill = "red", stroke = NA, pch = 21, size = 0.6,
      show.legend = FALSE
    ) +
    labs(
      y = expression(bold(bolditalic("-log")[10] ~ "P")),
      x = paste("Chromosome", chr_labels[which(chr_names == chr)])
    ) +
    scale_alpha_manual(values = c(0, 1)) +
    scale_y_continuous(
      breaks = seq(
        from = 0,
        to = ceiling(constant_y_max * 0.8),
        length.out = ceiling(height_mm / 7.5)
      ),
      expand = expansion(0),
      guide = guide_axis(cap = TRUE),
      limits = c(constant_y_min, ceiling(constant_y_max) * 1.1)
    ) +
    scale_x_continuous(
      breaks = seq(
        from = region_start - cand_padding,
        to = region_end + cand_padding,
        length.out = ceiling(width_mm / 80)
      ),
      expand = expansion(add = cand_padding / 10)
    ) +
    theme_common +
    theme(
      axis.text.x = element_text(colour = "grey45", size = 4),
      axis.ticks.x = element_line(colour = "grey45", linewidth = 0.1),
      axis.ticks.length.x = unit(0.5, "mm"),
      axis.title.x = element_text(
        colour = "grey45", size = 4, face = "bold",
        margin = margin(2, 0, 0, 0, unit = "mm")
      )
    )

  candhattan <- add_scan_labels(candhattan, scan_labels)

  outname <- paste(
    tools::file_path_sans_ext(input_scans), "/",
    chr, "_", format(region_start, scientific = FALSE),
    "_", format(region_end, scientific = FALSE), ".png",
    sep = ""
  )
  ggsave(
    candhattan,
    filename = outname,
    create.dir = TRUE,
    dpi = 1600,
    units = "mm",
    width = width_mm,
    height = height_mm * length(unique(parsed_scans$SCAN))
  )
  print(
    paste(
      "Candidate region plot ",
      i, "/", nrow(parsed_cands),
      " saved to ", outname,
      sep = ""
    )
  )

}

print("Done")
