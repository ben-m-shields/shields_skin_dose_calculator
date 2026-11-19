# =====================================================================================================
# Peak Skin Dose Estimate - Word (.docx) builder
# =====================================================================================================
# Replaces report.Rmd (rmarkdown -> pandoc -> word_document). The document's structure, headings, wording
# and references are unchanged - it is a clinical report a physicist signs, so it is reproduced faithfully
# rather than redesigned - but it is now built directly with `officer` + `flextable`, which means:
#
#   * no pandoc, and no rmarkdown render of a whole document per download;
#   * no `htmlwidgets::saveWidget` (report.Rmd's own comment: "this consumes lots of RAM"), no Google
#     Chrome and no `webshot`. Figure 2 comes from skin_dose_map.R, which rasterises the same mesh in R;
#   * the two tables are REAL Word tables rather than a `kable` and a screenshot, so they are selectable,
#     searchable and restyleable by the recipient.
#
# Four things are deliberately better than the Rmd, and each is a defect fixed rather than a preference:
#
#   * **Table 1 is keyed by COLUMN NAME.** The Rmd did `t(select(summary, -Patient_ID, ...))` and then
#     assigned 13 `rownames` positionally, so the day a column was added, reordered or renamed upstream
#     every row of the patient's dose table would be silently mislabelled - a wrong number against a wrong
#     name in a signed clinical document. SDR_SUMMARY_ROWS maps name -> label, so a missing column drops
#     out and nothing is ever misattributed.
#   * **Table 2 is the skin_reactions DATA, not skin_reactions.png.** The tibble is already in
#     environment.RData; the report was embedding a picture of it. As a real table it also lets the row
#     matching THIS patient's peak dose be highlighted, so the reader is not left to find the band by eye.
#   * **A page-numbered footer and a running header**, with live Word PAGE/NUMPAGES fields, because this
#     is a document that gets printed and filed.
#   * **A ruled signature block** in place of the literal text "Insert signature here".
#
# Depends on officer + flextable + ggplot2 + skin_dose_map.R. It reads app.R nothing: the caller passes the
# phantom, the summary frame and the patient identifiers in.

if (!exists("skin_dose_map_png", mode = "function")) source("skin_dose_map.R")

# Picked up from MIRI's branding.R when this app is folded into it; standalone defaults otherwise, so the
# file works in both trees without an edit.
if (!exists("BRAND_FONT")) BRAND_FONT <- "Aptos"     # Word 2024+ default body font
if (!exists("APP_ACCENT")) APP_ACCENT <- "#a6093d"

SDR_TITLE <- "Peak Skin Dose Estimate"

# Table 1's rows: stored column -> printed label, IN PRINT ORDER. Keyed by name on purpose - see the header.
SDR_SUMMARY_ROWS <- c(
  Study_Date                          = "Study Date",
  Accession_Number                    = "Accession Number",
  Study_Description                   = "Study Description",
  Performing_Physician                = "Performing Physician",
  Referring_Physician                 = "Referring Physician",
  Manufacturer                        = "Manufacturer",
  Model                               = "Model",
  Serial_Number                       = "Serial Number",
  Fluoro_Dose_RP_Total                = "Total Fluoro Dose (Gy)",
  Fluoro_Dose_Area_Product_Total      = "Total Fluoro DAP (Gy.cm2)",
  Total_Fluoro_Time                   = "Total Fluoro Time (s)",
  Acquisition_Dose_RP_Total           = "Total Acquisition Dose (Gy)",
  Acquisition_Dose_Area_Product_Total = "Total Acquisition DAP (Gy.cm2)",
  Total_Acquisition_Time              = "Total Acquisition Time (s)")

# --- small officer helpers (the same shapes MIRI's report_docx.R uses, so the two merge cleanly) --------
sdr_h1   <- function(doc, txt) officer::body_add_par(doc, txt, style = "heading 1")
sdr_body <- function(doc, txt) officer::body_add_fpar(doc, officer::fpar(
  officer::ftext(txt, officer::fp_text(font.family = BRAND_FONT, font.size = 11, color = "#333333"))))
sdr_caption <- function(doc, txt) officer::body_add_fpar(doc, officer::fpar(
  officer::ftext(txt, officer::fp_text(font.family = BRAND_FONT, font.size = 10, italic = TRUE, color = "#555555")),
  fp_p = officer::fp_par(text.align = "center", padding.top = 4)))

# A page break that takes up NO ROOM ITSELF. `officer::body_add_break()` adds a full-height paragraph
# carrying the break, and when the preceding content happens to fill the page exactly that paragraph no
# longer fits either - so it lands on the NEXT page and its break then pushes the following content one
# page further, leaving a completely blank page behind. That is not a hypothetical tuning problem: Table 1
# grows a row per procedure, so any given report can land on that knife edge. A near-zero line height means
# the break paragraph always fits wherever it is, and the break does only what it was asked to do.
sdr_page_break <- function(doc) officer::body_add_fpar(doc, officer::fpar(
  officer::run_pagebreak(), fp_p = officer::fp_par(line_spacing = 0.05, padding = 0)))

sdr_header_fpar <- function() officer::fpar(
  officer::ftext(SDR_TITLE, officer::fp_text(font.family = BRAND_FONT, font.size = 9, color = "#666666")),
  fp_p = officer::fp_par(text.align = "left", padding.bottom = 3,
                         border.bottom = officer::fp_border(color = APP_ACCENT, width = 1)))

sdr_footer_fpar <- function() {
  g <- officer::fp_text(font.family = BRAND_FONT, font.size = 9, color = "#666666")
  officer::fpar(officer::ftext("Page ", g), officer::run_word_field("PAGE", prop = g),
                officer::ftext(" of ", g), officer::run_word_field("NUMPAGES", prop = g),
                fp_p = officer::fp_par(text.align = "center", padding.top = 3,
                                       border.top = officer::fp_border(color = APP_ACCENT, width = 0.5)))
}

# A branded flextable squeezed to the A4 content width. autofit() (never set_table_properties(layout =
# "autofit", width = 1), which emits an EMPTY <w:tblGrid> that Word rejects as malformed) then proportional
# width() (never fit_to_width(), which fits by SHRINKING THE FONT below the 11pt floor).
sdr_flextable <- function(df, width = 6.2, header = TRUE, fontsize = 11, widths = NULL) {
  ft <- flextable::flextable(df)
  if (!header) ft <- flextable::delete_part(ft, part = "header")
  ft <- flextable::font(ft, fontname = BRAND_FONT, part = "all")
  ft <- flextable::fontsize(ft, size = fontsize, part = "all")
  if (header) {
    ft <- flextable::bg(ft, bg = APP_ACCENT, part = "header")
    ft <- flextable::color(ft, color = "white", part = "header")
    ft <- flextable::bold(ft, part = "header")
  }
  ft <- flextable::border_inner(ft, border = officer::fp_border(color = "#d9d9d9", width = 0.5), part = "body")
  ft <- flextable::valign(ft, valign = "top", part = "body")
  ft <- flextable::padding(ft, padding = 3, part = "all")
  ft <- flextable::autofit(ft)
  # An explicit width vector wins outright. Content-proportional widths are wrong for a table of PARALLEL
  # categories (the four reaction-window columns), where autofit sizes each by how wordy that window's
  # prose happens to be and the narrow ones then break words mid-word ("desquamatio / n").
  if (!is.null(widths)) return(flextable::width(ft, width = widths * (width / sum(widths))))
  w <- as.numeric(dim(ft)$widths)
  if (length(w) && sum(w) > width) {
    w <- w * (width / sum(w))
    # A PROPORTIONAL squeeze can crush a narrow column below the width of one character, and Word then
    # wraps its text ONE LETTER PER LINE. That is not hypothetical: the skin-reaction table's first column
    # ("Skin dose (Gy)", "0 to 2") sits beside four columns of 200-character prose, so autofit gave it a
    # tiny share and the squeeze took that to nothing - the header rendered as a vertical letter ladder.
    # Floor every column, then take the shortfall back off whichever columns still have slack.
    fl    <- min(0.75, width / length(w))
    short <- sum(pmax(fl - w, 0))
    slack <- pmax(w - fl, 0)
    w     <- pmax(w, fl) - if (short > 0 && sum(slack) > 0) short * slack / sum(slack) else 0
    ft <- flextable::width(ft, width = w)
  }
  ft
}

# Plain text prompts, one per line - NOT a table. A ruled two-column table was tried and dropped: it is
# five short labels the author types beside, so the table added borders, cell padding and a fixed row height
# without adding anything the labels did not already say. Each prompt is followed by a blank line so there
# is room to type or sign under it.
SDR_SIGNATURE_PROMPTS <- c("Signature", "Name", "Position", "Location", "Date")

sdr_signature_block <- function(doc) {
  blank <- officer::fpar(officer::ftext("", officer::fp_text(font.family = BRAND_FONT, font.size = 11)))
  for (p in SDR_SIGNATURE_PROMPTS) {
    doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(
      paste0(p, ":"), officer::fp_text(font.family = BRAND_FONT, font.size = 11, color = "#333333"))))
    doc <- officer::body_add_fpar(doc, blank)
  }
  doc
}

# officer 0.7.6 writes <w:spacing w:line="240"/> with NO w:lineRule - at EVERY `fp_par(line_spacing =)`
# value, and there is no argument to change it. Both Word and LibreOffice then treat the line box as an
# EXACT 12pt, so any run larger than that overflows its own line: the 18pt title bled ABOVE the top margin
# (measured at 0.92in against a 0.98in margin) and Word, which clips at the text-area boundary instead of
# letting it bleed, sliced the tops off the letters.
#
# Adding w:lineRule="auto" makes 240 mean "single spacing", so the line box grows to fit its run. It is
# applied to EVERY paragraph, not just the title, because the same box is too tight for the 11pt body text
# as well (11pt Aptos wants ~13.2pt of leading, and was getting 12) - descenders were being shaved
# document-wide, just too subtly to notice next to a decapitated title.
#
# It patches the IN-MEMORY xml tree, not the file on disk: officer re-serialises document.xml from that
# tree inside print(), so an edit to the unzipped file would simply be overwritten (unlike settings.xml
# below, which officer does not regenerate).
sdr_fix_line_rule <- function(doc) {
  b <- tryCatch(officer::docx_body_xml(doc), error = function(e) NULL)
  if (is.null(b)) return(doc)
  n <- xml2::xml_find_all(b, "//w:spacing[@w:line and not(@w:lineRule)]")
  if (length(n)) xml2::xml_set_attr(n, "w:lineRule", "auto")

  # ...EXCEPT the page-break paragraphs, which get an EXACT ~1pt line box instead.
  #
  # `run_pagebreak()` takes no arguments, so its run carries the template's 11pt font, and under
  # lineRule="auto" a renderer gives that paragraph a FULL line (~0.19in) however small a line_spacing it
  # was asked for. When the preceding content had less than that left on the page, the break paragraph was
  # pushed to the next page and its break then pushed the real content one page further - a completely
  # BLANK page in the middle of the report, which is exactly what happened between Table 1 and Figure 1.
  # An exact 1pt box makes the paragraph negligible whatever font the run inherits; clipping it is harmless
  # because the run draws nothing but the break itself.
  brk <- xml2::xml_find_all(b, "//w:p[.//w:br[@w:type='page']]/w:pPr/w:spacing")
  if (length(brk)) {
    xml2::xml_set_attr(brk, "w:line", "20")          # 20 twentieths of a point = 1pt
    xml2::xml_set_attr(brk, "w:lineRule", "exact")
    xml2::xml_set_attr(brk, "w:before", "0")
    xml2::xml_set_attr(brk, "w:after", "0")
  }
  doc
}

# Tell Word to recompute PAGE / NUMPAGES when the reader opens the file. Without it officer writes the
# fields with an empty cached result and every page reads "Page  of".
sdr_update_fields_on_open <- function(doc) {
  sp <- tryCatch(file.path(doc$package_dir, "word", "settings.xml"), error = function(e) NULL)
  if (is.null(sp) || !file.exists(sp)) return(doc)
  x <- paste(readLines(sp, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (grepl("w:updateFields", x, fixed = TRUE)) return(doc)
  x2 <- sub("(<w:settings\\b[^>]*>)", "\\1<w:updateFields w:val=\"true\"/>", x, perl = TRUE)
  if (!identical(x2, x)) writeLines(x2, sp, useBytes = TRUE)
  doc
}

# --- content builders ---------------------------------------------------------------------------------
# Numbers straight out of the dose engine carry full double precision ("0.123456789012"), which is noise in
# a printed table. Round for display only; strings and dates pass through untouched.
sdr_fmt_cell <- function(v) {
  if (is.numeric(v)) return(ifelse(is.na(v), "-", formatC(signif(v, 4), format = "g", drop0trailing = TRUE)))
  v <- as.character(v); v[is.na(v)] <- "-"; v
}

# Table 1: one row per metric, one column per procedure included in the estimate.
sdr_summary_table <- function(summary_df) {
  have <- names(SDR_SUMMARY_ROWS)[names(SDR_SUMMARY_ROWS) %in% names(summary_df)]
  cells <- lapply(have, function(cn) sdr_fmt_cell(summary_df[[cn]]))
  out <- data.frame(Metric = unname(SDR_SUMMARY_ROWS[have]), check.names = FALSE, stringsAsFactors = FALSE)
  n <- nrow(summary_df)
  labs <- if (n == 1L) "Value" else paste("Procedure", seq_len(n))
  for (k in seq_len(n)) out[[labs[k]]] <- vapply(cells, function(x) x[[k]], character(1))
  out
}

# Figure 1: the dose-metric bar chart, unchanged from the Rmd (same facets, same fill, same theme).
sdr_summary_figure_data <- function(summary_df) {
  s <- summary_df
  data.frame(
    Value = c(sum(s$Fluoro_Dose_RP_Total), sum(s$Acquisition_Dose_RP_Total),
              sum(s$Fluoro_Dose_Area_Product_Total), sum(s$Acquisition_Dose_Area_Product_Total),
              sum(s$Total_Fluoro_Time), sum(s$Total_Acquisition_Time)),
    Metric = c("Dose at Reference Point (Gy)", "Dose at Reference Point (Gy)",
               "Dose Area Product (Gy.cm2)", "Dose Area Product (Gy.cm2)",
               "Screening Time (s)", "Screening Time (s)"),
    Irradiation_Type = c("Fluoroscopy", "Acquisition", "Fluoroscopy", "Acquisition",
                         "Fluoroscopy", "Acquisition"),
    check.names = FALSE, stringsAsFactors = FALSE)
}

sdr_summary_plot <- function(d) {
  ggplot2::ggplot(d, ggplot2::aes(x = .data$Irradiation_Type, y = .data$Value,
                                  fill = .data$Irradiation_Type, group = .data$Metric)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ Metric, ncol = 1, scales = "free") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none",
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", colour = APP_ACCENT))
}

# The Balter et al. band a peak dose falls in. Bands are "0 to 2" / "2 to 5" / ... / "15+", so the lower
# bound is the leading number; the match is the LAST band whose lower bound the dose reaches.
sdr_reaction_band <- function(skin_reactions, peak) {
  lo <- suppressWarnings(as.numeric(sub("^\\s*([0-9.]+).*$", "\\1", skin_reactions[[1]])))
  if (!is.finite(peak) || all(is.na(lo))) return(NA_integer_)
  hit <- which(!is.na(lo) & peak >= lo)
  if (!length(hit)) NA_integer_ else max(hit)
}

# --- the document -------------------------------------------------------------------------------------
# `patient` is a list(Patient_ID, DOB, Patient_Name, Gender). Writes `path`; returns it invisibly.
skin_dose_report_docx <- function(path, patient, summary_df, phantom, skin_reactions,
                                  map_scale = SDM_VIRIDIS, map_ncol = 4L) {
  peak <- max(phantom$cumulative_skin_dose, na.rm = TRUE)
  dob  <- tryCatch(format(as.Date(patient$DOB), "%d/%m/%Y"), error = function(e) as.character(patient$DOB))
  if (is.na(dob) || !nzchar(dob)) dob <- as.character(patient$DOB)

  doc <- officer::read_docx()
  # The title is a styled PARAGRAPH, not a "heading 1". officer's default template numbers heading 1, so a
  # title set as one comes out "1. Peak Skin Dose Estimate" and pushes Background to 2 - the document's own
  # name reading as its first numbered section.
  doc <- officer::body_add_fpar(doc, officer::fpar(
    officer::ftext(SDR_TITLE, officer::fp_text(font.family = BRAND_FONT, font.size = 18,
                                               bold = TRUE, color = APP_ACCENT)),
    fp_p = officer::fp_par(padding.bottom = 6,
                           border.bottom = officer::fp_border(color = APP_ACCENT, width = 1.5))))
  doc <- officer::body_add_fpar(doc, officer::fpar(
    officer::ftext(sprintf("Exposure of %s (%s), DOB: %s, Sex: %s",
                           patient$Patient_Name, patient$Patient_ID, dob, patient$Gender),
                   officer::fp_text(font.family = BRAND_FONT, font.size = 11, bold = TRUE))))

  # No blank spacer paragraphs between these sections: the "heading 1" style carries its own space-before,
  # and once sdr_fix_line_rule() gave every paragraph its proper leading the extra blanks pushed Table 1's
  # last row onto a second page.
  # ---- Background ----
  doc <- sdr_h1(doc, "Background")
  doc <- sdr_body(doc, paste(
    "The above patient has recently been exposed to a level of ionising radiation in excess of the",
    "Substantial Radiation Dose Level (SRDL) stipulated within the facility's Radiation Safety and",
    "Protection Plan (RSPP) during one or more fluoroscopically guided procedures. This report will aim to",
    "estimate the radiation dose to the patient's skin and identify the expected skin reaction."))

  # ---- Results ----
  doc <- sdr_h1(doc, "Results")
  doc <- officer::body_add_fpar(doc, officer::fpar(
    officer::ftext("The peak absorbed skin dose resulting from the above-mentioned procedure(s) is estimated to be approximately ",
                   officer::fp_text(font.family = BRAND_FONT, font.size = 11)),
    officer::ftext(sprintf("%s Gy", formatC(signif(peak, 2), format = "g", drop0trailing = TRUE)),
                   officer::fp_text(font.family = BRAND_FONT, font.size = 11, bold = TRUE, color = APP_ACCENT)),
    officer::ftext(".", officer::fp_text(font.family = BRAND_FONT, font.size = 11))))
  doc <- officer::body_add_par(doc, "")
  doc <- flextable::body_add_flextable(doc, sdr_flextable(sdr_summary_table(summary_df)))
  doc <- sdr_caption(doc, "Table 1. Summary of the procedure(s) included in this peak skin dose estimate.")
  doc <- sdr_page_break(doc)

  # ---- Figure 1: dose metrics ----
  doc <- officer::body_add_gg(doc, sdr_summary_plot(sdr_summary_figure_data(summary_df)),
                              width = 6.0, height = 7.2, res = 200)
  doc <- sdr_caption(doc, "Figure 1. Dose metrics included in this peak skin dose estimate.")
  doc <- sdr_page_break(doc)

  # ---- Figure 2: the skin dose map (rasterised in R - no browser; see skin_dose_map.R) ----
  png <- tempfile(fileext = ".png")
  on.exit(unlink(png), add = TRUE)
  fig <- skin_dose_map_png(phantom, png, width_in = 6.2, res = 200,
                           scale = map_scale, ncol = map_ncol)
  doc <- officer::body_add_img(doc, png, width = fig$width_in, height = fig$height_in)
  doc <- sdr_caption(doc, "Figure 2. Visualisation of the skin dose to the patient.")
  doc <- sdr_page_break(doc)

  # ---- Discussion (starts on its own page) ----
  doc <- sdr_h1(doc, "Discussion")
  doc <- sdr_body(doc, "See below table for an indication of the expected skin reaction.")
  doc <- officer::body_add_par(doc, "")
  sr  <- as.data.frame(skin_reactions, check.names = FALSE)
  # Dose band narrow, the four time windows equal - they are parallel categories, not competing for space.
  ft  <- sdr_flextable(sr, fontsize = 9, widths = c(0.85, rep(1.34, max(1, ncol(sr) - 1))))
  band <- sdr_reaction_band(sr, peak)
  if (!is.na(band)) {                                     # highlight the band THIS patient falls in
    ft <- flextable::bg(ft, i = band, bg = "#fdeaee", part = "body")
    ft <- flextable::bold(ft, i = band, j = 1, part = "body")
    ft <- flextable::color(ft, i = band, j = 1, color = APP_ACCENT, part = "body")
  }
  doc <- flextable::body_add_flextable(doc, ft)
  doc <- sdr_caption(doc, paste(
    "Table 2. Tissue reactions from single-delivery radiation dose to skin of the neck, torso, pelvis,",
    "buttocks, or arms (Balter et al, 2010).",
    if (!is.na(band)) "The shaded row is the band this patient's peak skin dose falls in." else ""))
  doc <- officer::body_add_par(doc, "")

  # ---- signature ----
  # The Rmd printed today's date under the block automatically; "Date" is one of the prompts now, so the
  # auto-filled line went with the table rather than sitting redundantly underneath its own prompt.
  doc <- sdr_signature_block(doc)
  doc <- sdr_page_break(doc)

  # ---- References (verbatim from report.Rmd) ----
  doc <- sdr_h1(doc, "References")
  for (r in SDR_REFERENCES) { doc <- sdr_body(doc, r); doc <- officer::body_add_par(doc, "") }

  # Running header + page-numbered footer on every page, and the A4 page the body was laid out for -
  # body_set_default_section REPLACES the section, so omitting page_size reflows the document to Letter.
  #
  # `header_first` is a BLANK block_list, which is what makes officer emit <w:titlePg/> and drop the header
  # from page 1 - the title there already reads "Peak Skin Dose Estimate", so the running header only
  # repeated it. The FOOTER is deliberately kept on page 1: it carries "Page 1 of N", which is the one thing
  # a printed clinical report needs on its first page.
  doc <- officer::body_set_default_section(doc, officer::prop_section(
    header_default = officer::block_list(sdr_header_fpar()),
    footer_default = officer::block_list(sdr_footer_fpar()),
    header_first   = officer::block_list(officer::fpar(officer::ftext(""))),
    footer_first   = officer::block_list(sdr_footer_fpar()),
    page_size    = officer::page_size(width = 8.27, height = 11.69, orient = "portrait"),
    page_margins = officer::page_mar(bottom = 0.98, top = 0.98, right = 0.98, left = 0.98,
                                     header = 0.49, footer = 0.49, gutter = 0)))
  # Must run after ALL content is added (it patches every paragraph in the body) and before print().
  doc <- sdr_fix_line_rule(doc)
  doc <- sdr_update_fields_on_open(doc)
  print(doc, target = path)
  invisible(path)
}

SDR_REFERENCES <- c(
  "Andersson, J. et al. (2021) 'Estimation of patient skin dose in fluoroscopy: Summary of a joint report by AAPM TG357 and EFOMP', Medical Physics, 48(7), pp. 671-696. doi:10.1002/mp.14910.",
  "Balter, S. et al. (2010) 'Fluoroscopically guided interventional procedures: A review of radiation effects on patients' skin and hair', Radiology, 254(2), pp. 326-341. doi:10.1148/radiol.2542082312.",
  "Benmakhlouf, H. et al. (2011) 'Backscatter factors and mass energy-absorption coefficient ratios for diagnostic radiology dosimetry', Physics in Medicine and Biology, 56(22), pp. 7179-7204. doi:10.1088/0031-9155/56/22/012.",
  "ICRP, 2020. Adult mesh-type reference computational phantoms. ICRP Publication 145. Ann. ICRP 49(3).",
  "ICRP, 2024. Paediatric mesh-type reference computational phantoms. ICRP Publication 156. Ann. ICRP 53(1-2).",
  "Khodadadegan, Y. et al. (2010) 'Automatic monitoring of localized skin dose with fluoroscopic and interventional procedures', Journal of Digital Imaging, 24(4), pp. 626-639. doi:10.1007/s10278-010-9320-7.",
  "Krajinovic, M., Krzanovic, N. and Ciraj-Bjelac, O. (2021) 'Vendor-independent skin dose mapping application for interventional radiology and Cardiology', Journal of Applied Clinical Medical Physics, 22(2), pp. 145-157. doi:10.1002/acm2.13167.",
  "Rana, V.K., Rudin, S. and Bednarek, D.R. (2016) 'A tracking system to calculate patient skin dose in real-time during neurointerventional procedures using a biplane x-ray imaging system', Medical Physics, 43(9), pp. 5131-5144. doi:10.1118/1.4960368.")
