# =====================================================================================================
# Rasterised skin-dose map - replaces the plotly mesh3d + headless-Chrome screenshot
# =====================================================================================================
# The report's Figure 2 used to be built by rendering a 4-trace plotly `mesh3d` to an HTML widget
# (`htmlwidgets::saveWidget`, whose own comment in report.Rmd read "this consumes lots of RAM"), booting
# Google Chrome, waiting for plotly.js to establish a WebGL context, screenshotting the page and cropping
# it. That is why the Dockerfile carried Chrome and fourteen X/GTK shared libraries: ~400 MB of image, and
# several seconds plus a browser process per report.
#
# None of it was needed. The camera sat at (0, 0, big) looking straight down z with every axis hidden, so
# the "3D" render was an orthographic projection of a triangle mesh - something base R draws directly.
# This file rasterises the same mesh with a PAINTER'S ALGORITHM: transform the vertices per view, order
# the triangles far-to-near, and hand all of them to one vectorised `polygon()` call. Measured on the
# ICRP145 adult male (41,002 vertices / 82,000 triangles): **0.14 s per view**, no browser, no WebGL, no
# GL context, no new system dependency of any kind.
#
# Three things about it are load-bearing:
#
#   * **No back-face culling, deliberately.** Culling would halve the triangles, but it depends on the mesh
#     having consistent outward winding, and two of the four views apply a REFLECTION (`-z`, `-y`), which
#     flips the winding sign. Getting that wrong punches holes in the phantom rather than failing loudly.
#     Painter's algorithm is correct whatever the winding - back faces are simply drawn first and painted
#     over - and at 0.14 s per view there is nothing to buy.
#
#   * **FLAT shading, and no lighting.** Each triangle takes one colour, the mean of its three vertices'
#     dose. plotly Gouraud-interpolates, but at 82,000 triangles across a human body a triangle is well
#     under a pixel at report resolution, so the difference is invisible. Lighting is deliberately absent:
#     this is a QUANTITATIVE dose map read against a colour bar, and shading a surface by its angle to a
#     light source makes the same dose render as two different colours. Form is conveyed by the silhouette
#     and the four views, not by shading.
#
#   * **Intensity is PER VERTEX.** The phantom data.frame packs vertices and faces into one frame: x/y/z are
#     populated for the first `nv` rows only, while i/j/k (0-indexed) run the full height. plotly's default
#     `intensitymode` is "vertex", so it reads the first `nv` values of `cumulative_skin_dose` and ignores
#     the rest - this file does the same, and `sdm_mesh()` is the single place that unpacking happens.
#
# Depends on nothing but base R + grDevices. `ragg` is used when present (faster, better text) and
# `grDevices::png(type = "cairo")` is the fallback, so this runs on any rocker image unchanged.

# --- The colour scale -------------------------------------------------------------------------------
# **VIRIDIS IS THE DEFAULT because it is what the existing reports are drawn in.** report.Rmd never set
# `colorscale`, and R's plotly binding defaults a `mesh3d` intensity to Viridis - confirmed by pulling the
# figure straight out of a report built by the old pipeline, not by reading the docs. (plotly.js's *own*
# default is the blue-grey-red ramp kept below as SDM_PLOTLY_DEFAULT; that is NOT what the R binding
# produced here, and defaulting to it would have silently re-coloured every report.)
#
# It is also the right scale on the merits - perceptually uniform and colour-blind-safe, which a
# blue-through-grey-to-red diverging ramp is not, and this is a quantitative map read against a bar.
SDM_VIRIDIS <- "viridis"

SDM_PLOTLY_DEFAULT <- list(
  stops  = c(0, 0.35, 0.5, 0.6, 0.7, 1),
  colours = c("#050AAC", "#6A89F7", "#BEBEBE", "#DCAA84", "#E6915A", "#B20A1C"))

# Returns a function mapping t in [0, 1] to a hex colour. `scale` is either the name of an hcl.colors
# palette ("viridis", "Inferno", ...) or a list of `stops` + `colours` like SDM_PLOTLY_DEFAULT.
sdm_colour_fn <- function(scale = SDM_PLOTLY_DEFAULT) {
  if (is.character(scale)) {
    pal <- grDevices::hcl.colors(256, scale)
    return(function(t) pal[pmax(1L, pmin(256L, as.integer(t * 255) + 1L))])
  }
  # colorRamp interpolates in RGB between arbitrary stops, which is what a plotly colorscale is.
  ramp <- grDevices::colorRamp(scale$colours, space = "Lab")
  function(t) {
    m <- ramp(pmax(0, pmin(1, stats::approx(scale$stops, seq(0, 1, length.out = length(scale$stops)),
                                            xout = t, rule = 2)$y)))
    grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
  }
}

# A dose for a colour-bar end label. `signif()` alone is not enough: the floor of a real dose field is a
# denormal like 7.74e-18, and formatC(format = "fg") renders that as "0.000000000000000007740" - twenty
# digits of noise, which is what the first cut of this bar printed. Anything below a ten-thousandth of the
# peak is zero as far as a skin-dose report is concerned, so say "0".
sdm_fmt <- function(v, ref = v) {
  if (!is.finite(v)) return("-")
  if (abs(v) < max(abs(ref), .Machine$double.eps) * 1e-4) return("0")
  formatC(signif(v, 3), format = "g", drop0trailing = TRUE)
}

# --- Unpack the phantom frame into a mesh -----------------------------------------------------------
# One place, because the packing is non-obvious (see the header) and three callers would each get it
# subtly wrong. Returns 1-INDEXED faces, ready for R's own subsetting.
sdm_mesh <- function(phantom) {
  nv <- sum(!is.na(phantom$x))
  if (nv < 3L || !nrow(phantom)) stop("phantom carries no vertices")
  V <- cbind(x = phantom$x[seq_len(nv)], y = phantom$y[seq_len(nv)], z = phantom$z[seq_len(nv)])
  Fc <- cbind(phantom$i, phantom$j, phantom$k) + 1L          # plotly's i/j/k are 0-indexed
  Fc <- Fc[stats::complete.cases(Fc) & apply(Fc, 1, max) <= nv, , drop = FALSE]
  d  <- phantom$cumulative_skin_dose[seq_len(nv)]
  d[!is.finite(d)] <- 0
  list(V = V, F = Fc, intensity = d)
}

# --- The four views ---------------------------------------------------------------------------------
# Exactly the four report.Rmd built, and in the same anatomical sense: each entry names the surface facing
# the reader, the in-plane horizontal axis, and the DEPTH axis (larger = nearer the camera, so drawing in
# increasing depth order paints near over far).
#
# The ARRANGEMENT is the caller's choice (`ncol`), and only the arrangement: the projection, the depth axis
# and the visible surface of each panel are the Rmd's. `ncol = 4` (the default) puts the four side by side,
# which on a portrait page draws each body TALLER than the old 2x2 did while using little more than half the
# page height - so Figure 2 and the Discussion now share a page. `ncol = 2` reproduces the old 2x2 exactly
# for anyone who wants the familiar layout.
SDM_VIEWS <- list(
  list(name = "Posterior",     h = function(V)  V[, "y"], depth = function(V)  V[, "z"]),
  list(name = "Right lateral", h = function(V)  V[, "z"], depth = function(V)  V[, "y"]),
  list(name = "Left lateral",  h = function(V) -V[, "z"], depth = function(V) -V[, "y"]),
  list(name = "Anterior",      h = function(V) -V[, "y"], depth = function(V) -V[, "z"]))

# --- The renderer -----------------------------------------------------------------------------------
# Writes a PNG and returns its path. The image's ASPECT is derived from the content, and the pixel height
# from `width_in` x `res` - so the caller gets a picture with no letterboxing to place at its natural size,
# which is the other half of why Figure 2 used to be mostly blank.
#
# `peak_label = TRUE` prints the peak dose in the corner, as the app's on-screen figure does via its
# plotly annotation. Panel captions are always drawn: the old figure had none, so a reader had to know the
# convention to tell an anterior view from a posterior one.
skin_dose_map_png <- function(phantom, file,
                              width_in = 6.2, res = 200,
                              max_height_in = 8.4,
                              scale = SDM_VIRIDIS,
                              ncol = 4L,
                              units_label = "Gy",
                              peak_label = TRUE,
                              gap_frac = 0.14) {
  m   <- sdm_mesh(phantom)
  V   <- m$V; Fc <- m$F; inten <- m$intensity
  rng <- range(inten, na.rm = TRUE)
  # A flat dose field (an unrun calculation, or a genuinely uniform one) has no range to normalise by;
  # map it all to the bottom of the ramp rather than dividing by zero.
  norm <- if (diff(rng) > 0) function(v) (v - rng[1]) / diff(rng) else function(v) rep(0, length(v))
  colf <- sdm_colour_fn(scale)
  # Per-face colour: the mean of the face's three vertex doses (flat shading - see the header).
  face_t   <- (inten[Fc[, 1]] + inten[Fc[, 2]] + inten[Fc[, 3]]) / 3
  face_col <- colf(norm(face_t))

  vert  <- V[, "x"]                                  # the phantom's long axis - drawn vertically, head up
  vr    <- range(vert)
  vspan <- diff(vr)
  cap_h <- vspan * 0.055                             # room under each panel for its view caption
  bar_h <- vspan * 0.075                             # room under the last row for the colour bar
  np    <- length(SDM_VIEWS)
  nc    <- max(1L, min(as.integer(ncol), np))
  nr    <- ceiling(np / nc)
  # Every panel is given the WIDEST panel's cell width so the bodies sit on a common scale: a lateral view
  # is genuinely narrower than a posterior one and has to LOOK narrower.
  cell_w <- max(vapply(SDM_VIEWS, function(vw) diff(range(vw$h(V))), numeric(1))) * (1 + gap_frac)
  cell_h <- vspan + cap_h * 1.6
  panels <- lapply(seq_len(np), function(p) {
    vw  <- SDM_VIEWS[[p]]
    hv  <- vw$h(V)
    col <- (p - 1L) %%  nc
    row <- (p - 1L) %/% nc
    list(name = vw$name,
         h    = hv - mean(range(hv)) + col * cell_w,     # centre each panel in its own cell
         v    = vert - mean(vr) - row * cell_h,
         d    = vw$depth(V),
         cx   = col * cell_w,
         vbot = -vspan / 2 - row * cell_h)
  })

  world_x <- c(-cell_w / 2, (nc - 1) * cell_w + cell_w / 2)
  # Headroom at the top for the peak-dose label, which hangs from world_y[2] - without it the text's
  # ascenders are clipped by the device edge.
  world_y <- c(-vspan / 2 - (nr - 1) * cell_h - cap_h - bar_h, vspan / 2 + vspan * 0.05)

  # The image's aspect follows the CONTENT, so it is placed at its natural size with no letterboxing - the
  # other half of why the old Figure 2 sat in a mostly-empty box. `max_height_in` is the page's usable
  # height: a tall arrangement (ncol = 2 stacks two full-length bodies, ~3:1) would otherwise come out
  # 19 in tall at 6.2 in wide and be scaled down by Word with no say in it.
  height_in <- width_in * diff(world_y) / diff(world_x)
  if (is.finite(max_height_in) && height_in > max_height_in) {
    width_in  <- width_in * max_height_in / height_in
    height_in <- max_height_in
  }
  dev_open <- function(f, w, h, r) {
    if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png(f, width = w, height = h, res = r, background = "white")
    else grDevices::png(f, width = w, height = h, res = r, type = "cairo", bg = "white")
  }
  dev_open(file, round(width_in * res), round(height_in * res), res)
  on.exit(grDevices::dev.off(), add = TRUE)

  op <- graphics::par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot.new()
  graphics::plot.window(world_x, world_y, asp = 1)

  for (pn in panels) {
    # Painter's algorithm: order the triangles by increasing depth so nearer ones paint over further ones.
    dep <- (pn$d[Fc[, 1]] + pn$d[Fc[, 2]] + pn$d[Fc[, 3]]) / 3
    ord <- order(dep)
    xs  <- rbind(pn$h[Fc[ord, 1]], pn$h[Fc[ord, 2]], pn$h[Fc[ord, 3]], NA_real_)
    ys  <- rbind(pn$v[Fc[ord, 1]], pn$v[Fc[ord, 2]], pn$v[Fc[ord, 3]], NA_real_)
    # ONE polygon() call for the whole panel. Called per triangle this is minutes, not milliseconds.
    #
    # `border` is the FILL COLOUR, not NA. With border = NA the device antialiases each triangle against
    # the white background independently, so every shared edge keeps a hairline of background - across
    # 82,000 triangles that reads as a honeycomb texture over the whole phantom, heaviest exactly where the
    # dose gradient is flat. Stroking each triangle in its own fill colour closes the seams and changes
    # nothing else.
    graphics::polygon(as.vector(xs), as.vector(ys), col = face_col[ord],
                      border = face_col[ord], lwd = 0.4)
    graphics::text(pn$cx, pn$vbot - cap_h * 0.55, pn$name, cex = 0.62, col = "#333333")
  }

  # --- colour bar -------------------------------------------------------------------------------------
  # The end labels sit BESIDE the bar, not under it. Underneath they were clipped by the bottom of the
  # device (the reserved band has to hold the bar AND two lines of text), and moving them alongside costs
  # nothing but a narrower bar.
  bx <- c(world_x[1] + diff(world_x) * 0.30, world_x[2] - diff(world_x) * 0.30)
  by <- c(world_y[1] + bar_h * 0.28, world_y[1] + bar_h * 0.60)
  n  <- 256L
  edges <- seq(bx[1], bx[2], length.out = n + 1L)
  graphics::rect(edges[-(n + 1L)], by[1], edges[-1L], by[2],
                 col = colf(seq(0, 1, length.out = n)), border = NA)
  graphics::rect(bx[1], by[1], bx[2], by[2], border = "#888888", lwd = 0.5)
  pad <- diff(world_x) * 0.012
  graphics::text(bx[1] - pad, mean(by), sdm_fmt(rng[1], rng[2]), cex = 0.58, col = "#333333", adj = c(1, 0.5))
  graphics::text(bx[2] + pad, mean(by), sdm_fmt(rng[2], rng[2]), cex = 0.58, col = "#333333", adj = c(0, 0.5))
  graphics::text(mean(bx), by[2] + bar_h * 0.14, units_label, cex = 0.6, col = "#333333", adj = c(0.5, 0))
  if (isTRUE(peak_label))
    graphics::text(world_x[2], world_y[2],
                   sprintf("Peak skin dose: %s %s", sdm_fmt(rng[2], rng[2]), units_label),
                   cex = 0.62, col = "#333333", adj = c(1, 1))
  # The size comes back so the caller can place the image at its true aspect without re-reading the PNG
  # (which would mean depending on the `png` package just to recover a number computed right here).
  invisible(list(file = file, width_in = width_in, height_in = height_in))
}
