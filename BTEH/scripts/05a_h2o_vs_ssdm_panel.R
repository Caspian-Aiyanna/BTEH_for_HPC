# ======================= BEFORE vs AFTER (H2O & SSDM) PANELS ==================
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(terra)
  library(sf)
  library(dplyr)
  library(grid)       # for unit()
})

# --------------------------- CONFIG (edit as needed) ---------------------------
# Root folders where your rasters live
h2o_root     <- Sys.getenv("H2O_ROOT",     unset = "results/H2O")
ssdm_root    <- Sys.getenv("SSDM_ROOT",    unset = "results/SSDM")
compare_root <- Sys.getenv("COMPARE_ROOT", unset = "results/compare/h2o_vs_ssdm")

# Optional AOI: a polygon/lines shapefile/GeoPackage; set to "" to disable
aoi_path     <- Sys.getenv("AOI_PATH", unset = "data/shp/HV20233.shp")

# Elephants to process (E? with A/B variants)
ele_ids <- c("E3","E4","E5")

# Quantiles for hotspot/GSL and metrics
qs <- c(0.25, 0.50, 0.75)

# ------------------------------ UTILITIES -------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

dir_ensure <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

# Convert SpatRaster (single-layer) to xy/value data.frame for ggplot
as_xy <- function(r) {
  stopifnot(inherits(r, "SpatRaster"))
  if (nlyr(r) > 1L) r <- r[[1]]
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x", "y", "val")
  df
}

# Jaccard/IoU on **binary SpatRasters** (TRUE/FALSE)
jaccard_binary <- function(b1, b2) {
  stopifnot(inherits(b1, "SpatRaster"), inherits(b2, "SpatRaster"))
  b1 <- b1 & !is.na(b1)
  b2 <- b2 & !is.na(b2)
  inter <- global(b1 & b2, "sum", na.rm = TRUE)[[1]]
  union <- global(b1 | b2, "sum", na.rm = TRUE)[[1]]
  if (is.na(union) || union == 0) return(NA_real_)
  inter / union
}

# Robust save: accept ggplot/patchwork OR a raw list of plots (wraps lists to 2x2)
save_plot_safe <- function(filename, plt, width, height, dpi = 300, bg = "white") {
  if (is.list(plt) && !inherits(plt, "ggplot") && !inherits(plt, "patchwork")) {
    plt <- patchwork::wrap_plots(plt, ncol = 2, byrow = TRUE)
  }
  ggplot2::ggsave(filename = filename, plot = plt, width = width, height = height, dpi = dpi, bg = bg)
}

# Align two rasters to the **same grid/CRS/extent**
align_to <- function(r1, r2, categorical = FALSE) {
  stopifnot(inherits(r1, "SpatRaster"), inherits(r2, "SpatRaster"))
  # project r2 to r1 CRS
  if (!terra::same.crs(r1, r2)) {
    r2 <- terra::project(r2, r1, method = if (categorical) "near" else "bilinear")
  }
  # intersect extent
  ex <- terra::intersect(ext(r1), ext(r2))
  r1c <- terra::crop(r1, ex)
  r2c <- terra::crop(r2, ex)
  # resample r2 to r1 grid
  r2a <- terra::resample(r2c, r1c, method = if (categorical) "near" else "bilinear")
  list(r1 = r1c, r2 = r2a)
}

# AOI to lines data.frame (x,y,grp) in raster CRS for simple overlay
aoi_to_df <- function(aoi_sf, template_rast) {
  if (is.null(aoi_sf)) return(NULL)
  aoi <- suppressWarnings(suppressMessages(sf::st_transform(aoi_sf, terra::crs(template_rast))))
  # Cast polygons to multilines for outline
  geom <- sf::st_geometry(aoi)
  types <- unique(sf::st_geometry_type(geom))
  if (any(grepl("POLYGON", types))) {
    geom <- sf::st_cast(geom, "MULTILINESTRING")
  }
  aoi2 <- sf::st_as_sf(geom)
  coords <- sf::st_coordinates(aoi2)
  df <- as.data.frame(coords)
  if (!"L1" %in% names(df)) df$L1 <- 1L
  df <- df |>
    dplyr::rename(x = X, y = Y) |>
    dplyr::mutate(grp = interaction(L1, L2, drop = TRUE))
  df[, c("x","y","grp")]
}

# ---------------------- METRIC + PLOT BUILDERS --------------------------------
# Compact one-line metrics strip (Jaccard@q and %Gain@q)
make_metrics_strip_time <- function(jac_vec, gain_pct_vec, qs) {
  jlab <- paste0("Jaccard@", sprintf("%.0f%%", qs*100), " ", sprintf("%.3f", jac_vec))
  alab <- paste0("%Gain@",  sprintf("%.0f%%", qs*100), " ", sprintf("%.1f", gain_pct_vec))
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste(c(jlab, alab), collapse = "   •   "),
             size = 4, hjust = 0.5, vjust = 0.5) +
    theme_void() +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE)
}

# Two-class GSL (red=Loss before-only, blue=Gain after-only); stable = NA
make_gsl2_plot_time <- function(b_before, b_after, title_txt, aoi_df = NULL) {
  cls <- (b_before & !b_after) * 1 + (b_after & !b_before) * 2
  df  <- as.data.frame(cls, xy = TRUE, na.rm = FALSE); names(df) <- c("x","y","cls")
  df$cls <- factor(df$cls, levels = c(1,2),
                   labels = c("Loss (Before only)", "Gain (After only)"))
  p <- ggplot(df, aes(x = x, y = y, fill = cls)) +
    geom_raster(na.rm = FALSE) + coord_equal() +
    scale_fill_manual(values = c("Loss (Before only)" = "#d7191c",
                                 "Gain (After only)"  = "#2c7fb8"),
                      name = "GSL", drop = TRUE) +
    labs(title = title_txt, x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())
  if (!is.null(aoi_df)) {
    p <- p + geom_path(data = aoi_df, aes(x = x, y = y, group = grp),
                       inherit.aes = FALSE, color = "black", linewidth = 0.4)
  }
  p
}

# Build a 3x3 Before/After panel for one framework (H2O or SSDM)
build_ba_panel <- function(rB, rA, framework, eid, aoi_df) {

  # Δ limits (symmetric)
  limD <- {
    dvals <- terra::values(rA - rB, mat = FALSE)
    rng <- range(dvals, na.rm = TRUE)
    max(abs(rng))
  }

  # ------------- Base maps (Before, After, Δ) -------------
  df_B <- as_xy(rB); df_A <- as_xy(rA); df_D <- as_xy(rA - rB)

  p_B <- ggplot(df_B, aes(x=x,y=y,fill=val)) +
    geom_raster() + coord_equal() +
    scale_fill_viridis_c(name = "Suitability") +
    labs(title = paste0(eid, " — ", framework, " (Before)")) +
    theme_minimal(base_size = 11) +
    theme(axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank())

  p_A <- ggplot(df_A, aes(x=x,y=y,fill=val)) +
    geom_raster() + coord_equal() +
    scale_fill_viridis_c(name = "Suitability") +
    labs(title = paste0(eid, " — ", framework, " (After)")) +
    theme_minimal(base_size = 11) +
    theme(axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank())

  p_D <- ggplot(df_D, aes(x=x,y=y,fill=val)) +
    geom_raster() + coord_equal() +
    scale_fill_gradientn(colors=c("#2c7fb8","#ffffbf","#d7191c"),
                         limits=c(-limD, limD), name=expression(Delta~"(A - B)")) +
    labs(title = expression(Delta~"(A - B)")) +
    theme_minimal(base_size = 11) +
    theme(axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank())

  # Optional AOI overlay (add to each plot explicitly)
  if (!is.null(aoi_df)) {
    p_B <- p_B + geom_path(data=aoi_df, aes(x=x,y=y,group=grp),
                           inherit.aes=FALSE, color="black", linewidth=0.4)
    p_A <- p_A + geom_path(data=aoi_df, aes(x=x,y=y,group=grp),
                           inherit.aes=FALSE, color="black", linewidth=0.4)
    p_D <- p_D + geom_path(data=aoi_df, aes(x=x,y=y,group=grp),
                           inherit.aes=FALSE, color="black", linewidth=0.4)
  }

  # ------------- GSL & Metrics -------------
  vB <- terra::values(rB, mat = FALSE); vA <- terra::values(rA, mat = FALSE)
  idx <- !is.na(vB) & !is.na(vA)

  # Handle degenerate case (no overlapping valid pixels)
  if (!any(idx)) {
    jacs <- rep(NA_real_, length(qs))
    gain_pct <- rep(NA_real_, length(qs))
    p_gsl <- replicate(length(qs), plot_spacer(), simplify = FALSE)
  } else {
    qB <- stats::quantile(vB[idx], probs = qs, na.rm = TRUE)
    qA <- stats::quantile(vA[idx], probs = qs, na.rm = TRUE)

    jacs <- gain_pct <- numeric(length(qs))
    p_gsl <- vector("list", length(qs))
    for (i in seq_along(qs)) {
      bB <- rB > qB[i]; bA <- rA > qA[i]
      jacs[i] <- jaccard_binary(bB, bA)
      n_valid <- global(!is.na(rB) & !is.na(rA), "sum", na.rm = TRUE)[[1]]
      n_gain  <- global(bA & !bB, "sum", na.rm = TRUE)[[1]]
      gain_pct[i] <- ifelse(n_valid > 0, 100 * n_gain / n_valid, NA_real_)
      p_gsl[[i]] <- make_gsl2_plot_time(bB, bA, paste0("GSL (thr=", qs[i], ")"), aoi_df)
    }
  }

  p_strip <- make_metrics_strip_time(jacs, gain_pct, qs)

  # ------------- Assemble 3×3 Layout -------------
  row1 <- wrap_elements(p_B) | wrap_elements(p_A) | wrap_elements(p_D)
  row2 <- (p_gsl[[1]] %||% plot_spacer()) |
          (p_gsl[[2]] %||% plot_spacer()) |
          (p_gsl[[3]] %||% plot_spacer())
  # pad metrics row to exactly 3 panels
  row3 <- wrap_elements(p_strip) | plot_spacer() | plot_spacer()

  panel <- (row1 / row2 / row3) +
    plot_layout(
      guides = "collect",
      widths = c(1, 1, 1),
      heights = c(1, 0.9, 0.28)
    ) +
    plot_annotation(
      title = paste0("Elephant ", eid, " — ", framework,
                     " | Before vs After: Δ (A−B), GSL, metrics")
    )

  panel & theme(
    legend.position   = "right",
    legend.box        = "vertical",
    legend.key.height = unit(10, "pt"),
    legend.key.width  = unit(10, "pt"),
    legend.title      = element_text(size = 10, face = "bold"),
    legend.text       = element_text(size = 9),
    plot.title        = element_text(face = "bold", hjust = 0)
  )
}

# ---------------------------- LOAD AOI (optional) -----------------------------
aoi_sf <- NULL
if (nzchar(aoi_path) && file.exists(aoi_path)) {
  aoi_sf <- tryCatch(sf::st_read(aoi_path, quiet = TRUE),
                     error = function(e) { message("[AOI] Failed to read: ", e$message); NULL })
  if (!is.null(aoi_sf)) message("[AOI] Loaded ", aoi_path)
} else {
  message("[AOI] Skipping AOI overlay (no path or file missing).")
}

# -------------------------- MAIN LOOP OVER E3/E4/E5 ---------------------------
out_dir <- file.path(compare_root, "04_panels_before_after")
dir_ensure(out_dir)

for (eid in ele_ids) {
  spB <- paste0(eid, "B"); spA <- paste0(eid, "A")

  # expected files
  f_B_h2o  <- file.path(h2o_root,  "B", spB, paste0("prediction_", spB, ".tif"))
  f_A_h2o  <- file.path(h2o_root,  "A", spA, paste0("prediction_", spA, ".tif"))
  f_B_ssdm <- file.path(ssdm_root, "B", spB, paste0("ESDM_",      spB, ".tif"))
  f_A_ssdm <- file.path(ssdm_root, "A", spA, paste0("ESDM_",      spA, ".tif"))

  have_h2o  <- file.exists(f_B_h2o)  && file.exists(f_A_h2o)
  have_ssdm <- file.exists(f_B_ssdm) && file.exists(f_A_ssdm)

  if (!have_h2o && !have_ssdm) {
    warning("[Before/After] Skipping ", eid, " — no complete H2O or SSDM pairs found.")
    next
  }

  # choose a raster to derive AOI CRS (prefer H2O A if available)
  aoi_df <- NULL
  if (!is.null(aoi_sf)) {
    if (have_h2o)  aoi_df <- aoi_to_df(aoi_sf, terra::rast(f_A_h2o))
    if (!have_h2o && have_ssdm) aoi_df <- aoi_to_df(aoi_sf, terra::rast(f_A_ssdm))
  }

  # H2O panel
  if (have_h2o) {
    r_B_h2o <- terra::rast(f_B_h2o); r_A_h2o <- terra::rast(f_A_h2o)
    al_h2o <- align_to(r_B_h2o, r_A_h2o, categorical = FALSE)
    p_h2o_panel <- build_ba_panel(al_h2o$r1, al_h2o$r2, "H2O", eid, aoi_df)
    out_dir_h2o <- file.path(out_dir, "H2O"); dir_ensure(out_dir_h2o)
    save_plot_safe(file.path(out_dir_h2o, paste0("panel_before_after_H2O_", eid, ".png")),
                   p_h2o_panel, width = 14.5, height = 9.5, dpi = 300, bg = "white")
  }

  # SSDM panel
  if (have_ssdm) {
    r_B_ssdm <- terra::rast(f_B_ssdm); r_A_ssdm <- terra::rast(f_A_ssdm)
    al_ssdm <- align_to(r_B_ssdm, r_A_ssdm, categorical = FALSE)
    p_ssdm_panel <- build_ba_panel(al_ssdm$r1, al_ssdm$r2, "SSDM", eid, aoi_df)
    out_dir_ssdm <- file.path(out_dir, "SSDM"); dir_ensure(out_dir_ssdm)
    save_plot_safe(file.path(out_dir_ssdm, paste0("panel_before_after_SSDM_", eid, ".png")),
                   p_ssdm_panel, width = 14.5, height = 9.5, dpi = 300, bg = "white")
  }

  # Combined (stack H2O over SSDM) if both exist
  if (have_h2o && have_ssdm) {
    combined <- (p_h2o_panel / p_ssdm_panel) +
      plot_layout(guides = "collect") +
      plot_annotation(title = paste0("Elephant ", eid, " — Before vs After (H2O & SSDM)"))
    combined <- combined & theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0)
    )
    out_dir_comb <- file.path(out_dir, "Combined"); dir_ensure(out_dir_comb)
    save_plot_safe(file.path(out_dir_comb, paste0("panel_before_after_Combined_", eid, ".png")),
                   combined, width = 14.5, height = 19.0, dpi = 300, bg = "white")
  }
}

message("[DONE] Panels written under: ", normalizePath(out_dir, winslash = "/"))
