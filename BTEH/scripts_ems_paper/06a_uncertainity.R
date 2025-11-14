#!/usr/bin/env Rscript
# =============================================================================
# 06a_uncertainty.R (folder-complete outputs)
# Uncertainty Decomposition (SSF–RSF vs H2O & SSDM), per elephant & run (A/B)
#
# Inputs (tries in order, first hit wins)
#   SP = <ELE><RUN>  e.g., E3A, E5B
#   H2O  : results/H2O/<RUN>/<SP>/prediction_<SP>.tif
#          results/H2O/<RUN>/prediction_<SP>.tif
#          (fallback: first *.tif under those dirs)
#   SSDM : results/SSDM/<RUN>/<SP>/ESDM_<SP>.tif
#          results/SSDM/<RUN>/ESDM_<SP>.tif
#          (fallback: first *.tif under those dirs)
#   SSF  : results/SSF/<RUN>/<SP>/<SP>_SSF_rsf_0to1.tif
#          results/SSF/<RUN>/<SP>_SSF_rsf_0to1.tif
#          results/SSF/<SP>_SSF_rsf_0to1.tif
#          (fallback: first *.tif under those dirs)
#
# Outputs:
#   paper_results/uncertainty/<ELE>/
#     H2O/      <ELE>_A/B_H2O.tif, <ELE>_A/B_H2O_q75mask.tif, <ELE>_delta_BminusA.tif
#     SSDM/     <ELE>_A/B_SSDM.tif, <ELE>_A/B_SSDM_q75mask.tif, <ELE>_delta_BminusA.tif
#     SSF/      <ELE>_A/B_SSF.tif,  <ELE>_A/B_SSF_q75mask.tif,  <ELE>_delta_BminusA.tif
#     combined/ mean/sd/cv/agreement/ diffs/ change classes (TIF)
#     figures/  ALL PNGs for maps
#     tables/   per-run stats, pairwise corrs, Jaccard, temporal Jaccard
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(viridisLite)
  library(rnaturalearth)
})

# ------------------------------- CLI -----------------------------------------
option_list <- list(
  make_option(c("-e","--elephant"), type="character", default=NULL,
              help="Elephants: E3 or E3,E4,E5 or ALL. (default set used if source() without args)"),
  make_option(c("-a","--tagA"), type="character", default="A", help="Before tag [default %default]"),
  make_option(c("-b","--tagB"), type="character", default="B", help="After tag [default %default]"),
  make_option(c("-r","--root"), type="character", default="results", help="Results root [default %default]"),
  make_option(c("-q","--qmask"), type="double", default=0.75, help="Quantile for hotspot masks [default %default]"),
  make_option(c("-o","--outdir"), type="character", default="paper_results/uncertainty", help="Output root [default %default]"),
  make_option(c("-s","--seed"), type="integer", default=20161113L, help="Seed [default %default]"),
  make_option(c("--skip_figs"), action="store_true", default=FALSE, help="Skip PNG rendering")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ---- Coastline toggles ----
opt$coast <- TRUE
if (exists("opt_disable_coast", inherits = TRUE) && isTRUE(opt_disable_coast)) opt$coast <- FALSE
if (!requireNamespace("rnaturalearth", quietly = TRUE) ||
    (!requireNamespace("rnaturalearthdata", quietly = TRUE) &&
     !requireNamespace("rnaturalearthhires", quietly = TRUE))) opt$coast <- FALSE

# Honor absolute root when sourcing()
if (exists("opt_root", inherits = TRUE) && nzchar(opt_root)) {
  opt$root <- opt_root
  message("[ROOT] Using opt_root = ", normalizePath(opt$root, mustWork = FALSE))
}

if (is.null(opt$elephant)) {
  message("No --elephant provided; assuming interactive run via source().")
  opt$elephant <- "E3,E4,E5,E1,E2,E6"
}

set.seed(opt$seed)
runA <- toupper(opt$tagA); runB <- toupper(opt$tagB)

# --------------------------- Path builders -----------------------------------
first_existing <- function(...) {
  cands <- c(...)
  for (p in cands) if (nzchar(p) && file.exists(p)) return(p)
  cands[1]
}
sp_tag <- function(ele, run) paste0(ele, run)

p_h2o <- function(ele, run) {
  SP <- sp_tag(ele, run)
  d1 <- file.path(opt$root, "H2O", "A", SP)
  d2 <- file.path(opt$root, "H2O", "A")
  first_existing(
    file.path(d1, sprintf("prediction_%s.tif", SP)),
    file.path(d2, sprintf("prediction_%s.tif", SP)),
    { tf <- if (dir.exists(d1)) list.files(d1, "\\.tif$", full.names=TRUE); if (length(tf)) tf[1] else "" },
    { tf <- if (dir.exists(d2)) list.files(d2, "\\.tif$", full.names=TRUE); if (length(tf)) tf[1] else "" }
  )
}
p_ssdm <- function(ele, run) {
  SP <- sp_tag(ele, run)
  d1 <- file.path(opt$root, "SSDM", "A", SP)
  d2 <- file.path(opt$root, "SSDM", "A")
  first_existing(
    file.path(d1, sprintf("ESDM_%s.tif", SP)),
    file.path(d2, sprintf("ESDM_%s.tif", SP)),
    { tf <- if (dir.exists(d1)) list.files(d1, "\\.tif$", full.names=TRUE); if (length(tf)) tf[1] else "" },
    { tf <- if (dir.exists(d2)) list.files(d2, "\\.tif$", full.names=TRUE); if (length(tf)) tf[1] else "" }
  )
}
p_ssf <- function(ele, run) {
  SP <- sp_tag(ele, run)
  d1 <- file.path(opt$root, "SSF", run, SP)
  d2 <- file.path(opt$root, "SSF", run)
  d3 <- file.path(opt$root, "SSF")
  first_existing(
    file.path(d1, sprintf("%s_SSF_rsf_0to1.tif", SP)),
    file.path(d2, sprintf("%s_SSF_rsf_0to1.tif", SP)),
    file.path(d3, sprintf("%s_SSF_rsf_0to1.tif", SP)),
    { tf <- if (dir.exists(d1)) list.files(d1, "\\.tif$", full.names=TRUE); if (length(tf)) tf[1] else "" },
    { tf <- if (dir.exists(d2)) list.files(d2, "\\.tif$", full.names=TRUE); if (length(tf)) tf[1] else "" },
    { tf <- if (dir.exists(d3)) list.files(d3, "\\.tif$", full.names=TRUE); if (length(tf)) tf[1] else "" }
  )
}
exists_all <- function(ele, run) {
  all(file.exists(p_h2o(ele, run)),
      file.exists(p_ssdm(ele, run)),
      file.exists(p_ssf(ele, run)))
}

# --------------------------- Output dirs -------------------------------------
mk_outdirs <- function(ele) {
  out_root <- file.path(opt$outdir, ele)
  dirs <- list(
    H2O=file.path(out_root,"H2O"),
    SSDM=file.path(out_root,"SSDM"),
    SSF=file.path(out_root,"SSF"),
    combined=file.path(out_root,"combined"),
    figures=file.path(out_root,"figures"),
    tables=file.path(out_root,"tables")
  )
  invisible(lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE))
  dirs
}

# --------------------------- Plot helpers ------------------------------------
.coast_cache <- new.env(parent = emptyenv())
get_coastline <- function(template_rast){
  if (!isTRUE(opt$coast)) return(NULL)
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) return(NULL)
  if (!requireNamespace("rnaturalearthdata", quietly = TRUE) &&
      !requireNamespace("rnaturalearthhires", quietly = TRUE)) return(NULL)
  key <- paste0("crs:", as.character(terra::crs(template_rast)))
  if (!exists(key, envir = .coast_cache)) {
    coast_sf <- tryCatch(rnaturalearth::ne_coastline(scale="medium", returnclass="sf"), error=function(e) NULL)
    if (is.null(coast_sf)) return(NULL)
    coast_v <- terra::vect(coast_sf)
    if (!is.na(terra::crs(template_rast))) {
      coast_v <- tryCatch(terra::project(coast_v, template_rast), error=function(e) coast_v)
    }
    assign(key, coast_v, envir = .coast_cache)
  }
  get(key, envir = .coast_cache)
}
save_raster_map <- function(r, title, fname,
                            palette=c("sequential","diverging","binary","tri"),
                            width_px=1600, height_px=1200, dpi=240){
  palette <- match.arg(palette)
  dir.create(dirname(fname), TRUE, TRUE)
  v <- values(r); v <- v[is.finite(v)]
  png(fname, width=width_px, height=height_px, res=dpi, bg="white")
  op <- par(no.readonly=TRUE); on.exit({par(op); dev.off()}, add=TRUE)
  par(mar=c(3.2,3.6,3.2,6))
  coast <- get_coastline(r)
  if (length(v)==0){ plot(r, main=title); if(!is.null(coast)) try(lines(crop(coast, ext(r)), lwd=0.4), silent=TRUE); return() }
  if (palette=="binary"){
    plot(r, main=title, col=c("#f0f0f0","#3b528b"),
         breaks=c(-0.5,0.5,1.5), axes=FALSE, box=TRUE, legend=TRUE, colNA="grey85")
  } else if (palette=="diverging"){
    vmax <- max(abs(range(v, na.rm=TRUE)))
    plot(r, main=title, col=colorRampPalette(c("#3b528b","#f7f7f7","#b40426"))(201),
         zlim=c(-vmax, vmax), axes=FALSE, box=TRUE, legend=TRUE, colNA="grey85")
  } else if (palette=="tri"){
    plot(r, main=title, col=c("#b40426","#f7f7f7","#3b528b"),
         breaks=c(0.5,1.5,2.5,3.5), axes=FALSE, box=TRUE, legend=TRUE, colNA="grey85")
  } else {
    plot(r, main=title, col=viridisLite::viridis(200),
         zlim=range(v, na.rm=TRUE), axes=FALSE, box=TRUE, legend=TRUE, colNA="grey85")
  }
  if(!is.null(coast)) try(lines(crop(coast, ext(r)), lwd=0.4, col="black"), silent=TRUE)
}

# ------------------------------- Utilities -----------------------------------
write_gt <- function(r, p){ dir.create(dirname(p), TRUE, TRUE); writeRaster(r, p, overwrite=TRUE) }
clamp01  <- function(x){ x <- ifel(x<0,0,x); ifel(x>1,1,x) }
qmask <- function(r, q = 0.75) {
  thr <- tryCatch(
    as.numeric(terra::global(r, fun = quantile, probs = q, na.rm = TRUE)[1, 1]),
    error = function(e) NA_real_
  )
  if (is.na(thr)) {
    v <- terra::values(r, mat = FALSE); v <- v[is.finite(v)]
    if (length(v) == 0L) return(r * NA)
    thr <- as.numeric(stats::quantile(v, probs = q, na.rm = TRUE, names = FALSE))
  }
  r > thr
}
rmse <- function(a,b) sqrt(mean((a-b)^2, na.rm=TRUE))
mae  <- function(a,b) mean(abs(a-b), na.rm=TRUE)

# Align SSF/H2O/SSDM to a common finest grid, clamp to [0,1]
harmonize <- function(SSF, H2O, SSDM){
  resos  <- c(SSF=mean(res(SSF)), H2O=mean(res(H2O)), SSDM=mean(res(SSDM)))
  target <- list(SSF=SSF, H2O=H2O, SSDM=SSDM)[[names(which.min(resos))]]
  conv <- function(r) {
    if (!compareGeom(r, target, stopOnError=FALSE)) {
      r <- project(r, crs(target))
      r <- resample(r, target, method="bilinear")
    }
    r
  }
  stk <- c(conv(SSF), conv(H2O), conv(SSDM)); names(stk) <- c("SSF","H2O","SSDM")
  clamp01(stk)
}

# ------------------------------- Per-run -------------------------------------
analyze_run <- function(ele, run, dirs){
  SP <- sp_tag(ele, run)
  message(sprintf("Reading %s (%s)…", ele, run))

  f_ssf  <- p_ssf(ele, run)
  f_h2o  <- p_h2o(ele, run)
  f_ssdm <- p_ssdm(ele, run)

  message(sprintf("[INPUT] SSF : %s (exists=%s)", normalizePath(f_ssf,  mustWork=FALSE), file.exists(f_ssf)))
  message(sprintf("[INPUT] H2O : %s (exists=%s)", normalizePath(f_h2o,  mustWork=FALSE), file.exists(f_h2o)))
  message(sprintf("[INPUT] SSDM: %s (exists=%s)", normalizePath(f_ssdm, mustWork=FALSE), file.exists(f_ssdm)))
  if (!all(file.exists(c(f_ssf, f_h2o, f_ssdm)))) {
    stop(sprintf("Missing inputs for %s (%s). Check paths above.", ele, run))
  }

  r_ssf  <- rast(f_ssf)
  r_h2o  <- rast(f_h2o)
  r_ssdm <- rast(f_ssdm)

  stk <- harmonize(r_ssf, r_h2o, r_ssdm)

  # ---- per-method outputs (rasters go into method folders) ----
  write_gt(stk[["H2O"]],  file.path(dirs$H2O,  sprintf("%s_%s_H2O.tif",  ele, run)))
  write_gt(stk[["SSDM"]], file.path(dirs$SSDM, sprintf("%s_%s_SSDM.tif", ele, run)))
  write_gt(stk[["SSF"]],  file.path(dirs$SSF,  sprintf("%s_%s_SSF.tif",  ele, run)))

  m_ssf  <- qmask(stk[["SSF"]],  opt$qmask)
  m_h2o  <- qmask(stk[["H2O"]],  opt$qmask)
  m_ssdm <- qmask(stk[["SSDM"]], opt$qmask)
  write_gt(m_h2o,  file.path(dirs$H2O,  sprintf("%s_%s_H2O_q%02dmask.tif",  ele, run, round(opt$qmask*100))))
  write_gt(m_ssdm, file.path(dirs$SSDM, sprintf("%s_%s_SSDM_q%02dmask.tif", ele, run, round(opt$qmask*100))))
  write_gt(m_ssf,  file.path(dirs$SSF,  sprintf("%s_%s_SSF_q%02dmask.tif",  ele, run, round(opt$qmask*100))))

  # ---- combined summaries (rasters into combined/) ----
  r_mean <- mean(stk, na.rm=TRUE); names(r_mean) <- "mean"
  r_sd   <- stdev(stk, na.rm=TRUE); names(r_sd)   <- "sd"
  r_cv   <- r_sd / (r_mean + 1e-9); names(r_cv)   <- "cv"
  r_alg  <- stk[["H2O"]] - stk[["SSDM"]]; names(r_alg) <- "diff_H2O_minus_SSDM"
  r_sdmM <- (stk[["H2O"]] + stk[["SSDM"]]) / 2
  r_beh  <- r_sdmM - stk[["SSF"]]; names(r_beh) <- "diff_SDMmean_minus_SSF"
  m_ag   <- (m_ssf + m_h2o + m_ssdm) / 3; names(m_ag) <- "agree_prop"

  write_gt(r_mean, file.path(dirs$combined, sprintf("%s_%s_mean.tif", ele, run)))
  write_gt(r_sd,   file.path(dirs$combined, sprintf("%s_%s_sd.tif",   ele, run)))
  write_gt(r_cv,   file.path(dirs$combined, sprintf("%s_%s_cv.tif",   ele, run)))
  write_gt(r_alg,  file.path(dirs$combined, sprintf("%s_%s_diff_H2O_minus_SSDM.tif", ele, run)))
  write_gt(r_beh,  file.path(dirs$combined, sprintf("%s_%s_diff_SDMmean_minus_SSF.tif", ele, run)))
  write_gt(m_ag,   file.path(dirs$combined, sprintf("%s_%s_agreement_prop_q%02d.tif", ele, run, round(opt$qmask*100))))

  # ---- PNGs go to figures/ ----
  if (!opt$skip_figs){
    # per-method surfaces
    save_raster_map(stk[["H2O"]],  sprintf("%s %s — H2O",  ele, run),
                    file.path(dirs$figures, sprintf("%s_%s_H2O.png",  ele, run)))
    save_raster_map(stk[["SSDM"]], sprintf("%s %s — SSDM", ele, run),
                    file.path(dirs$figures, sprintf("%s_%s_SSDM.png", ele, run)))
    save_raster_map(stk[["SSF"]],  sprintf("%s %s — SSF",  ele, run),
                    file.path(dirs$figures, sprintf("%s_%s_SSF.png",  ele, run)))
    # masks
    save_raster_map(m_h2o,  sprintf("%s %s — H2O Hotspots (Q%02d)",  ele, run, round(opt$qmask*100)),
                    file.path(dirs$figures, sprintf("%s_%s_H2O_q%02dmask.png",  ele, run, round(opt$qmask*100))), "binary")
    save_raster_map(m_ssdm, sprintf("%s %s — SSDM Hotspots (Q%02d)", ele, run, round(opt$qmask*100)),
                    file.path(dirs$figures, sprintf("%s_%s_SSDM_q%02dmask.png", ele, run, round(opt$qmask*100))), "binary")
    save_raster_map(m_ssf,  sprintf("%s %s — SSF Hotspots (Q%02d)",  ele, run, round(opt$qmask*100)),
                    file.path(dirs$figures, sprintf("%s_%s_SSF_q%02dmask.png",  ele, run, round(opt$qmask*100))), "binary")
    # combined summaries
    save_raster_map(r_mean, sprintf("%s %s — Mean", ele, run),
                    file.path(dirs$figures, sprintf("%s_%s_mean.png", ele, run)))
    save_raster_map(r_sd,   sprintf("%s %s — SD",   ele, run),
                    file.path(dirs$figures, sprintf("%s_%s_sd.png",   ele, run)))
    save_raster_map(r_cv,   sprintf("%s %s — CV",   ele, run),
                    file.path(dirs$figures, sprintf("%s_%s_cv.png",   ele, run)))
    save_raster_map(r_alg,  sprintf("%s %s — H2O − SSDM", ele, run),
                    file.path(dirs$figures, sprintf("%s_%s_diff_H2O_minus_SSDM.png", ele, run)), "diverging")
    save_raster_map(r_beh,  sprintf("%s %s — SDM̄ − SSF", ele, run),
                    file.path(dirs$figures, sprintf("%s_%s_diff_SDMmean_minus_SSF.png", ele, run)), "diverging")
    save_raster_map(m_ag,   sprintf("%s %s — Agreement (Q%02d)", ele, run, round(opt$qmask*100)),
                    file.path(dirs$figures, sprintf("%s_%s_agreement_q%02d.png", ele, run, round(opt$qmask*100))))
  }

  # ---- tables ----
  df <- as.data.frame(stk, xy=FALSE, na.rm=TRUE)
  glb <- tibble(
    elephant = ele, run = run, metric = c("mean","sd","cv"),
    SSF = c(global(stk[["SSF"]],"mean",na.rm=TRUE)[1,1],
            global(stk[["SSF"]],"sd",  na.rm=TRUE)[1,1],
            global(stk[["SSF"]],"sd",  na.rm=TRUE)[1,1]/(global(stk[["SSF"]],"mean",na.rm=TRUE)[1,1]+1e-9)),
    H2O = c(global(stk[["H2O"]],"mean",na.rm=TRUE)[1,1],
            global(stk[["H2O"]],"sd",  na.rm=TRUE)[1,1],
            global(stk[["H2O"]],"sd",  na.rm=TRUE)[1,1]/(global(stk[["H2O"]],"mean",na.rm=TRUE)[1,1]+1e-9)),
    SSDM= c(global(stk[["SSDM"]],"mean",na.rm=TRUE)[1,1],
            global(stk[["SSDM"]],"sd",  na.rm=TRUE)[1,1],
            global(stk[["SSDM"]],"sd",  na.rm=TRUE)[1,1]/(global(stk[["SSDM"]],"mean",na.rm=TRUE)[1,1]+1e-9))
  )

  cor_tbl <- tibble(
    elephant = ele, run = run,
    pair     = c("H2O~SSDM","SSF~H2O","SSF~SSDM","SSF~SDMmean"),
    pearson  = c(cor(df$H2O, df$SSDM, use="pairwise"),
                 cor(df$SSF, df$H2O,  use="pairwise"),
                 cor(df$SSF, df$SSDM, use="pairwise"),
                 cor(df$SSF, (df$H2O+df$SSDM)/2, use="pairwise")),
    spearman = c(cor(df$H2O, df$SSDM, use="pairwise", method="spearman"),
                 cor(df$SSF, df$H2O,  use="pairwise", method="spearman"),
                 cor(df$SSF, df$SSDM, use="pairwise", method="spearman"),
                 cor(df$SSF, (df$H2O+df$SSDM)/2, use="pairwise", method="spearman")),
    rmse     = c(rmse(df$H2O, df$SSDM),
                 rmse(df$SSF, df$H2O),
                 rmse(df$SSF, df$SSDM),
                 rmse(df$SSF, (df$H2O+df$SSDM)/2)),
    mae      = c(mae(df$H2O, df$SSDM),
                 mae(df$SSF, df$H2O),
                 mae(df$SSF, df$SSDM),
                 mae(df$SSF, (df$H2O+df$SSDM)/2))
  )

  # Pairwise Jaccard@Q (numbers for tables)
  jac <- {
    inter <- global((m_ssf & m_h2o), "sum", na.rm=TRUE)[1,1]; uni <- global((m_ssf | m_h2o), "sum", na.rm=TRUE)[1,1]
    j1 <- ifelse(uni==0, NA_real_, inter/uni)
    inter <- global((m_ssf & m_ssdm), "sum", na.rm=TRUE)[1,1]; uni <- global((m_ssf | m_ssdm), "sum", na.rm=TRUE)[1,1]
    j2 <- ifelse(uni==0, NA_real_, inter/uni)
    inter <- global((m_h2o & m_ssdm), "sum", na.rm=TRUE)[1,1]; uni <- global((m_h2o | m_ssdm), "sum", na.rm=TRUE)[1,1]
    j3 <- ifelse(uni==0, NA_real_, inter/uni)
    tibble(elephant=ele, run=run, q=opt$qmask,
           pair=c("SSF~H2O","SSF~SSDM","H2O~SSDM"),
           jaccard_q=c(j1,j2,j3))
  }

  list(stack=stk, summaries=glb, cors=cor_tbl, jacc=jac,
       masks=list(H2O=m_h2o, SSDM=m_ssdm, SSF=m_ssf))
}

# ----------------------------- Temporal A vs B --------------------------------
temporal_compare <- function(ele, Astk, Bstk, dirs){
  methods <- c("SSF","H2O","SSDM")
  jac <- list()
  for (m in methods){
    rA <- Astk[[m]]; rB <- Bstk[[m]]
    d  <- rB - rA
    # raster to method dir
    write_gt(d, file.path(dirs[[m]], sprintf("%s_delta_BminusA.tif", ele)))
    # png to figures
    if (!opt$skip_figs)
      save_raster_map(d, sprintf("%s Δ (B−A) — %s", ele, m),
                      file.path(dirs$figures, sprintf("%s_delta_BminusA_%s.png", ele, m)), "diverging")

    # hotspot change classes
    mA <- qmask(rA, opt$qmask); mB <- qmask(rB, opt$qmask)
    inter <- global(mA & mB, "sum", na.rm=TRUE)[1,1]
    uni   <- global(mA | mB, "sum", na.rm=TRUE)[1,1]
    jac[[m]] <- tibble(elephant=ele, method=m, q=opt$qmask, jaccard_AvB=ifelse(uni==0, NA_real_, inter/uni))

    loss   <- (mA==1) & (mB==0)
    stable <- (mA==1) & (mB==1)
    gain   <- (mA==0) & (mB==1)
    change <- (loss*1) + (stable*2) + (gain*3)
    names(change) <- paste0("change_", m, "_Q", round(opt$qmask*100))
    write_gt(change, file.path(dirs$combined, sprintf("%s_change_%s_q%02d.tif", ele, m, round(opt$qmask*100))))
    if (!opt$skip_figs)
      save_raster_map(change, sprintf("%s Change — %s (Q%02d)", ele, m, round(opt$qmask*100)),
                      file.path(dirs$figures, sprintf("%s_change_%s_q%02d.png", ele, m, round(opt$qmask*100))), "tri")
  }
  bind_rows(jac)
}

# -------------------------- DEBUG: show resolved paths ------------------------
debug_inputs <- function(ele, run) {
  h <- p_h2o(ele, run); s <- p_ssdm(ele, run); f <- p_ssf(ele, run)
  cat("\n[DEBUG]", ele, run,
      "\n H2O :",  normalizePath(h, mustWork = FALSE), "exists?", file.exists(h),
      "\n SSDM:",  normalizePath(s, mustWork = FALSE), "exists?", file.exists(s),
      "\n SSF :",  normalizePath(f, mustWork = FALSE), "exists?", file.exists(f), "\n")
}

# --------------------------------- MAIN ---------------------------------------
elephants <- toupper(gsub("[^A-Z0-9,]", "", opt$elephant))
if (identical(elephants, "ALL")) elephants <- c("E1","E2","E3","E4","E5","E6") else elephants <- strsplit(elephants, ",", fixed = TRUE)[[1]]
elephants <- unique(elephants[nzchar(elephants)])

# strip accidental trailing run letters in inputs like E3A
normalize_ele <- function(s) {
  s <- trimws(s)
  if (nchar(s) >= 2) {
    last <- substr(s, nchar(s), nchar(s))
    if (last %in% c("A","B")) return(substr(s, 1, nchar(s)-1))
  }
  s
}
elephants <- unique(vapply(elephants, normalize_ele, FUN.VALUE = character(1)))

for (ele in elephants) {
  dirs <- mk_outdirs(ele)

  message(sprintf("\n== Uncertainty decomposition: %s ==", ele))
  debug_inputs(ele, runA)
  debug_inputs(ele, runB)

  hasA <- exists_all(ele, runA)
  hasB <- exists_all(ele, runB)

  if (!hasA && !hasB) {
    warning(sprintf("Skipping %s: neither %s nor %s inputs found.", ele, runA, runB))
    next
  }

  A <- NULL; B <- NULL
  if (hasA) {
    A <- analyze_run(ele, runA, dirs)
  } else message(sprintf("No %s run for %s — skipping A.", runA, ele))

  if (hasB) {
    B <- analyze_run(ele, runB, dirs)
  } else message(sprintf("No %s run for %s — skipping B.", runB, ele))

  # ----- Tables per run -----
  tbl_dir <- file.path(opt$outdir, ele, "tables"); dir.create(tbl_dir, TRUE, TRUE)
  if (!is.null(A)) {
    write_csv(A$summaries, file.path(tbl_dir, sprintf("%s_global_stats_%s.csv", ele, runA)))
    write_csv(A$cors,      file.path(tbl_dir, sprintf("%s_pairwise_correlations_%s.csv", ele, runA)))
    write_csv(A$jacc,      file.path(tbl_dir, sprintf("%s_jaccard_q%02d_%s.csv", ele, round(opt$qmask*100), runA)))
  }
  if (!is.null(B)) {
    write_csv(B$summaries, file.path(tbl_dir, sprintf("%s_global_stats_%s.csv", ele, runB)))
    write_csv(B$cors,      file.path(tbl_dir, sprintf("%s_pairwise_correlations_%s.csv", ele, runB)))
    write_csv(B$jacc,      file.path(tbl_dir, sprintf("%s_jaccard_q%02d_%s.csv", ele, round(opt$qmask*100), runB)))
  }

  # ----- Temporal only if both runs exist -----
  if (!is.null(A) && !is.null(B)) {
    TEMP <- temporal_compare(ele, A$stack, B$stack, dirs)
    write_csv(TEMP, file.path(tbl_dir, sprintf("%s_temporal_jaccard_q%02d.csv", ele, round(opt$qmask*100))))
  } else {
    message(sprintf("Temporal skipped for %s (need both %s and %s).", ele, runA, runB))
  }

  message(sprintf("Finished %s.", ele))
}

message("\nAll done.\n")
