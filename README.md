# ems_paper
# Run SSDM for all species in both runs A and B
for (RUN in c("A","B")) {
  system2(file.path(R.home("bin"), "Rscript"),
    c("scripts/04_ssdm_train.R",
      "--run", RUN,
      "--mode", "FAST")
  )
}
# Run H2O AutoML for all species in both runs (A and B)
for (RUN in c("A","B")) {
  system2(file.path(R.home("bin"), "Rscript"),
    c("scripts/03_h2o_train.R",
      "--run", RUN,
      "--mode", "FAST",
      "--max_models", "40",
      "--max_runtime_secs", "300")
  )
}
