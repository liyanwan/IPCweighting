folder_path = "~/IPCweighting/ML_SurvSL/R/"
r_files = list.files(folder_path, pattern = "\\.R$", full.names = TRUE)
sapply(r_files, source)
source("~/IPCweighting/ML_SurvSL/run_simulation_SL.R")

