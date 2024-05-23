# Set-up OM (sigmaR = 0.4)
remotes::install_github(
  repo = "Bai-Li-NOAA/Age_Structured_Stock_Assessment_Model_Comparison"
)

working_dir <- getwd()

maindir <- tempdir()
model_input <- ASSAMC::save_initial_input()
FIMS_C0_estimation <- ASSAMC::save_initial_input(
  base_case = TRUE,
  input_list = model_input,
  maindir = maindir,
  om_sim_num = 1,
  keep_sim_num = 1,
  figure_number = 1,
  seed_num = 9924,
  case_name = "FIMS_C0_estimation"
)

# generate om_input, om_output, and em_input
# using function from the model comparison project
ASSAMC::run_om(input_list = FIMS_C0_estimation)

on.exit(unlink(maindir, recursive = TRUE), add = TRUE)

setwd(working_dir)
on.exit(setwd(working_dir), add = TRUE)

save(om_input, om_output, em_input,
     file = test_path("fixtures", "integration_test_data.RData"))
