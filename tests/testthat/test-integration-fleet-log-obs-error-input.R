test_that("single fleet log_obs_error input works", {
  load(test_path("fixtures", "integration_test_data.RData"))
  log_obs_error <- c(log(sqrt(log(0.05 ^ 2 + 1))), -2.9, 0.5,
                     1, -1, -1.1)

  for (i in seq_along(log_obs_error)){
    log_obs_error_env <- setup_fims(
      om_input = om_input,
      om_output = om_output,
      em_input = em_input
    )

    # Change fleet log_obs_error from a vector of values to a single value
    log_obs_error_env$fims$fishing_fleet$log_obs_error <- log_obs_error[i]
    # Set-up TMB
    log_obs_error_env$fims$CreateTMBModel()
    # Create parameter list from Rcpp modules
    parameters <- list(p = log_obs_error_env$fims$get_fixed())
    obj <- TMB::MakeADFun(data = list(), parameters, DLL = "FIMS")

    opt <- with(obj, optim(par, fn, gr,
                           method = "BFGS",
                           control = list(maxit = 1000000, reltol = 1e-15)
    ))

    validate_integration_test_results(tmb_obj = obj)

    log_obs_error_env$fims$clear()
  }

})
