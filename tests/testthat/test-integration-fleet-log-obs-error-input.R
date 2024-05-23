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
    log_obs_error_env$fims$fishing_fleet$log_obs_error <- log_obs_error
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

test_that("wrong fleet log_obs_error input dimenstion triggers failure", {
  load(test_path("fixtures", "integration_test_data.RData"))

  # true nyr = 30, set up log_obs_error input values with wrong dimension
  log_obs_error_long <- rep(log(sqrt(log(em_input$cv.L$fleet1^2 + 1))), 100)

  for (i in seq_along(log_obs_error_long)){
    log_obs_error_long_env <- setup_fims(
      om_input = om_input,
      om_output = om_output,
      em_input = em_input
    )

    # Change fleet log_obs_error from a vector of values to a single value
    log_obs_error_long_env$fims$fishing_fleet$log_obs_error <- log_obs_error_long
    # Set-up TMB
    log_obs_error_long_env$fims$CreateTMBModel()
    # Create parameter list from Rcpp modules
    parameters <- list(p = log_obs_error_long_env$fims$get_fixed())
    obj <- TMB::MakeADFun(data = list(), parameters, DLL = "FIMS")

    opt <- with(obj, optim(par, fn, gr,
                           method = "BFGS",
                           control = list(maxit = 1000000, reltol = 1e-15)
    ))

    validate_integration_test_results(tmb_obj = obj)

    log_obs_error_long_env$fims$clear()
  }

  # true nyr = 30, set up log_obs_error input values with wrong dimension
  log_obs_error_short <- rep(log(sqrt(log(em_input$cv.L$fleet1^2 + 1))), 10)

  for (i in seq_along(log_obs_error_short)){
    log_obs_error_short_env <- setup_fims(
      om_input = om_input,
      om_output = om_output,
      em_input = em_input
    )

    # Change fleet log_obs_error from a vector of values to a single value
    log_obs_error_short_env$fims$fishing_fleet$log_obs_error <- log_obs_error_short
    # Set-up TMB
    log_obs_error_short_env$fims$CreateTMBModel()
    # Create parameter list from Rcpp modules
    parameters <- list(p = log_obs_error_short_env$fims$get_fixed())
    obj <- TMB::MakeADFun(data = list(), parameters, DLL = "FIMS")

    opt <- with(obj, optim(par, fn, gr,
                           method = "BFGS",
                           control = list(maxit = 1000000, reltol = 1e-15)
    ))

    validate_integration_test_results(tmb_obj = obj)

    log_obs_error_short_env$fims$clear()
  }

})
