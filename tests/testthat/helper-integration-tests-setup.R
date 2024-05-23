# Set-up Rcpp modules and fix parameters to "true"
# om_input, om_output, and em_input generated from ASSAMC::run_om
setup_fims <- function(om_input, om_output, em_input) {
  # create new R environment in which to conduct tests
  test_env <- new.env(parent = emptyenv())
  load(test_path("fixtures", "integration_test_data.RData"))
  # initialize FIMS model
  test_env$fims <- Rcpp::Module("fims", PACKAGE = "FIMS")

  # Recruitment
  # create new module in the recruitment class (specifically Beverton-Holt,
  # when there are other options, this would be where the option would be chosen)
  test_env$recruitment <- methods::new(test_env$fims$BevertonHoltRecruitment)

  # NOTE: in first set of parameters below (for recruitment),
  # $is_random_effect (default is FALSE) and $estimated (default is FALSE)
  # are defined even if they match the defaults in order to provide an example
  # of how that is done. Other sections of the code below leave defaults in
  # place as appropriate.

  # set up logR_sd
  # logR_sd is NOT logged. It needs to enter the model logged b/c the exp() is
  # taken before the likelihood calculation
  test_env$recruitment$log_sigma_recruit$value <- log(om_input$logR_sd)
  test_env$recruitment$log_sigma_recruit$is_random_effect <- FALSE
  test_env$recruitment$log_sigma_recruit$estimated <- FALSE
  # set up log_rzero (equilibrium recruitment)
  test_env$recruitment$log_rzero$value <- log(om_input$R0)
  test_env$recruitment$log_rzero$is_random_effect <- FALSE
  test_env$recruitment$log_rzero$estimated <- TRUE
  # set up logit_steep
  test_env$recruitment$logit_steep$value <- -log(1.0 - om_input$h) + log(om_input$h - 0.2)
  test_env$recruitment$logit_steep$is_random_effect <- FALSE
  test_env$recruitment$logit_steep$estimated <- FALSE
  # turn on estimation of deviations
  test_env$recruitment$estimate_log_devs <- TRUE
  # recruit deviations should enter the model in normal space.
  # The log is taken in the likelihood calculations
  # alternative setting: recruitment$log_devs <- rep(0, length(om_input$logR.resid))
  test_env$recruitment$log_devs <- om_input$logR.resid[-1]


  # Data
  test_env$catch <- em_input$L.obs$fleet1
  # set fishing fleet catch data, need to set dimensions of data index
  # currently FIMS only has a fleet module that takes index for both survey index and fishery catch
  test_env$fishing_fleet_index <- new(test_env$fims$Index, length(test_env$catch))
  test_env$fishing_fleet_index$index_data <- test_env$catch
  # set fishing fleet age comp data, need to set dimensions of age comps
  test_env$fishing_fleet_age_comp <- new(test_env$fims$AgeComp, length(test_env$catch), om_input$nages)
  test_env$fishing_fleet_age_comp$age_comp_data <- c(t(em_input$L.age.obs$fleet1)) * em_input$n.L$fleet1

  # repeat for surveys
  test_env$survey_index <- em_input$surveyB.obs$survey1
  test_env$survey_fleet_index <- new(test_env$fims$Index, length(test_env$survey_index))
  test_env$survey_fleet_index$index_data <- test_env$survey_index
  test_env$survey_fleet_age_comp <- new(test_env$fims$AgeComp, length(test_env$survey_index), om_input$nages)
  test_env$survey_fleet_age_comp$age_comp_data <- c(t(em_input$survey.age.obs$survey1)) * em_input$n.survey$survey1

  # Growth
  test_env$ewaa_growth <- new(test_env$fims$EWAAgrowth)
  test_env$ewaa_growth$ages <- om_input$ages
  test_env$ewaa_growth$weights <- om_input$W.mt

  # Maturity
  test_env$maturity <- new(test_env$fims$LogisticMaturity)
  test_env$maturity$inflection_point$value <- om_input$A50.mat
  test_env$maturity$inflection_point$is_random_effect <- FALSE
  test_env$maturity$inflection_point$estimated <- FALSE
  test_env$maturity$slope$value <- om_input$slope
  test_env$maturity$slope$is_random_effect <- FALSE
  test_env$maturity$slope$estimated <- FALSE

  # Fleet
  # Create the fishing fleet
  test_env$fishing_fleet_selectivity <- new(test_env$fims$LogisticSelectivity)
  test_env$fishing_fleet_selectivity$inflection_point$value <- om_input$sel_fleet$fleet1$A50.sel1
  test_env$fishing_fleet_selectivity$inflection_point$is_random_effect <- FALSE
  # turn on estimation of inflection_point
  test_env$fishing_fleet_selectivity$inflection_point$estimated <- TRUE
  test_env$fishing_fleet_selectivity$slope$value <- om_input$sel_fleet$fleet1$slope.sel1
  # turn on estimation of slope
  test_env$fishing_fleet_selectivity$slope$is_random_effect <- FALSE
  test_env$fishing_fleet_selectivity$slope$estimated <- TRUE

  test_env$fishing_fleet <- new(test_env$fims$Fleet)
  test_env$fishing_fleet$nages <- om_input$nages
  test_env$fishing_fleet$nyears <- om_input$nyr
  test_env$fishing_fleet$log_Fmort <- log(om_output$f)
  test_env$fishing_fleet$estimate_F <- TRUE
  test_env$fishing_fleet$random_F <- FALSE
  test_env$fishing_fleet$log_q <- log(1.0)
  test_env$fishing_fleet$estimate_q <- FALSE
  test_env$fishing_fleet$random_q <- FALSE
  test_env$fishing_fleet$log_obs_error <- rep(log(sqrt(log(em_input$cv.L$fleet1^2 + 1))), om_input$nyr)
  test_env$fishing_fleet$estimate_obs_error <- FALSE
  # Modules are linked together using module IDs
  # Each module has a get_id() function that returns the unique ID for that module
  # Each fleet uses the module IDs to link up the correct module to the correct fleet
  # Note: Likelihoods not yet set up as a stand-alone modules, so no get_id()
  test_env$fishing_fleet$SetAgeCompLikelihood(1)
  test_env$fishing_fleet$SetIndexLikelihood(1)
  test_env$fishing_fleet$SetSelectivity(test_env$fishing_fleet_selectivity$get_id())
  test_env$fishing_fleet$SetObservedIndexData(test_env$fishing_fleet_index$get_id())
  test_env$fishing_fleet$SetObservedAgeCompData(test_env$fishing_fleet_age_comp$get_id())

  # Create the survey fleet
  test_env$survey_fleet_selectivity <- new(test_env$fims$LogisticSelectivity)
  test_env$survey_fleet_selectivity$inflection_point$value <- om_input$sel_survey$survey1$A50.sel1
  test_env$survey_fleet_selectivity$inflection_point$is_random_effect <- FALSE
  # turn on estimation of inflection_point
  test_env$survey_fleet_selectivity$inflection_point$estimated <- TRUE
  test_env$survey_fleet_selectivity$slope$value <- om_input$sel_survey$survey1$slope.sel1
  test_env$survey_fleet_selectivity$slope$is_random_effect <- FALSE
  # turn on estimation of slope
  test_env$survey_fleet_selectivity$slope$estimated <- TRUE

  test_env$survey_fleet <- new(test_env$fims$Fleet)
  test_env$survey_fleet$is_survey <- TRUE
  test_env$survey_fleet$nages <- om_input$nages
  test_env$survey_fleet$nyears <- om_input$nyr
  test_env$survey_fleet$estimate_F <- FALSE
  test_env$survey_fleet$random_F <- FALSE
  test_env$survey_fleet$log_q <- log(om_output$survey_q$survey1)
  test_env$survey_fleet$estimate_q <- TRUE
  test_env$survey_fleet$random_q <- FALSE
  test_env$survey_fleet$log_obs_error <- rep(log(sqrt(log(em_input$cv.survey$survey1^2 + 1))), om_input$nyr)
  test_env$survey_fleet$estimate_obs_error <- FALSE
  test_env$survey_fleet$SetAgeCompLikelihood(1)
  test_env$survey_fleet$SetIndexLikelihood(1)
  test_env$survey_fleet$SetSelectivity(test_env$survey_fleet_selectivity$get_id())
  test_env$survey_fleet$SetObservedIndexData(test_env$survey_fleet_index$get_id())
  test_env$survey_fleet$SetObservedAgeCompData(test_env$survey_fleet_age_comp$get_id())

  # Population
  test_env$population <- new(test_env$fims$Population)
  test_env$population$log_M <- rep(log(om_input$M.age[1]), om_input$nyr * om_input$nages)
  test_env$population$estimate_M <- FALSE
  test_env$population$log_init_naa <- log(om_output$N.age[1, ])
  test_env$population$estimate_init_naa <- TRUE
  test_env$population$nages <- om_input$nages
  test_env$population$ages <- om_input$ages
  test_env$population$nfleets <- sum(om_input$fleet_num, om_input$survey_num)
  test_env$population$nseasons <- 1
  test_env$population$nyears <- om_input$nyr
  test_env$population$SetMaturity(test_env$maturity$get_id())
  test_env$population$SetGrowth(test_env$ewaa_growth$get_id())
  test_env$population$SetRecruitment(test_env$recruitment$get_id())

  # end of setup_fims function, returning test_env
  return(test_env)
}

# integration test results validation function
validate_integration_test_results <- function(tmb_obj){
  load(test_path("fixtures", "integration_test_data.RData"))
  obj <- tmb_obj
  # Call report using MLE parameter values
  # obj$report() requires parameter list to avoid errors
  report <- obj$report(obj$env$last.par.best)
  sdr <- TMB::sdreport(obj)
  sdr_report <- summary(sdr, "report")
  # Numbers at age
  # Estimates and SE for NAA
  sdr_naa <- sdr_report[which(rownames(sdr_report) == "NAA"), ]
  naa_are <- rep(0, length(c(t(om_output$N.age))))
  for (i in 1:length(c(t(om_output$N.age)))) {
    naa_are[i] <- abs(sdr_naa[i, 1] - c(t(om_output$N.age))[i])
  }
  # Expect 95% of absolute error to be within 2*SE of NAA
  expect_lte(
    sum(naa_are > qnorm(.975) * sdr_naa[1:length(c(t(om_output$N.age))), 2]),
    0.05 * length(c(t(om_output$N.age)))
  )

  # Biomass
  sdr_biomass <- sdr_report[which(rownames(sdr_report) == "Biomass"), ]
  biomass_are <- rep(0, length(om_output$biomass.mt))
  for (i in 1:length(om_output$biomass.mt)) {
    biomass_are[i] <- abs(sdr_biomass[i, 1] - om_output$biomass.mt[i]) # / om_output$biomass.mt[i]
    # expect_lte(biomass_are[i], 0.15)
  }
  expect_lte(
    sum(biomass_are > qnorm(.975) * sdr_biomass[1:length(om_output$biomass.mt), 2]),
    0.05 * length(om_output$biomass.mt)
  )

  # Spawning biomass
  sdr_sb <- sdr_report[which(rownames(sdr_report) == "SSB"), ]
  sb_are <- rep(0, length(om_output$SSB))
  for (i in 1:length(om_output$SSB)) {
    sb_are[i] <- abs(sdr_sb[i, 1] - om_output$SSB[i]) # / om_output$SSB[i]
    # expect_lte(sb_are[i], 0.15)
  }
  expect_lte(
    sum(sb_are > qnorm(.975) * sdr_sb[1:length(om_output$SSB), 2]),
    0.05 * length(om_output$SSB)
  )

  # Recruitment
  fims_naa <- matrix(report$naa[[1]][1:(om_input$nyr * om_input$nages)],
                     nrow = om_input$nyr, byrow = TRUE
  )
  sdr_naa1_vec <- sdr_report[which(rownames(sdr_report) == "NAA"), 2]
  sdr_naa1 <- sdr_naa1_vec[seq(1, om_input$nyr * om_input$nages, by = om_input$nages)]
  fims_naa1_are <- rep(0, om_input$nyr)
  for (i in 1:om_input$nyr) {
    fims_naa1_are[i] <- abs(fims_naa[i, 1] - om_output$N.age[i, 1]) # /
    # om_output$N.age[i, 1]
    # expect_lte(fims_naa1_are[i], 0.25)
  }
  expect_lte(
    sum(fims_naa1_are > qnorm(.975) * sdr_naa1[1:length(om_output$SSB)]),
    0.05 * length(om_output$SSB)
  )

  expect_equal(
    fims_naa[, 1],
    report$recruitment[[1]][1:om_input$nyr]
  )

  # recruitment log deviations
  # the initial value of om_input$logR.resid is dropped from the model
  sdr_rdev <- sdr_report[which(rownames(sdr_report) == "LogRecDev"), ]
  rdev_are <- rep(0, length(om_input$logR.resid) - 1)

  for (i in 1:(length(report$log_recruit_dev[[1]]) - 1)) {
    rdev_are[i] <- abs(report$log_recruit_dev[[1]][i] - om_input$logR.resid[i + 1]) # /
    #   exp(om_input$logR.resid[i])
    # expect_lte(rdev_are[i], 1) # 1
  }
  expect_lte(
    sum(rdev_are > qnorm(.975) * sdr_rdev[1:length(om_input$logR.resid) - 1, 2]),
    0.05 * length(om_input$logR.resid)
  )

  # F (needs to be updated when std.error is available)
  sdr_F <- sdr_report[which(rownames(sdr_report) == "FMort"), ]
  f_are <- rep(0, length(om_output$f))
  for (i in 1:length(om_output$f)) {
    f_are[i] <- abs(sdr_F[i, 1] - om_output$f[i])
  }
  # Expect 95% of absolute error to be within 2*SE of Fmort
  expect_lte(
    sum(f_are > qnorm(.975) * sdr_F[1:length(om_output$f), 2]),
    0.05 * length(om_output$f)
  )

  # Expected fishery catch and survey index
  fims_index <- sdr_report[which(rownames(sdr_report) == "ExpectedIndex"), ]
  fims_catch <- fims_index[1:om_input$nyr, ]
  fims_survey <- fims_index[(om_input$nyr + 1):(om_input$nyr * 2), ]

  # Expected fishery catch - om_output
  catch_are <- rep(0, length(om_output$L.mt$fleet1))
  for (i in 1:length(om_output$L.mt$fleet1)) {
    catch_are[i] <- abs(fims_catch[i, 1] - om_output$L.mt$fleet1[i])
  }
  # Expect 95% of absolute error to be within 2*SE of fishery catch
  expect_lte(
    sum(catch_are > qnorm(.975) * fims_catch[, 2]),
    0.05 * length(om_output$L.mt$fleet1)
  )

  # Expected fishery catch - em_input
  catch_are <- rep(0, length(em_input$L.obs$fleet1))
  for (i in 1:length(em_input$L.obs$fleet1)) {
    catch_are[i] <- abs(fims_catch[i, 1] - em_input$L.obs$fleet1[i])
  }
  # Expect 95% of absolute error to be within 2*SE of fishery catch
  expect_lte(
    sum(catch_are > qnorm(.975) * fims_catch[, 2]),
    0.05 * length(em_input$L.obs$fleet1)
  )


  # Expected fishery catch number at age
  sdr_cnaa <- sdr_report[which(rownames(sdr_report) == "CNAA"), ]
  cnaa_are <- rep(0, length(c(t(om_output$L.age$fleet1))))
  for (i in 1:length(c(t(om_output$L.age$fleet1)))) {
    cnaa_are[i] <- abs(sdr_cnaa[i, 1] - c(t(om_output$L.age$fleet1))[i])
  }
  # Expect 95% of absolute error to be within 2*SE of CNAA
  expect_lte(
    sum(cnaa_are > qnorm(.975) * sdr_cnaa[, 2]),
    0.05 * length(c(t(om_output$L.age$fleet1)))
  )


  # Expected catch number at age in proportion
  # fims_cnaa <- matrix(report$cnaa[1:(om_input$nyr*om_input$nages), 1],
  #                     nrow = om_input$nyr, byrow = TRUE)
  # fims_cnaa_proportion <- fims_cnaa/rowSums(fims_cnaa)
  # for (i in 1:length(c(t(em_input$L.age.obs$fleet1)))){
  #
  #   if (c(t(em_input$L.age.obs$fleet1))[i] == 0) {
  #     expect_lte(abs(c(t(fims_cnaa_proportion))[i] - (c(t(em_input$L.age.obs$fleet1))[i]+0.001))/
  #                  (c(t(em_input$L.age.obs$fleet1))[i]+0.001), 18)
  #   } else {
  #     expect_lte(abs(c(t(fims_cnaa_proportion))[i] - c(t(em_input$L.age.obs$fleet1))[i])/
  #                  c(t(em_input$L.age.obs$fleet1))[i], 3) # Inf when i = 299; landings was 0.
  #   }
  # }


  # Expected survey index - om_output
  index_are <- rep(0, length(om_output$survey_index_biomass$survey1))
  for (i in 1:length(om_output$survey_index_biomass$survey1)) {
    index_are[i] <- abs(fims_survey[i, 1] - om_output$survey_index_biomass$survey1[i])
  }
  # Expect 95% of absolute error to be within 2*SE of survey index
  expect_lte(
    sum(index_are > qnorm(.975) * fims_survey[, 2]),
    0.05 * length(om_output$survey_index_biomass$survey1)
  )

  # Expected survey index - em_input
  index_are <- rep(0, length(em_input$surveyB.obs$survey1))
  for (i in 1:length(em_input$surveyB.obs$survey1)) {
    index_are[i] <- abs(fims_survey[i, 1] - em_input$surveyB.obs$survey1[i])
  }
  # Expect 95% of absolute error to be within 2*SE of survey index
  # expect_lte(
  #   sum(index_are > qnorm(.975) * fims_survey[, 2]),
  #   0.05 * length(em_input$surveyB.obs$survey1)
  # )

  for (i in 1:length(em_input$surveyB.obs$survey1)) {
    expect_lte(abs(fims_survey[i, 1] - em_input$surveyB.obs$survey1[i]) /
                 em_input$surveyB.obs$survey1[i], 0.25)
  }

  # Expected survey number at age
  # for (i in 1:length(c(t(om_output$survey_age_comp$survey1)))){
  #   expect_lte(abs(report$cnaa[i,2] - c(t(om_output$survey_age_comp$survey1))[i])/
  #                c(t(om_output$survey_age_comp$survey1))[i], 0.001)
  # }

  # Expected catch number at age in proportion
  # fims_cnaa <- matrix(report$cnaa[1:(om_input$nyr*om_input$nages), 2],
  #                     nrow = om_input$nyr, byrow = TRUE)
  # fims_cnaa_proportion <- fims_cnaa/rowSums(fims_cnaa)
  #
  # for (i in 1:length(c(t(em_input$survey.age.obs)))){
  #   expect_lte(abs(c(t(fims_cnaa_proportion))[i] - c(t(em_input$L.age.obs$fleet1))[i])/
  #                c(t(em_input$L.age.obs$fleet1))[i], 0.15)
  # }
}
