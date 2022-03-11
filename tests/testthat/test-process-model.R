context('process-model')

test_that('components have correct dimensions', {
  model <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations,
    sensitivities
  )
  for (name in c(
    'control_emissions',
    'perturbations',
    'control_mole_fraction',
    'a_prior',
    'w_prior'
  )) {
    expect_false(is.null(model[[name]]))
  }
  for (name in c('a', 'w', 'eta')) {
    expect_true(is.null(model[[name]]))
  }
  n_alpha <- length(regions) * length(month_starts)
  expect_equal(dim(model$H), c(nrow(control_mole_fraction), n_alpha))
  expect_equal(nrow(model$Psi), nrow(control_mole_fraction))
  expect_equal(dim(model$eta_prior_precision), rep(ncol(model$Psi), 2))
})

test_that('samples have correct dimensions', {
  model <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations,
    sensitivities
  )
  process_sample <- generate(model)
  expect_length(process_sample$a, 2)
  expect_length(process_sample$alpha, ncol(model$H))
  expect_length(process_sample$eta, ncol(model$Psi))
})

test_that('update works', {
  model <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations,
    sensitivities
  )
  updated_model <- update(
    model,
    lag = months(1),
    eta = 1,
    alpha = 1,
    sensitivities = sensitivities
  )

  expect_false(max(abs(updated_model$H - model$H)) == 0)
  expect_equal(updated_model$alpha, rep(1, ncol(model$H)))
  expect_equal(updated_model$eta, rep(1, ncol(model$Psi)))
})

test_that('Error if control mole fraction is not ordered by model_id', {
  f <- function() {
    flux_process_model(
    control_emissions,
    control_mole_fraction %>% arrange(desc(time)),
    perturbations,
    sensitivities
  )}

  expect_error(f(),
    'control_mole_fraction must be ordered by model_id',
    fixed = TRUE)
})

test_that('H is always ordered the same regardless of order in control ems', {
  process_model <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations,
    sensitivities
  )

  process_model_arranged <- flux_process_model(
    control_emissions %>% arrange(desc(month_start)),
    control_mole_fraction,
    perturbations,
    sensitivities
  )

  expect_equal(process_model$H, process_model_arranged$H)
})

test_that('H is always ordered the same regardless of order in perturbations', {
  process_model <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations,
    sensitivities
  )

  process_model_arranged <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations %>% arrange(desc(from_month_start)),
    sensitivities
  )

  expect_equal(process_model$H, process_model_arranged$H)
})

test_that('H is always ordered the same regardless of order in sensitivities', {
  process_model <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations,
    sensitivities
  )

  process_model_arranged <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations,
    sensitivities %>% arrange(desc(from_month_start))
  )

  expect_equal(process_model$H, process_model_arranged$H)
})

test_that('H is ordered by the model id', {
  process_model <- flux_process_model(
    control_emissions,
    control_mole_fraction,
    perturbations,
    sensitivities
  )

  # row index should be model_id
  row_id <- sensitivities$model_id[[which(sensitivities$co2_sensitivity == 2)]]
  # work out col index
  nregion <- length(unique(sensitivities$region))
  mask <- which(sensitivities$co2_sensitivity == 2)
  month_2 <- month(sensitivities$from_month_start[[mask]]) -
             month(min(sensitivities$from_month_start))
  region_2 <- sensitivities$region[[which(sensitivities$co2_sensitivity == 2)]]
  col_id <- month_2 * nregion + region_2

  expect_equal(c(row_id, col_id),
               as.numeric(which(process_model$H == 2, arr.ind = TRUE)))
})

test_that('indexing is correct in output', {
  # 2 in sensitivities means third value should be +1
  # second value in output should be +1
  control_mole_fraction_2 <- data.frame(control_mole_fraction)
  control_mole_fraction_2$co2[[2]] <- 1

  # build process model
  model <- flux_process_model(
    control_emissions,
    control_mole_fraction_2,
    perturbations,
    sensitivities
  )

  process_sample <- generate(model)
  # as all sensitivities (bar the one 2) are 1,
  # expect not altered output to be nregions*nmonths*1
  process_sample$alpha <- rep(1, dim(model$H)[[2]])
  # eta makes things confusing, turn it off
  process_sample$eta <- rep(0, length(process_sample$eta))

  # work out expected_value
  expected_value <- rep(dim(model$H)[[2]], dim(control_mole_fraction)[[1]])
  expected_value[[2]] <- dim(model$H)[[2]] + 1
  expected_value[[3]] <- dim(model$H)[[2]] + 1

  expect_equal(calculate(model, 'Y2', process_sample), expected_value)

  })
