library(microsimulation)
context("callFhcrc")

## To run all the test either:
## R CMD check
## or (faster)
## test_dir("/home/andkar/src/ki/microsimulation/tests/.")

test_that("Check the structure of the returned fhcrc object", {
    expect_is(sim <- callFhcrc(), "fhcrc")
    expect_output(str(sim), "List of 13")
    expect_output(str(sim), "$ n",                     fixed = TRUE)
    expect_output(str(sim), "$ screen",                fixed = TRUE)
    expect_output(str(sim), "$ enum",                  fixed = TRUE)
    expect_output(str(sim), "$ lifeHistories",         fixed = TRUE)
    expect_output(str(sim), "$ parameters",            fixed = TRUE)
    expect_output(str(sim), "$ summary",               fixed = TRUE)
    expect_output(str(sim), "$ healthsector.costs",    fixed = TRUE)
    expect_output(str(sim), "$ societal.costs",        fixed = TRUE)
    expect_output(str(sim), "$ psarecord",             fixed = TRUE)
    expect_output(str(sim), "$ diagnoses",             fixed = TRUE)
    expect_output(str(sim), "$ cohort",                fixed = TRUE)
    expect_output(str(sim), "$ simulation.parameters", fixed = TRUE)
    expect_output(str(sim), "$ falsePositives",        fixed = TRUE)
})


test_that("Check that all scenarios run", {

    ## scenarios <-  c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70", "fourYearlyScreen50to70", "screen50", "screen60", "screen70", "screenUptake", "stockholm3_goteborg", "stockholm3_risk_stratified", "goteborg", "risk_stratified", "mixed_screening", "regular_screen", "single_screen")

  ## ## This does not work as expected
  ## ## http://stackoverflow.com/questions/10826365/how-to-test-that-an-error-does-not-occur
  ## expect_error(for (screen in scenarios){expect_error(callFhcrc(), NA, info = paste("screen =", screen))})

  ## expect_error(lapply(c("noScreening", "screen500"), function(x) callFhcrc(screen = x)))

  expect_failure(expect_error(callFhcrc(screen = "noScreening")))
  expect_failure(expect_error(callFhcrc(screen = "randomScreen50to70")))
  expect_failure(expect_error(callFhcrc(screen = "twoYearlyScreen50to70")))
  expect_failure(expect_error(callFhcrc(screen = "fourYearlyScreen50to70")))
  expect_failure(expect_error(callFhcrc(screen = "screen50")))
  expect_failure(expect_error(callFhcrc(screen = "screen60")))
  expect_failure(expect_error(callFhcrc(screen = "screen70")))
  expect_failure(expect_error(callFhcrc(screen = "screenUptake")))
  expect_failure(expect_error(callFhcrc(screen = "stockholm3_goteborg")))
  expect_failure(expect_error(callFhcrc(screen = "stockholm3_risk_stratified")))
  expect_failure(expect_error(callFhcrc(screen = "goteborg")))
  expect_failure(expect_error(callFhcrc(screen = "risk_stratified")))
  expect_failure(expect_error(callFhcrc(screen = "mixed_screening")))
  expect_failure(expect_error(callFhcrc(screen = "regular_screen")))
  expect_failure(expect_error(callFhcrc(screen = "single_screen")))

})

test_that("Check that input on the R side variables are equal on the C++ side", {
# TBA
})

## expect_output(print(sim),
##               "FHCRC prostate cancer model with 10 individual(s) under scenario 'noScreening'.")
