library(microsimulation)
context("callFhcrc")

## To run all the test either:
## R CMD check
## or (faster)
## test_dir("/home/andkar/src/ki/microsimulation/tests/.")

test_returned_object_structure <- function(obj = obj) {
    test_that(paste("Check the structure of the returned fhcrc object with scenario:", try(obj$screen)), {
        expect_is(obj, "fhcrc")
        expect_output(str(obj), "List of 13")
        expect_output(str(obj), "$ n",                     fixed = TRUE)
        expect_output(str(obj), "$ screen",                fixed = TRUE)
        expect_output(str(obj), "$ enum",                  fixed = TRUE)
        expect_output(str(obj), "$ lifeHistories",         fixed = TRUE)
        expect_output(str(obj), "$ parameters",            fixed = TRUE)
        expect_output(str(obj), "$ summary",               fixed = TRUE)
        expect_output(str(obj), "$ healthsector.costs",    fixed = TRUE)
        expect_output(str(obj), "$ societal.costs",        fixed = TRUE)
        expect_output(str(obj), "$ psarecord",             fixed = TRUE)
        expect_output(str(obj), "$ diagnoses",             fixed = TRUE)
        expect_output(str(obj), "$ cohort",                fixed = TRUE)
        expect_output(str(obj), "$ simulation.parameters", fixed = TRUE)
        expect_output(str(obj), "$ falsePositives",        fixed = TRUE)
    })
}

test_scenario <- function(screen){
    test_that(paste("Check microsimulation scenario:", screen), {

        ## N.b. double negation expect_failure() expect_error() <=> expect no error
        expect_failure(expect_error(sim <- callFhcrc(screen = screen)))
        ## Nested check of return object for each scenario
        test_returned_object_structure(sim)

    })
}

scenarios <-  c("noScreening", "randomScreen50to70", "twoYearlyScreen50to70", "fourYearlyScreen50to70", "screen50", "screen60", "screen70", "screenUptake", "stockholm3_goteborg", "stockholm3_risk_stratified", "goteborg", "risk_stratified", "mixed_screening", "regular_screen", "single_screen")

lapply(scenarios, function(x) test_scenario(screen = x))


test_that("Check that input on the R side variables are equal on the C++ side", {
    ## TBA
})

## WIP
## expect_output(print(sim),
##               "FHCRC prostate cancer model with 10 individual(s) under scenario 'noScreening'.")
