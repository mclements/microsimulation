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

test_speed <- function(time_str){
    test_that("Check that the execution speed was not doubled", {
        ## As adviced we skip timing check on CRAN:
        ## http://r-pkgs.had.co.nz/tests.html
        ## To run locally: Sys.setenv(NOT_CRAN='true')
        skip_on_cran()

        ## Cut-of value arbitarily set to double the execution speed,
        ## 0.011s, measured 2017-01-26. Note that this could fail on
        ## really slow systems.
        expect_true(as.double(unlist(
            strsplit(time_str[2], "[[:space:]]"))[9]) < 2 * 0.011)
    })
}

test_scenario <- function(screen){
    test_that(paste("Check microsimulation scenario:", screen), {

        ## Make sure no errors are returned. N.b. double negation
        ## expect_failure() expect_error() <=> expect no error
        expect_failure(expect_error(time_str <- capture.output(
                                        sim <- callFhcrc(screen = screen))))

        ## Nested check of return object for each scenario
        test_returned_object_structure(sim)

        ## Nested check of execution time
        test_speed(time_str)
    })
}

scenarios <-  c("noScreening", "randomScreen50to70",
               "twoYearlyScreen50to70", "fourYearlyScreen50to70",
               "screen50", "screen60", "screen70", "screenUptake",
               "stockholm3_goteborg", "stockholm3_risk_stratified",
               "goteborg", "risk_stratified", "mixed_screening",
               "regular_screen", "single_screen")

lapply(scenarios, function(x) test_scenario(screen = x))


## test_that("Check that input on the R side variables are equal on the C++ side", {
##     ## TBA
## })
