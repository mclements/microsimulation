library(microsimulation)
context("callIllnessDeath")

## To run all the test either:
## R CMD check
## or (faster)
## test_dir("/home/andkar/src/ki/microsimulation/tests/.")

test_returned_object_structure <- function(obj = obj) {
    test_that(paste("Check the structure of the Illness-Death object:", try(obj$screen)), {
        expect_is(obj, "list")
        expect_output(str(obj), "List of 4")
        expect_output(str(obj), "$ pt",                      fixed = TRUE)
        expect_output(str(obj), "$ ut",                      fixed = TRUE)
        expect_output(str(obj), "$ events",                  fixed = TRUE)
        expect_output(str(obj), "$ prev",                    fixed = TRUE)
    })
}

test_speed <- function(time_str){
    test_that("Check that the execution speed was not doubled", {
        ## As adviced we skip timing check on CRAN:
        ## http://r-pkgs.had.co.nz/tests.html
        ## To run locally: Sys.setenv(NOT_CRAN='true')
        skip_on_cran()

        ## Cut-of value arbitarily set. Note that this could fail on
        ## really slow systems.
        expect_true(time_str[2] <  0.1)
    })
}

test_callIllnessDeath <- function(){
    test_that("Check Illness-Death model:", {

        ## Make sure no errors are returned. N.b. double negation
        ## expect_failure() expect_error() <=> expect no error
        expect_failure(expect_error(time_str <- system.time(
                                        sim <- callIllnessDeath())))

        ## Check return object
        test_returned_object_structure(sim)

        ## Nested check of execution time
        test_speed(time_str)
    })
}

test_callIllnessDeath()
