library(testthat)
library(microsimulation)

test_check("microsimulation")

## Information on how to write test:
## http://r-pkgs.had.co.nz/tests.html


## What to test

##     Whenever you are tempted to type something into a print
##     statement or a debugger expression, write it as a test
##     instead. — Martin Fowler

## There is a fine balance to writing tests. Each test that you write
## makes your code less likely to change inadvertently; but it also
## can make it harder to change your code on purpose. It’s hard to
## give good general advice about writing tests, but you might find
## these points helpful:

##   * Focus on testing the external interface to your functions - if
##     you test the internal interface, then it’s harder to change the
##     implementation in the future because as well as modifying the
##     code, you’ll also need to update all the tests.

##   * Strive to test each behaviour in one and only one test. Then if
##     that behaviour later changes you only need to update a single
##     test.

##   * Avoid testing simple code that you’re confident will
##     work. Instead focus your time on code that you’re not sure
##     about, is fragile, or has complicated interdependencies. That
##     said, I often find I make the most mistakes when I falsely
##     assume that the problem is simple and doesn’t need any tests.

##   * Always write a test when you discover a bug. You may find it
##     helpful to adopt the test-first philosphy. There you always
##     start by writing the tests, and then write the code that makes
##     them pass. This reflects an important problem solving strategy:
##     start by establishing your success critieria, how you know if
##     you’ve solved the problem.
