#' S3 priority queue implementation using C++
#'
#' This provides a priority queue that is sorted by the priority and entry order. The priority is assumed to be numeric. The events can be of any type. As an extension, events can be cancelled if they satisfy a certain predicate. Note that the inactive events are not removed, rather they are marked as cancelled and will not be available to be popped.
#'
#' @param lower boolean to determine whether to give priority to lower values (default=TRUE)
#' or higher values
#' 
#' @return a list with
#' \description{
#' \item{push}{function with arguments priority (numeric) and event (SEXP). Pushes an event with a given priority}
#' \item{pop}{function to return a list with a priority (numeric) and an event (SEXP). This pops the first active event.}
#' \item{cancel}{function that takes a predicate (or R function) for a given event and returns a logical that indicates whether to cancel that event or not. This may cancel some events that will no longer be popped.}
#' \item{empty}{function that returns whether the priority queue is empty (or has no active events).}
#' \item{clear}{function to clear the priority queue.}
#' \item{ptr}{XPtr value}
#' }
#' 
#' @export
#'
#' @examples
#' pq = pqueue()
#' pq$push(3,"Clear drains")
#' pq$push(4, "Feed cat")
#' pq$push(5, "Make tea")
#' pq$push(1, "Solve RC tasks")
#' pq$push(2, "Tax return")
#' while(!pq$empty())
#'   print(pq$pop())
#' 
pqueue <- function(lower = TRUE) {
    stopifnot(is.logical(lower), length(lower)==1)
    ptr <- .Call("pqueue__new", as.logical(lower), PACKAGE="microsimulation")
    push <- function(priority, event) {
        stopifnot(is.numeric(priority), length(priority)==1)
        .Call("pqueue__push", ptr, as.double(priority), event, PACKAGE="microsimulation")
        invisible()
    }
    pop <- function() .Call("pqueue__pop", ptr, PACKAGE="microsimulation")
    cancel <- function(predicate) {
        stopifnot(is.function(predicate))
        .Call("pqueue__cancel", ptr, predicate, PACKAGE="microsimulation")
        invisible()
    }
    empty <- function() {
        .Call("pqueue__empty", ptr, PACKAGE="microsimulation")
    }
    clear <- function() {
        .Call("pqueue__clear", ptr, PACKAGE="microsimulation")
        invisible()
    }
    structure(list(push=push, pop=pop, cancel=cancel, empty=empty, clear=clear, ptr=ptr),
              class="pqueue")
}
