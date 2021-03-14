pqueue <- function() {
    ptr <- .Call("pqueue__new", PACKAGE="microsimulation")
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
