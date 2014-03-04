\name{EventQueue}
\docType{methods}
\alias{class-EventQueue}
\alias{pop,EventQueue-method}
\alias{clear,EventQueue-method}
\alias{insert,EventQueue-method}
\alias{remove,EventQueue-method}
\alias{empty,EventQueue-method}
\title{
  R5 class for an event queue
}
\description{
  R5 class for an event queue. In includes methods to:
  \code{insert(time,event)} at a given \code{time} and \code{event};
  \code{pop()} from the queue and return the value; check with the queue
  is \code{empty()}; and \code{remove(predicate)} for events that
  satisfy the \code{predicate}.
  
}
\usage{
\S4method{insert}{EventQueue}(time, event)
\S4method{clear}{EventQueue}()
\S4method{pop}{EventQueue}()
\S4method{remove}{EventQueue}(predicate)
\S4method{empty}{EventQueue}()
}
\arguments{
  \item{time}{Time for an event (numeric).}
  \item{event}{An arbitrary object representing an event in the queue.}
  \item{predicate}{A predicate that takes an event and returns a Boolean.}
}
\section{Methods}{
  \describe{
    \item{\code{insert(time,event)}}{
      Insert and event at a given time into the event queue.
    }
    \item{\code{clear()}}{
      Remove all of the events from the queue.
    }
    \item{\code{pop()}}{
      Remove the first event in the queue and return its value.
    }
    \item{\code{remove(predicate)}}{
      Applies the predicate to all events in the queue to determine
      whether to remove the event from the queue.
    }
    \item{\code{empty()}}{
      Boolean to say whether the queue is empty or not.
    }
  }
}
\details{

  This event queue is simple and useful for pedagogic purposes.

  The algorithm for inserting values into the queue is computationally
  very simple: simply rank the times using \code{order()} and re-order
  times and events. This approach is probably of acceptable performance
  for smaller queue. A more computationally efficient approach for
  insert into larger queues would be to use a binary search (e.g. using
  \code{findInterval()}).

}
% \value{
%   The output class can be
%   }
\examples{
\dontrun{
queue <- EventQueue$new()
queue$insert(2,"B")
queue$insert(1,"A")
queue$insert(3,"C")
while (!queue$empty())
  print(queue$pop())
}
}
\seealso{class-BaseDiscreteEventSimulation}