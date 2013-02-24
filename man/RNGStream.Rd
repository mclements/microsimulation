\name{RNGStream}
\alias{RNGStream}
\alias{with.RNGStream}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function returning a closure with different methods for using random
  number streams.
}
\description{
  
  \code{RNGStream()} function returns a closure that includes methods to
  reimplements most of the RngStreams library functions. The function
  builds on the lightweight random number streams offered by
  \code{parallel::nextRNGStream()} and
  \code{parallel::nextRNGSubStream()}.

}
\usage{
RNGStream(nextStream = TRUE, iseed = NULL)

\method{with}{RNGStream}(stream, expr, ...)
}
\arguments{
  \item{nextStream}{Boolean to say whether to use the next stream
    [default] or the current stream.}
  \item{iseed}{specify the random number seed as an integer vector of
    length six.}
  \item{stream}{RNGStream to use for evaluating an expression involving
    random numbers.}
  \item{expr}{expression to evaluate}
  \item{...}{arguments to be passed to \code{eval()}}
}
% \details{
% }
\value{
\code{RNGStream()} returns a closure with methods:
 \item{\code{resetRNGkind}}{Return the \code{RNGkind} and \code{.Random.seed} to the
   state when \code{RNGStream} was first called.}
 \item{\code{open}}{Open (or select) the stream for use.}
 \item{\code{close}}{Close (or deselect) the stream for use.}
 \item{\code{resetStream}}{Reset the stream to its beginning.}
 \item{\code{resetSubStream}}{Reset the sub-stream to its beginning.}
 \item{\code{nextStream}}{Move to the start of the next stream.}
 \item{\code{nextSubStream}}{Move to the start of the next sub-stream.}
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Mark Clements.
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
## set up one stream
s1 <- RNGStream()
s1$open()
rnorm(1)
s1$nextSubStream()
rnorm(1)
## reset the stream
s1$resetStream()
rnorm(2)
s1$nextSubStream()
rnorm(2)

## now do with two streams
s1$resetStream()
s2 <- RNGStream()
with(s1,rnorm(1))
with(s2,rnorm(1))
s1$nextSubStream()
with(s1,rnorm(1))
## now reset the streams and take two samples each time
s1$resetStream()
s2$resetStream()
with(s1,rnorm(2))
with(s2,rnorm(2))
s1$nextSubStream()
with(s1,rnorm(2))

## tidy up
rm(s1,s2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
