\name{RNGstate}
\alias{RNGstate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function returning a closure for a random state.
}
\description{
  
  This function returns a simple closure that includes the old seed and
  a method to reset the \code{RNGkind} and \code{.Random.seed} to the
  values that were current when the function was called.

}
\usage{
RNGstate()
}
%- maybe also 'usage' for other objects documented here.
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  Returns a closure with:
  \item{oldseed}{a vector with the previous random number seed}
  \item{reset(}}{a method to reset the \code{RNGkind} and \code{.Random.seed} to the
  values that were current when the function was called.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function () 
% {
%     oldseed <- if (exists(".Random.seed", envir = .GlobalEnv, 
%         inherits = FALSE)) 
%         get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
%     else NULL
%     reset <- function() {
%         if (!is.null(oldseed)) 
%             assign(".Random.seed", oldseed, envir = .GlobalEnv)
%         else rm(.Random.seed, envir = .GlobalEnv)
%     }
%     list(oldseed = oldseed, reset = reset)
%   }
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
