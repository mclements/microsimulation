\name{enum}
\alias{enum<-}
\alias{enum}
\title{
  Use 0-based \code{enum} factors.
}
\description{
  Use 0-based \code{enum} factors, particularly for use with \code{C++} code.
  
}
\usage{
enum(obj, labels)
enum<-(obj, value)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
    an integer vector representing the levels of the 0-based factor
}
  \item{labels}{
    a vector representing the labels of the factor
}
  \item{value}{
    as per \code{labels} (\code{value} name required for the assignment)
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
y <- enum(rep(0:3,each=2), labels = "")
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
