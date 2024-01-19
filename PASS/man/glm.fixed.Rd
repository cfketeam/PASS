\name{glm.fixed}
\alias{glm.fixed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feature subset selection with a specified model size
%%  ~~function to do ... ~~
}
\description{
Feature subset selection with a specified model size for logistic regression model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
glm.fixed(x, y, n, sam.time=10, weight=2/3, pweight=10, method="mean", sim.percent=0.5, seed=101)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Input matrix of explanatory variables, where rows are the observations and columns are the features.}
  \item{y}{Response variable, which should be a factor with two levels.}
  \item{n}{Specified model size.}
  \item{sam.time}{Number of sampling times. The default is 10.}
  \item{weight}{The percentage of observations to be sampled. The default is 0.8.}
  \item{p.weight}{Prediction weight in subset frequency score (SFS). The default is 10.}
  \item{method}{Prediction performance evaluation method: “mean” or “Tmean” (for further descriptions, please refer to details). The default is “mean”.}
  \item{sim.percent}{A parameter to define the percentage of included feature subsets using the “mean” method. These feature subsets, most similar to the optimal feature subset, are used to calculate the mean of their prediction performance as the prediction performance of the optimal feature subset. The default is 0.5.}
  \item{seed}{The random seed to keep the repeatability. The default is 101.}
}
\details{To evaluate each subset objectively, we define the subset frequency score (SFS) as the following:
\deqn{SFS_m=\sum_{i=1}^{sam.time}{\sum_{j=1}^{n}{I^{m}_{x_i \, in \, subset_n} + pweight * weight}}}
where the weight is the predicted performance of each model in the test set, pweight is a parameter controlling the rate in SFS. The best subset has the highest SFS. \cr The prediction performance is estimated in two different metrics: (1) The mean of the prediction ability index, calculated by excluding the maximum and minimum prediction values among all subsets and then averaging the remaining values, is used as the indicator for evaluating the prediction ability (named as Tmean). (2) the average prediction performance of several candidate subsets that are most similar to the highest consistency scoring subset is used as an indicator for evaluating the prediction ability (named as mean).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{n }{The specified model size.}
\item{subset }{The optimal feature subset.}
\item{bestmodel }{The fitted model based on the optimal feature subset, the class of which is "glm".}
\item{performance }{The prediction performance of the best fitted model, which is shown as the AUC value.}
\item{other }{The prediction performance of the best fitted model estimated through the other prediction evaluation method, which is shown as the AUC value.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chaofu Ke, Luying Wu
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

k <-10
sigma <- 1
rho <- 0.5
a <- simu.data(n=200,p=1000,family = "binomial", k, rho, sigma,alpha=3)
fit <- glm.fixed(seed=101, x=a$x, y=a$y, sam.time=10, n=6 )
fit

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
