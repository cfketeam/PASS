\name{pass.cox}
\alias{pass.cox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feature subset selection with an automatically defined model size.
%%  ~~function to do ... ~~
}
\description{
Feature subset selection with an automatically defined model size for Cox proportional hazards regression model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pass.cox(x,y ,bar=0.01,method="mean",n.max= NULL,sam.time=10,alpha=3,weight=0.8,sim.percent=0.5,seed=101)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Input matrix of explanatory variables, where rows are the observations and columns are the features.}
  \item{y}{Response variable, , which should be a two-column matrix with columns named 'time' and 'status'.}
  \item{method}{Prediction performance evaluation method: “mean” or “Tmean” (for further descriptions, please refer to details). The default is “mean”.}
  \item{bar}{A parameter to define the significant improvement in prediction performance that makes sense. The default is 0.01.}
  \item{n.max}{The maximum model size allowed. Default is \eqn{min(p, n/log(n))}.}
  \item{sam.time}{Number of sampling times. The default is 10.}
  \item{weight}{The percentage of observations to be sampled. The default is 0.8.}
  \item{sim.percent}{A parameter to define the percentage of included feature subsets  using the “mean” method. These feature subsets, most similar to the optimal feature subset, are used to calculate the mean of their prediction performance as the prediction performance of the optimal feature subset. The default is 0.5.}
  \item{alpha}{A parameter to define the highest prediction value. This parameter indicates the number of subsequent prediction values smaller than significant improvement of the current prediction value. The default is 3.}
  \item{seed}{The random seed to keep the repeatability. The default is 101.}

}
\details{
To evaluate each subset objectively, we define the subset frequency score (SFS) as the following:
\deqn{SFS_m=\sum_{i=1}^{sam.time}{\sum_{j=1}^{n}{I^{m}_{x_i \, in \, subset_n} + pweight * weight}}}
where the weight is the predicted performance of each model in the test set, pweight is a parameter controlling the rate in SFS. The best subset has the highest SFS. \cr The prediction performance is estimated in two different metrics: (1) The mean of the prediction ability index, calculated by excluding the maximum and minimum prediction values among all subsets and then averaging the remaining values, is used as the indicator for evaluating the prediction ability (named as Tmean). (2) the average prediction performance of several candidate subsets that are most similar to the highest consistency scoring subset is used as an indicator for evaluating the prediction ability (named as mean). \cr In order to automatically search for the optimal model size, we evaluate all prediction performance at different model sizes and adopt a parameter to define the highest prediction value. This parameter indicates the number of subsequent prediction values smaller than significant improvement of the current prediction value. Additionally, we introduce a parameter bar to determine whether an improvement is meaningful, which encourages the method to favor smaller subset sizes.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{n}{The optimal model size.}
\item{subset}{The optimal feature subset. }
\item{bestmodel }{The fitted model based on the optimal feature subset, the class of which is "coxph".}
\item{performance}{The prediction performance of the optimal feature subset fitted model, which is shown as the C-index.}
\item{performance.list}{The list of prediction performance by C-index values at different model size. }
\item{othersubet}{The optimal feature subset defined by the other criterion. }
\item{otherperformane}{The prediction performance of the optimal feature subset fitted model  defined by the other prediction evaluation method, which is shown as the C-index.}
\item{other.list}{The list of prediction performance by C-index values at different model size,  estimated through the other prediction evaluation method.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chaofu Ke; Luying wu.
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
k <-3
sigma <- 1
rho <- 0.5
a <- simu.data(n=100,p=1000,family = "cox", k, rho,scal=10, sigma,alpha=3)
fit <- pass.cox(x=a$x, y=a$y,sam.time=10,seed=101,alpha=3,bar=0.01,weight=0.9,n.max=5)
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
