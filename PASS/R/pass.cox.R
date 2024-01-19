#' logistic variable selection
#'
#'
#' @param x Input matrix,of dimension n × p; each row is an observation vector.
#' @param y Response variable,of length n, should be a two-column matrix with columns named 'time' and 'status'.
#' @param sam.time The parameter determine the the number of repeating time.
#' @param seed The random seed to keep the repeatability.
#' @param alpha A parameter to define the highest peak value.
#' @param bar A parameter to recoginze the improvement whether make sense.
#' @param weight The sample percent of the training set.
#' @param n.max The max model size
#' @param sim.percent A parameter to define the mean of the the nearest result using how many percent in the repeating time.
#'
#' @return
#' @returns re_mean   The result searched by the standard mean
#' @returns re_tmean  The result searched by the standard Trim mean
#' @returns auc_mean  The list of mean of each number size
#' @returns auc_tmean The list of Trim mean of each number size
#' @export
#'
#' @examples
#' K <-10
#' sigma <- 1
#' rho <- 0.5
#' a <- simu.data(n=200,p=1000,family = "cox", K, rho,scal=10, sigma,alpha=4)
#' fit <- cox.varselect(x=a$x, y=a$y,sam.time=10,seed=101,alpha=3,bar=0.01,weight=0.9)
#' fit
#'

pass.cox <- function(x,y ,bar=0.01,method="mean",n.max= NULL,sam.time=10,alpha=3,weight=0.8,sim.percent=0.5,seed=101)
{
  re <- list()
  rlist <- c()
  numlist <- c()
  k <- c(1)
  l=0
  varl <- 1:ncol(x)
  n_min <- 1
  if (is.null(n.max))
    n.max=round(nrow(x)/log(nrow(x)))
  n_max <- min(ncol(x), round(nrow(x)/log(nrow(x))),n.max)
  data <- data.frame(y,x)
  auc_mean <- c()
  auc_tmean <- c()
  re_r <- list()


  while(k[1]<n_max)
  {
    kr <- n_min
    rr <- cox.fixed(x=x,y=y,n=n_min,sam.time=sam.time,weight=weight,sim.percent=sim.percent)
    auc_mean <- c(rr$performance)
    auc_tmean <- c(rr$other)

    ms <- paste0("Model size: 1")
    if(method=="mean")
      pf <- paste0("  C-index: ",format(auc_mean[length(auc_mean)], nsmall = 3,digits = 3))
    if(method=="Tmean")
      pf <- paste0("  C-index: ",format(auc_mean[length(auc_mean)], nsmall = 3,digits = 3))
    cat(paste0(ms,pf),"\n")

    re <- list()
    re_mean <- list()
    re_tmean <- list()
    k_mean <- c()
    k_tmean <- c()
    ris_mean <- 0
    ris_tmean <- 0
    nris_mean <- 0
    nris_tmean <- 0
    re_k <- c(1)
    re_mean[[1]] <- rr
    re_tmean[[1]] <- rr
    t_mean <- rr
    t_tmean <- rr
    t_mean$performance <- 0
    t_tmean$other <- 0
    p_mean <- 1
    p_tmean <- 1
    s_tmean <- 0
    s_mean <- 0

    while(min(s_tmean,s_mean)<2 && kr<n_max)
    {
      rl <- rr
      kl <- kr
      kr <- kr+1
      rr <- cox.fixed(x=x,y=y,sam.time=sam.time,n=kr,weight=weight,sim.percent=sim.percent,seed=seed)
      auc_mean <- c(auc_mean, rr$performance)
      auc_tmean <- c(auc_tmean, rr$other)
      if(rr$performance - re_mean[[p_mean]]$performance<= bar || ris_tmean>=alpha*3 || ris_mean>=alpha*3)
      {
        ris_mean <- ris_mean+1

        if(rr$performance < t_mean$performance)
          nris_mean <- nris_mean+1

        if(((ris_mean>=alpha && nris_mean>=alpha)||ris_mean>=alpha*3 || ris_tmean>=alpha*3 ) && (t_mean$performance < (re_mean[[p_mean]]$performance + bar)))
        {
          s_mean <- 2
        }
      }
      if(rr$other - re_tmean[[p_tmean]]$other<= bar || ris_tmean>=alpha*3 || ris_mean>=alpha*3 )
      {
        ris_tmean <- ris_tmean+1
        if(rr$other < t_tmean$other)
          nris_tmean <- nris_tmean+1
        if(((ris_tmean>=alpha && nris_tmean>=alpha)||ris_tmean>=alpha*3 || ris_mean>=alpha*3 ) && (t_tmean$other < (re_tmean[[p_tmean]]$other + bar)))
        {
          s_tmean <- 2
        }
      }
      if((rr$performance-rl$performance) > 0)
      {
        if(ris_mean<alpha && (rr$performance-re_mean[[p_mean]]$performance) > bar)
        {
          re_mean[[p_mean]] <- rr
          ris_mean <- 0
          nris_mean <- 0
          t_mean$performance <- 0
        }
        else if(ris_mean>=alpha && (rr$performance-re_mean[[p_mean]]$performance) > bar)
        {
          p_mean <- p_mean+1
          re_mean[[p_mean]] <- rr
          ris_mean <- 0
          nris_mean <- 0
          t_mean$performance <- 0
        }
        else if(rr$performance-t_mean$performance > 0 && rl$performance < re_mean[[p_mean]]$performance)
        {
          t_mean <- rr
          nris_mean <- 0
        }

      }

      if((rr$other-rl$other) > 0)
      {
        if(ris_tmean < alpha && (rr$other-re_tmean[[p_tmean]]$other) > bar)
        {
          re_tmean[[p_tmean]] <- rr
          ris_tmean <- 0
          nris_tmean <- 0
          t_tmean$other <- 0
        }
        else if(ris_tmean >= alpha && (rr$other-re_tmean[[p_tmean]]$other) > bar)
        {
          p_tmean <- p_tmean+1
          re_tmean[[p_tmean]] <- rr
          ris_tmean <- 0
          nris_tmean <- 0
          t_tmean$other <- 0
        }
        else if(rr$other-t_tmean$other > 0 &&  rl$other < re_tmean[[p_tmean]]$other)
        {
          t_tmean <- rr
          nris_tmean <- 0
        }

      }
      re_k <- c(re_k,kr)
      ms <- paste0("Model size: ",re_k[length(re_k)])
      if(method=="mean")
        pf <- paste0("  C-index: ",format(auc_mean[length(auc_mean)], nsmall = 3,digits = 3))
      if(method=="Tmean")
        pf <- paste0("  C-index: ",format(auc_mean[length(auc_mean)], nsmall = 3,digits = 3))
      cat(paste0(ms,pf),"\n")

    }


    if(method=="mean")
    {
      best <- re_mean[[length(re_mean)]]
      other <- re_tmean[[length(re_tmean)]]
      performance <- re_mean[[length(re_mean)]]$performance
      performance.list <- auc_mean
      otherperformance <- other$other
      other.list <- auc_tmean
    }

    if(method=="Tmean")
    {
      best <- re_tmean[[length(re_tmean)]]
      other <- re_mean[[length(re_mean)]]
      performance <- re_mean[[length(re_mean)]]$other
      performance.list <- auc_tmean
      otherperformance <- other$performance
      other.list <- auc_mean
    }
    var <- best$subset
    fullformula <- as.formula(paste0("Surv(",colnames(data)[1],",",colnames(data)[2],")", "~", paste(colnames(data)[var+2], collapse = " + ")))
    bestmodel <- coxph(fullformula, data = data.frame(y,x))

    re <- list(n=length(var),subset=best$subset, bestmodel=bestmodel,performance=performance,
               performance.list=performance.list,othersubset=other$subset,otherperformance=otherperformance,
               other.list=other.list)
    return(re)
    break


  }


}
