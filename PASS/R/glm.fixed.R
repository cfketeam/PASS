#' logistic variable select
#'
#' @param seed The random seed to keep the repeatability
#' @param x Input matrix,of dimension n × p; each row is an observation vector.
#' @param y Response variable,of length n, should be a factor with two levels.
#' @param sam.time number of repeat times
#' @param n model size
#' @param weight sample percent
#' @param pweight prediction weight in SFS
#' @param sim.percent A parameter to define the mean of the the nearest result using how many percent in the repeating time.
#'
#' @return
#' @returns mean;tmean  The predict performance is shown through the AUC.
#' @returns subset      The best subset column number.
#' @returns allsubset   The subset result of each sample.
#' @returns predictions The prediction performance of each subset.
#' @export
#'
#' @examples
#' K <-10
#' sigma <- 1
#' rho <- 0.5
#' a <- simu.data(n=200,p=1000,family = "binomial", K, rho, sigma,alpha=4)
#' fit <- mbess.glm(seed=101, x=a$x, y=a$y, sam.time=10, n=6 )
#' fit
#'
glm.fixed <- function(x, y, n, sam.time=10, weight=0.8, pweight=10, method="mean", sim.percent=0.5, seed=101)
{
  re <- list()
  all <- c()
  b <- c()
  var_rate <- c()
  var_all <- c()
  set.seed(seed)
  rate <- c()
  data <- data.frame(y,x)
  simr <- rep(0,sam.time)
  value_re <- rep(0,sam.time)
  data <- na.omit(data)

  b_one <- function(y,x,n,i,weight)
  {
    if(is.na(n))
      return(warning("Lack model size"))
    i <- i-100
    a <- 0
    while(a==0)
    {
      i <- i+100
      df <- Bootstrap(data.frame(y,x),i,weight = weight)
      x_t <- df[[1]][,-1]
      y_t <- df[[1]][,1]
      lb <- tryCatch(bess.one(x_t,y_t, family = "binomial", s = n), error = function(e) e)

      if(!is.null(lb$beta))
      {
        lb <- lb$beta
        var <- c(1:ncol(x_t))[lb !=0]
        fullformula <- as.formula(paste0(colnames(data.frame(y,x))[1], "~", paste(colnames(data.frame(y,x))[var+1], collapse = " + ")))
        fit_glm <- glm(fullformula, family = binomial, data = df[[1]])
        prob <- predict.glm(fit_glm, df[[2]], family = binomial, type = "response")
        auc <- tryCatch(roc(df[[2]][,1], prob,quiet=T)$auc, error = function(e) e)
        if (is.numeric(auc) && !is.na(auc))
        {
          a <- 1
        }
      }
    }
    return(re <- list(var=var,auc=auc))
  }


  Bootstrap <- function(df, seed, weight=weight) {
    set.seed(seed)
    re <- list()
    c <- 1:nrow(df)
    a <- sample(c, round(nrow(df)*weight), replace = F)
    t <- c[!(c %in% unique(a))]
    train <- df[unique(a), ]
    re[[1]] <- train
    test <- df[t, ]
    re[[2]] <- test
    return(re)
  }

  cl<- makeCluster(4)
  registerDoParallel(cl)
  varfit <- foreach(i=1:sam.time+seed,.packages = c("BeSS", "pROC") ) %dopar% b_one(y,x,n=n,i=i,weight = weight)
  stopCluster(cl)

  for(i in 1:sam.time)
  {
    all <- c(all,as.numeric(varfit[[i]]$auc))
    var_all <- c(var_all,varfit[[i]]$var)
    var_rate <- cbind(var_rate,varfit[[i]]$var)
  }


  for(i in 1:max(table(var_all)))
  {
    value <- c()
    rate <- names(table(var_all))[table(var_all)==i]
    for(m in 1:sam.time)
    {
      value_t <- sum(var_rate[,m] %in% rate)*(i - 1)
      value <- c(value,value_t)
    }
    value_re <- value_re + value
  }
  value_re <- value_re+pweight*all
  var_re <- var_rate[,c(1:sam.time)[value_re == max(value_re)][1]]

  for (i in unique(table(var_all[var_all %in% var_re])))
  {
    value <- c()
    rate <- names(table(var_all[var_all %in% var_re]))[table(var_all[var_all %in% var_re]) == i]
    for (m in 1:sam.time) {
      value_t <- sum(var_rate[, m] %in% rate) * (i - 1)
      value <- c(value, value_t)
    }
    simr <- simr + value
  }
  rate <- all[order(simr, decreasing = T)][1:(sam.time * sim.percent)]

  mean <- mean(rate)
  tmean <- (sum(all)-min(all)-max(all))/(sam.time-2)
  if (method=="mean")
  {
    performance=mean
    other=tmean
  }
  if (method=="Tmean")
  {
    performance=tmean
    other=mean
  }
  subset <- var_re
  var <- subset
  fullformula <- as.formula(paste0(colnames(data.frame(y,x))[1], "~", paste(colnames(data.frame(y,x))[var+1], collapse = " + ")))
  bestmodel <- glm(fullformula, family = binomial, data = data.frame(y,x))

  re <- list(n=length(subset),subset=subset,bestmodel=bestmodel,performance=performance,other=other)
  return(re)
}

