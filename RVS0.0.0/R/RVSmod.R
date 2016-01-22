## in this file, I tried to define classes and methods
#' method for RVS_asy
RVSmod <- function (x, ... ) UseMethod ('RVSmod')
RVSmod.default <- function(x,y,p,method='RVS', ...)
{
   x<-as.numeric(x)
   y<-as.numeric(y)
   p<-as.numeric(p)
   p.rvs <- RVS_asy(y,x,p,method)
   p.rvs$call<-match.call()
   class(p.rvs) <-'RVSmod'
   p.rvs
}

print.RVSmod <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n estimates \t pvalue:\n")
  print(c(x$coefficients,x$p_rvs))
}


summary.linmod <- function(object, ...)
{
  se <- sqrt(diag(object$vcov))
  tval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df=object$df))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.linmod"
  res }


print.summary.linmod <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

linmod.formula <- function(formula, data=list(), ...)
{
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  est <- linmod.default(x, y, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}
