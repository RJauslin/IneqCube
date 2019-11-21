#' Title
#'
#' @param X
#' @param pik
#' @param deepness
#' @param EPS
#'
#' @return
#' @export
#'
#' @examples
fast.flight.cube<-function(X,pik,deepness=1,EPS=0.0000001)
{
  A=as.matrix(X/pik)
  if(nrow(A)==0) A=matrix(0,c(length(pik),1))
  if(is.null(X)) prof=1 else{if(is.matrix(X))prof=ncol(X)+1 else prof=2 }
  TEST=(EPS<pik) & (pik<1-EPS)
  prof2=min(sum(TEST),prof)
  if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
  AR=matrix(A[TEST,],c(sum(TEST),sum(A[TEST,])/sum(TEST)))[1:prof2,];
  NN=Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))));
  a=ncol(NN)}
  while(a>0.5)
  {
    u=NN[,1]
    l1=min(pmax((1-pikR)/u,-pikR/u))
    l2=min(pmax((pikR-1)/u,pikR/u))
    pik[TEST][1:prof2]=if (runif(1)<l2/(l1+l2)) pikR+l1*u else pikR-l2*u
    TEST=(EPS<pik) & (pik<1-EPS)
    prof2=min(sum(TEST),prof)
    if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
    AR=matrix(A[TEST,],c(sum(TEST),sum(A[TEST,])/sum(TEST)))[1:prof2,];
    NN=Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))));
    a=ncol(NN)}
  }
  pik
}
