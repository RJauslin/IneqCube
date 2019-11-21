#' Title
#'
#' @param X
#' @param pik
#' @param EPS
#'
#' @return
#' @export
#'
#' @examples
flight_cube<-function(X,pik,EPS=0.0000001)
{
  A=as.matrix(X/pik)
  if(nrow(A)==0) A=matrix(0,c(length(pik),1))
  TEST=(EPS<pik) & (pik<1-EPS)
  if(sum(TEST)==0) a=0 else {pikR=pik[TEST];AR=A[TEST,];NN=Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))));a=ncol(NN)}
  while(a>0.5)
  {
    u=NN[,1]
    l1=min(pmax((1-pikR)/u,-pikR/u))
    l2=min(pmax((pikR-1)/u,pikR/u))
    pik[TEST] = if (runif(1)<l2/(l1+l2)) pikR+l1*u else pikR-l2*u
    TEST=(EPS<pik) & (pik<1-EPS)
    if(sum(TEST)==0) a=0 else {pikR=pik[TEST];AR=A[TEST,];NN=Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))));a=ncol(NN)}
  }
  pik
}
