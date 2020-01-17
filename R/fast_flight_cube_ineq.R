#' Title
#'
#' @param X
#' @param pik
#' @param B
#' @param r
#' @param deepness
#' @param EPS
#'
#' @return
#' @export
#'
#' @useDynLib IneqCube, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' rm(list = ls())
#' EPS=0.0000001
#' N=1000
#' n=300
#' p=1
#' q=7
#' z=runif(N)
#' #z=rep(1,N)
#' pik=sampling::inclusionprobabilities(z,n)
#' X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
#' A=X/pik
#' Z=cbind(matrix(rbinom(N*q,1,1/2),c(N,q)))
#' B=cbind(Z,-Z)
#' r=c(ceiling(pik%*%B))
#' r[abs(pik%*%B-round(pik%*%B))<EPS]=round(pik%*%B)[abs(pik%*%B-round(pik%*%B))<EPS]
#' s=fast.flight.cube.ineq(X,pik,B,r,deepness=1,EPS=0.000001)
#' round(s,3)
#' c(t(B)%*%pik)
#' c(r)
#' c(t(B)%*%s)
#' colSums(X)
#' c(t(X)%*%(s/pik))
#'
#' piks=fast.flight.cube(cbind(X,B*pik),pik,deepness=1,EPS=0.000001)
fast.flight.cube.ineq<-function(X,pik,B,r,deepness=1,EPS=0.000001){
  A=as.matrix(X/pik)
  if(nrow(A)==0){
    A=matrix(0,c(length(pik),1))
  }
  ######
  B=as.matrix(B/rep(1,length(pik)))
  TT =! (nrow(B)==0)
  if(TT){
    c=c(t(B)%*%pik)
  }
  ######
  if(TT){
    TE=abs(c(t(B)%*%pik)-r)<=EPS
    A=cbind(A,B[,TE])
    r=r[!TE]
    B=B[,!TE]
    c=c[!TE]
  }
  B=as.matrix(B/rep(1,length(pik)))
  TT =! (nrow(B)==0)
  if(is.null(A)){
    prof=deepness
  }else{
    if(is.matrix(A)){
      prof=ncol(A)+deepness
    }else{
      prof=1+deepness
    }
  }
  TEST=(EPS<pik) & (pik<1-EPS)
  prof2=min(sum(TEST),prof)
  if(TT & prof2!=0){
    BR = matrix(B[TEST,],c(sum(TEST),length(B[TEST,])/sum(TEST)))[1:prof2,];
    if(ncol(B)!=0){
      BR=matrix(BR,c(length(BR)/ncol(B),ncol(B)));
    }
  }
  if(prof2==0){
    a=0
  }else{
    pikR=pik[TEST][1:prof2];
    AR=matrix(A[TEST,],c(sum(TEST),length(A[TEST,])/sum(TEST)))[1:prof2,];
    AR=matrix(AR,c(length(AR)/ncol(A),ncol(A)));
    if(nrow(AR)==1 & sum(AR)==0){
      NN=matrix(1,c(1,1))
    }else{
      NN=MASS::Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))));
    }
    a=ncol(NN)
  }

  # print(AR);
  # print(BR);

  while(a>0.5)
  {
    #####
    #####
    u=NN[,1]
    l1=min(pmax((1-pikR)/u,-pikR/u))
    l2=min(pmax((pikR-1)/u,pikR/u))
    if(TT ){
      c=c(t(B)%*%pik)
      # cR=c(t(BR)%*%pikR)
      pet=(r-c)/c(t(BR)%*%u)
      # print(pet);
      # cat("\n");
      if(sum(pet>0)>0) nu1=min(pet[pet>0])  else nu1=100000000000
      if(sum(pet<0)>0) nu2=min(-pet[pet<0]) else nu2=100000000000
      l1=min(l1,nu1)
      l2=min(l2,nu2)
    }
    pik[TEST][1:prof2]=if (runif(1)<l2/(l1+l2)) pikR+l1*u else pikR-l2*u
    TEST=(EPS<pik) & (pik<1-EPS)
    #####
    if(TT)
    {
      TE=abs(c(t(B)%*%pik)-r)<=EPS
      A=cbind(A,B[,TE])
      r=r[!TE]
      B=B[,!TE]
      c=c[!TE]
    }
    B=as.matrix(B/rep(1,length(pik)))
    TT=!(nrow(B)==0)
    ###
    prof=ncol(A)+deepness
    prof2=min(sum(TEST),prof)
    if(TT & prof2!=0) {BR=matrix(B[TEST,],c(sum(TEST),length(B[TEST,])/sum(TEST)))[1:prof2,];
    if(ncol(B)!=0) BR=matrix(BR,c(length(BR)/ncol(B),ncol(B)));}
    if(prof2==0) a=0 else {pikR=pik[TEST][1:prof2];
    AR=matrix(A[TEST,],c(sum(TEST),sum(A[TEST,])/sum(TEST)))[1:prof2,];
    AR=matrix(AR,c(length(AR)/ncol(A),ncol(A)));
    if(sum(abs(AR))==0) NN=matrix(1,c(nrow(AR),1)) else NN=MASS::Null(matrix(AR,c(length(AR)/ncol(A),ncol(A))));
    a=ncol(NN)}

    ####
    ######
  }
  pik
}
