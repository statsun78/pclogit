\name{sel.pclogit}
\alias{sel.pclogit}
\title{selection probabilities of regression coefficients}
\description{
    The selection probability of each regression coefficient is computed 
    based on resamplings. 
    }
\usage{sel.pclogit(x,y,stra=NULL,...,psub=0.5,N.lam=5,K=100)}
\arguments{
    \item{x}{The design matrix (\eqn{n} x \eqn{p}) without an intercept.
        \code{pclogit} standardizes the data by default, but the coefficients
        are always returned on the original scale.}
    \item{y}{The response variable coded as 1 for cases and 0 for the 
        matched controls.}
    \item{stra}{A vector of consecutive integers indicating the stratum of
        each observation. Each stratum must have exact one case and at least
        one control. If not specified, \code{sel.pclogit} fits an ordinary 
        logistic regression.}
    \item{\dots}{Other arguments that can be passed to \code{pclogit}.}
    \item{psub}{The proportion of subsamples used for resamplings, and 
        \code{psub}\eqn{\in[0.5,1)}. The defalut is 0.5.}
    \item{N.lam}{The number of \code{lambda} values used for resamplings,
        and default is 5.}
    \item{K}{The number of resamplings, and default is 100.}    
        }
\details{
        The half of the strata \code{stra} are randomly selected without
        replacement \code{K} times. For each replication, the paired
        (\code{x},\code{y}) in the selected strata are only
        used for \code{pclogit} to find non-zero coefficients along with
        \code{N.lam} \code{lambda} values. The selection probability of each
        coefficient is then computed based on the proportion of non-zeros
        out of \code{K} replciations. In an ordinary logistic model, the half
        of cases and controls are selected each time.      
        }
\value{
    \item{beta}{The selection prbabilities (\eqn{p} x \code{N.lam})}
    \item{maxsel}{The maximum selection probability of each coefficient
      are listed in descending order along with the corresponding variable.}
    \item{lambda}{The actual sequence of \code{lambda} values used}
    \item{K}{The actual number of resamplings used}
     }
\references{
    H. Sun and S. Wang (2012)
    \emph{Penalized Logistic Regression for High-dimensional DNA Methylation
    Data with Case-Control Studies}, Bioinformatics 28(10), 1368--1375 \cr
    \cr
    H. Sun and S. Wang (2012)
    \emph{Network-based Regularization for Matched Case-Control
    Analysis of High-dimensional DNA Methylation Data}, Statistics in Medicine
    32(12), 2127--2139 
    }
\author{Hokeun Sun <hsun@pusan.ac.kr>}
\examples{
    n<-200
    p<-1000
    x<-matrix(rnorm(n*p),n,p)
    
    # one-to-one matched set
    y<-c(rep(0,n/2),rep(1,n/2))
    st<-rep(seq(n/2),2)
    
    # one-to-four matched set
    y<-c(rep(0,4*n/5),rep(1,n/5))
    st<-c(rep(seq(n/5),rep(4,n/5)),rep(seq(n/5),1))
     
    # a total of 100 groups each of which consists of 5, 10, or 20 members 
    gr<-c(rep(5,40),rep(10,40),rep(20,20))
    
    # an example of adjacency matrix
    adjm<-cov(x)
    diag(adjm)<-0
    adjm[abs(adjm)<=0.3]<-0
    adjm[abs(adjm)>0.3]<-1
    
    # an example of signs of coefficients
    sg<-sign(rnorm(p))
    
    # Lasso
    g1<-sel.pclogit(x,y,st)
    
    # Elastic-net
    g2<-sel.pclogit(x,y,st,alpha=0.1)
    
    # Ring network of grouped covariates
    g3<-sel.pclogit(x,y,st,alpha=0.1,group=gr)
    
    # Fully connected network of grouped covariates
    g4<-sel.pclogit(x,y,st,alpha=0.1,group=gr,type="fcon")
    
    # Graph-constrained covariates
    g5<-sel.pclogit(x,y,st,alpha=0.1,group=adjm)
    
    # Adaptive graph-constrained covariates
    g6<-sel.pclogit(x,y,st,alpha=0.1,group=adjm,sgnc=sg)
}
   

