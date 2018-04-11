\name{pclogit}
\alias{pclogit}
\title{penalized conditional (unconditional) logistic regression for 
    grouped or graph-constrained variables}
\description{
    Fit a regularization path of conditional (unconditional) logistic
    regression model for a matched (unmatched) case-control response 
    at a grid of values for regularization parameter lambda.\cr 
    \cr
    When predictors are correlated witin either a group or a network graph, 
    Laplacian matrix can be imposed into the regularization procedure to
    efficiently select relevant variables.
    }
\usage{pclogit(x,y,stra=NULL,alpha=1.0,nlam=100,lambda=NULL,eps=NULL,
    maxit=100000,maxp=p,thre=1e-6,group=NULL,type=c("ring","fcon"),
    sgnc=NULL)
    }
\arguments{
    \item{x}{The design matrix (\eqn{n} x \eqn{p}) without an intercept.
        \code{pclogit} standardizes the data by default, but the coefficients
        are always returned on the original scale.}
    \item{y}{The response variable coded as 1 for cases and 0 for the 
        matched controls.}
    \item{stra}{A vector of consecutive integers indicating the stratum of
        each observation. Each stratum must have exact one case and at least
        one control. If not specified, \code{pclogit} fits an ordinary 
        logistic regression.}
    \item{alpha}{The penalty mixing parameter with \eqn{0\le\alpha\le 1} and
        default is 1. See details.}
    \item{nlam}{The number of \code{lambda} values and default is 100.}
    \item{lambda}{A user supplied sequence of \code{lambda} values. Typically,
        this is left unspecified, and the program automatically computes its
        own \code{lambda} sequence based on \code{nlam} and \code{eps}.}
    \item{eps}{The smallest value for \code{lambda} as a fraction
        of \code{lambda.max}. The value should be greater than \code{1E-5}.
        The default is .0001 if \eqn{n > p} and .05 if \eqn{n \le p}.}
    \item{maxit}{Maximum number of passes over the data for all \code{lambda}
        values, and default is 10^5.}
    \item{maxp}{Limit the maximum number of variables ever to be nonzero.}
    \item{thre}{Convergence threshold for coordinate descent algorithm.
        The default value is \code{1E-6}.}
    \item{group}{Either an integer vector of group sizes or a symmetric 
        adjacency matrix. \code{group} describes either grouped or graph 
        structure of predictors \code{x}. If no information between predictors
        is available, i.e., not specified, the \code{pclogit} performs an 
        elastic-net regularization procedure. See details.}
    \item{type}{A type of grouping network when \code{group} is defined as 
        a vector of group sizes. "\code{ring}" and "\code{fcon}" represent
        a ring and fully connected network, respectively. Default is 
        "\code{ring}". See details.}     
    \item{sgnc}{Signs of regression coefficients. This can be provided only
        if \code{group} is specified as either a list of a group size or an 
        adjacency matrix. The estimated signs of ridge regression for 
        \eqn{n\le p} or ordinary regression for \eqn{n > p} can be used
        for an adaptive network-based regularization procedure.
        See details.}
        }
\details{
       The penalty function of \code{pclogit} is defined as
       \deqn{\alpha||\beta||_1+(1-\alpha)(\beta^{T}S^{T}LS\beta)/2,}
       where \eqn{S} is a \eqn{p} dimensional diagonal matrix with 
       estimated signs of regression coefficients on its diagonal entries,
       and \eqn{L} is a Laplacian matrix describing a graph structure of
       covariates. This penalty is equivalent to the Lasso penalty if
       \code{alpha=1}. When \code{group} and \code{sgnc} are not defined,
       \eqn{L} and \eqn{S} in the penalty function are replaced by
       an identity matrix, respectively. In this case, \code{pclogit}
       performs an elastic-net regularization procedure since the second
       term of the penalty simply reduces to the \eqn{l_2} norm of 
       \eqn{\beta}. \cr
       \cr
       If group sizes of predictors are listed in \code{group}, it is assumed
       that all variables of the same groups are linked with each other like
       a ring or a fully connected network. In this case, the Laplacian matrix
       forms a block-wise diagonal matrix. The signs of regression coefficients 
       \code{sgnc} can provide more accuarate estimates in case some variables 
       either in the same group or linked with each other have different signs
       of thier regression coefficients, where the coefficients are not expected
       to be loccally smooth. 
       }
\value{
    \item{b0}{Intercept sequence of length of \code{lambda}. This is present
       only if an ordinary logistic regression is fit, i.e., \code{stra} was not
       defined.}
    \item{strata}{The strata of obsevations if \code{stra} was defined.}
    \item{beta}{The coefficient matrix with a dimension (\eqn{p} x \code{nlam}),
       stored in sparse column format (\code{"CsparseMatrix"})}
    \item{lambda}{The actual sequence of \code{lambda} values used}
    \item{df}{The number of nonzero coefficients for each value of 
      \code{lambda}}
    \item{nobs}{The number of observations, \eqn{n}}
    \item{alpha}{The value of \code{alpha} used}
    \item{iterations}{Total passes over the data summed over all lambda values}
     \item{jerr}{The error flag, for warnings and errors (largely for internal
      debugging)}
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
    g1<-pclogit(x,y,st)
    
    # Elastic-net
    g2<-pclogit(x,y,st,alpha=0.1)
    
    # Ring network of grouped covariates
    g3<-pclogit(x,y,st,alpha=0.1,group=gr)
    
    # Fully connected network of grouped covariates
    g4<-pclogit(x,y,st,alpha=0.1,group=gr,type="fcon")
    
    # Graph-constrained covariates
    g5<-pclogit(x,y,st,alpha=0.1,group=adjm)
    
    # Adaptive graph-constrained covariates
    g6<-pclogit(x,y,st,alpha=0.1,group=adjm,sgnc=sg)
}
    



