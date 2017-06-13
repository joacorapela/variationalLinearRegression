
# t: dependent variable (length N)
# phi: independent variables (dim N x M)
# a0, b0: initial value for coefs precision
# c0, d0: initial value for noise precision
#
# Note: both t and phi should be standarized

variationalLinearRegression <- function(t, phi, a0, b0, c0, d0, 
                                           maxIter,
                                           convergenceTol=1e-5) {
    computeLowerBound <- function() {
        ePyGwt <- N/2*(digamma(aN)-log(bN)-log(2*pi))-
                   .5*(tr(phiTPhi%*%vN)+aN/bN*l2Norm2(t-phi%*%mN))
        ePwtGa <-M/2*(digamma(aN)-log(bN)+digamma(cN)-log(dN)-log(2*pi))-
                  .5*cN/dN*(aN/bN*l2Norm2(mN)+tr(vN))-
                  lgamma(a0)+a0*log(b0)+
                  (a0-1)*(digamma(aN)-log(bN))-b0*aN/bN
        ePa <- -lgamma(c0)+d0*log(c0)+
                (c0-1)*(digamma(cN)-log(dN))-d0*cN/dN
        eQwt <- M/2*(digamma(aN)-log(bN)-log(2*pi)-1)-.5*logdet(vN)-
                 lgamma(aN)+aN*log(bN)+
                 (aN-1)*(digamma(aN)-log(bN))-aN
        eQa <- -lgamma(cN)+(cN-1)*digamma(cN)+log(dN)-cN
        lowerBound <- ePyGwt+ePwtGa+ePa-eQwt-eQa
        return(lowerBound)
    }

    N <- nrow(phi)
    M <- ncol(phi)
    phiT <- t(phi)
    phiTPhi <- phiT%*%phi
    idM <- diag(M)
    lowerBound <- -.Machine$double.xmax
    converged <- FALSE
    lowerBounds <- c()
    eAlpha <- c0/d0
    aN <- a0; bN <- b0; cN <- c0; dN <- d0
    for(i in 1:maxIter) {
        show(sprintf("Iteration %d: eAlpha=%f, lowerBound=%f", 
                     i, eAlpha, lowerBound))
        vN <- solve(eAlpha*idM+phiTPhi)
        mN <- as.vector(vN%*%phiT%*%t)
        lowerBoundOld <- lowerBound
        lowerBound <- computeLowerBound()
        if(lowerBoundOld>lowerBound) {
            warning(sprintf("Lower bound decreased from %f to %f", 
                            lowerBoundOld, lowerBound))
        }
        if(lowerBound-lowerBoundOld<abs(convergenceTol*lowerBoundOld)) {
            converged <- TRUE
            break
        }
        lowerBounds <- c(lowerBounds, lowerBound)
        aN <- a0+N/2
        bN <- b0+.5*(l2Norm2(t-phi%*%mN)+eAlpha*l2Norm2(mN))
        eTauWTW2 <- tr(vN)+l2Norm2(mN)*aN/bN

        cN <- c0+M/2
        dN <- d0+eTauWTW2/2
        eAlpha <- cN/dN
    }
    return(list(mN=mN, vN=vN, aN=aN, bN=bN, cN=cN, dN=dN,
                       lowerBounds=lowerBounds, converged=converged))
}

