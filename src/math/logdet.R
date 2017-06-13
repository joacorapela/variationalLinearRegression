logdet <- function(m) {
    return(2*sum(log(diag(chol(m)))))
}

