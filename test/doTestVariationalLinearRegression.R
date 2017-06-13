
source("variationalLinearRegression.R")
source("../../math/tr.R")
source("../../math/l2Norm2.R")
source("../../math/logdet.R")

processAll <- function() {
    set.seed(0)
    wlims <- c(-5, 5)
    ylims <- c(-11, 11)

    D <- 100
    regressorsNames <- sprintf("%d", 1:D)
    N <- 150
    N_test <- 50
    w <- rnorm(D)
    x <- matrix(rnorm(N*D)-.5, nrow=N)
    colnames(x) <- regressorsNames
    x_test <- matrix(rnorm(N_test*D)-.5, nrow=N_test)
    colnames(x_test) <- regressorsNames
    y <- x%*%w+rnorm(N)
    y_test <- x_test%*%w+rnorm(N_test)


    resVLM <- variationalLinearRegression(t=y, phi=x, 
                                               a0=1e-2, b0=1e-4,
                                               c0=1e-2, d0=1e-4, 
                                               maxIter=500,
                                               convergenceTol=1e-4)
    w_VB <- resVLM$mN
    V_VB <- resVLM$sN
    y_VB <- x%*%w_VB
    y_test_VB <- x_test%*%w_VB

    res_LM <- lm(y~x, data=list(y=y, x=I(x)))
    w_ML <- coefficients(res_LM)
    wint_ML <- confint(res_LM)
    resPredict_LM <- predict(res_LM, interval="predict")
    y_ML <- resPredict_LM[,1]
    y_ML_confint <- resPredict_LM[, c(2,3)]
    resPredict_LM_test <- predict(res_LM, newdata=list(x=I(x_test)), 
                                          interval="predict")
    y_test_ML <- resPredict_LM_test[,1]
    y_test_ML_confint <- resPredict_LM_test[, c(2,3)]

    show(sprintf("Training set MSE: ML = %f, VB = %f", 
                 mean((y-y_ML)^2), mean((y-y_VB)^2)))
    show(sprintf("Test     set MSE: ML = %f, VB = %f", 
                 mean((y_test-y_test_ML)^2), mean((y_test-y_test_VB)^2)))
}

processAll()

rm(processAll)
