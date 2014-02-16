## ando.R -- Testing scripts for Ando & Zhang's structural learning algorithm
##
## by Artem Sokolov

test1.data <- function()
  {
    p <- 100	## Raw space dimensionality
    h <- 10	## Shared space dimensionality
    l <- 20	## Number of problems
    n <- 200	## Number of samples in each problem

    ## Generate a (random) common linear map Theta
    ## Theta %*% t(Theta) must be an h-by-h orthonormal matrix, a constraint
    ##   satisfied by the SVD bases
    M <- matrix( runif(p*h), h, p )
    Theta <- t(svd(M)$v)

    ## Generate the "true" w and v weights for the l prediction problems
    W <- matrix( runif( p*l, -1, 1 ), p, l )
    V <- matrix( runif( h*l, -1, 1 ), h, l )
    
    ## Generate the l prediction problems
    X <- list()
    y <- list()
    for( i in 1:l )
      {
        ## Use a Gaussian "white" ball as the input data
        m <- runif( p, -l, l )		## Mean of the Gaussian ball
        X[[i]] <- matrix( rnorm( n*p, mean=m ), n, p )

        ## Compute the response
        w <- W[,i]
        v <- V[,i]
        y[[i]] <- X[[i]] %*% W[,i] + X[[i]] %*% t(Theta) %*% V[,i]
      }

    list( X.list = X, y.list = y, Theta = Theta, W = W, V = V )
  }

## Evaluates the objective for a single prediction problem at (w.hat, v.hat, Theta.hat)
f.obj1 <- function( X, y, w.hat, v.hat, Theta.hat, lambda=1 )
  {
    ## Compute the predictions using parameter estimates
    z <- X %*% w.hat + X %*% t(Theta.hat) %*% v.hat

    ## Compute the loss (MSE)
    loss <- mean( (y - z)^2 )

    ## Compute the regularization term
    reg <- sum( w.hat^2 )

    ## The joint objective value is loss + lambda * regularizer
    loss + lambda * reg
  }

## Evaluates the joint objective function at (W.hat, V.hat, Theta.hat)
f.obj <- function( X.list, y.list, W.hat, V.hat, Theta.hat, lambda=1 )
  {
    l <- length(X.list)

    ## Compute problem-specific objective values
    val <- rep( 0, l )
    for( i in 1:l )
      val[i] <- f.obj1( X.list[[i]], y.list[[i]], W.hat[,i], V.hat[,i],
                       Theta.hat, lambda )
    
    ## Check the orthonormality constraint
    eiv <- eigen( Theta.hat %*% t( Theta.hat ) )$values
    cat( "The orthonormality constraint is violated by", max(abs(eiv - 1)), "\n" )

    sum(val)
  }

## Perturbation experiment
pert.sol <- function( X, y, w.hat, v.hat, Theta.hat, lambda=1, tol=1e-5 )
  {
    f0 <- f.obj1( X, y, w.hat, v.hat, Theta.hat, lambda )

    ## Perturb w
    for( j in 1:length(w.hat) )
      {
        ## Positive perturbation
        w1 <- w.hat; w1[j] <- w.hat[j] + tol
        f1 <- f.obj1( X, y, w1, v.hat, Theta.hat, lambda )
        if( f0 - f1 > tol*abs(f1) )
          cat( "Perturbing w[", j, "] in the positive direction yields the following improvement:", f0-f1, "\n" )

        ## Negative perturbation
        w2 <- w.hat; w2[j] <- w.hat[j] - tol
        f2 <- f.obj1( X, y, w2, v.hat, Theta.hat, lambda )
        if( f0 - f2 > tol*abs(f2) )
          cat( "Perturbing w[", j, "] in the negative direction yields the following improvement:", f0-f2, "\n" )
      }

    ## Perturb v
    for( j in 1:length(v.hat) )
      {
        ## Positive perturbation
        v1 <- v.hat; v1[j] <- v.hat[j] + tol
        f1 <- f.obj1( X, y, w.hat, v1, Theta.hat, lambda )
        if( f0 - f1 > tol )
          cat( "Perturbing v[", j, "] in the positive direction yields the following improvement:", f0-f1, "\n" )
        
        ## Negative perturbation
        v2 <- v.hat; v2[j] <- v.hat[j] - tol
        f2 <- f.obj1( X, y, w.hat, v2, Theta.hat, lambda )
        if( f0 - f2 > tol )
          cat( "Perturbing v[", j, "] in the negative direction yields the following improvement:", f0-f2, "\n" )
      }
  }

main <- function()
  {
    mydata <- test1.data()
	source("multitask.R")
	h <- 10
	output <- ando_test_output(mydata, h)
	W.hat <- output$W.hat
	V.hat <- output$V.hat
	Theta.hat <- output$Theta.hat

	#W.hat <- w_min(

    ## Train a model on (mydata$X.list[[i]], mydata$y.list[[i]]) prediction problems
    ## W.hat <- ...
    ## V.hat <- ...
    ## Theta.hat <- ...

    ## Evaluate the objective at the solution
    val <- f.obj( mydata$X.list, mydata$y.list, W.hat, V.hat, Theta.hat, 1 )
    cat( "Joint objective value =", val, "\n" )

    ## Perform the perturbation experiments
    for( i in 1:length(mydata$X.list) )
      {
        cat( "=== Testing problem", i, " ===\n" )
        X <- mydata$X.list[[i]]
        y <- mydata$y.list[[i]]
        w.hat <- W.hat[,i]
        v.hat <- V.hat[,i]
        pert.sol( X, y, w.hat, v.hat, Theta.hat, 1 )
      }
  }
