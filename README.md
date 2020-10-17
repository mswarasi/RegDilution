
This repository contains R and FORTRAN programs associated with the article, “Group testing regression models with dilution submodels.”

We provide two main funcitons that work for logistic regression and any dilution submodels.

1. Dilution.Obs.MPT.R – main R program to fit regression models with initial pooled responses.

2. Dilution.EM.Dorfman.R – main R program to fit regression models with Dorfman two-stage pooling data.

Note: There are 12 compiled FORTRAN subroutines and sources codes that work under the main R programs as a means of improving computing efficiency. Documentation and user instructions are provided within the R programs.

Reference:

Warasi, M., McMahan, C., Tebbs, J., and Bilder, C. (2017). Group testing regression models with dilution submodels. Statistics in Medicine, 36, 4860-4872.


################### Simulation with master pool testing ###################

## Specify the working directory
setwd(dir = "C:\\programs")

## Import source files
source("SupportPrograms.txt")
source("MultistageHierarchicalData.txt")
source("TwoStageArrayData.txt")

# Specify the regression link function g(.). This program works with any g, not limited to the logit link we used in the paper.
g <- function(t){
	res <- 1/(1+exp(-t))
	return(res)	
}

# Specify the dilution submodel h(.,.,.). In this example, h is derived from the logistic CDF. The code works with any submodel, h. The hessian matrix is calculated by the optim(). 
h <- function(gamma,k,cj,Se){
	num <- exp(log(Se/(1-Se))-exp(gamma)*(cj-k)/cj)
	res <- num/(1+num)
	return(res)
}


## Results from the dilution model.
res[[1]]
$par
[1] -3.030533  1.849901  1.307539  0.963596

$value
[1] 413.2407

$counts
function gradient 
     127       NA 

$convergence
[1] 0

$message
NULL

$hessian
          [,1]      [,2]      [,3]      [,4]
[1,] 190.06302 115.91705  32.89459 -66.32533
[2,] 115.91705  93.28644  16.82929 -45.52714
[3,]  32.89459  16.82929  14.94246 -12.61639
[4,] -66.32533 -45.52714 -12.61639  28.41934

$p.value
[1] 0

## Results from the constant Se model in Vansteelandt, Goetghebeur, and Verstraeten (2000).
res[[2]]
$par
[1] -3.089595  1.729393  1.222612

$value
[1] 413.5733

$counts
function gradient 
      62       NA 

$convergence
[1] 0

$message
NULL

$hessian
          [,1]     [,2]     [,3]
[1,] 222.28961 139.7245 39.60555
[2,] 139.72445 116.3728 21.54830
[3,]  39.60555  21.5483 18.08751




################### Simulation with Dorfman testing ###################

## Specify the regression link function 'g'. 
g <- function(t){
	res <- 1/(1+exp(-t))
	return(res)	
}

## Specify the dilution submodel 'h' and its 1st and 2nd derivatives, dh and d2h. The code works with any submodel h:
h <- function(gamma,k,cj,Se){
	num <- exp(log(Se/(1-Se))-exp(gamma)*(cj-k)/cj)
	res <- num/(1+num)
	return(res)
}

dh <- function(gamma,k,cj,Se){
	H <- h(gamma,k,cj,Se)
	res <- H*(1-H)*exp(gamma)*(k-cj)/cj
	return(res)
}

d2h <- function(gamma,k,cj,Se){
	H <- h(gamma,k,cj,Se)
	DH <- dh(gamma,k,cj,Se)
	res <- DH*( DH/H - DH/(1-H) + 1 )
	return(res)
}

## Specify the true parameter values for simulation.
lambda.t <- 3.8		     	
gamma.t <- log(lambda.t) 	  
beta.t <- c(-3,2,1)		 	# Regression parameter.	
theta.t <- c(beta.t,gamma.t)  	        # Parameters to estimate.
Se.t <- 0.99	 			# Individual testing sensitivity.
Sp.t <- 0.99		 		# Individual testing specificity.

## Specify the vector of pool sizes.
c.vec <- c(rep(5,1000))		        
N <- sum(c.vec)                

## Set a seed for reproducible results, and simulate the covariate values.
set.seed(123)
x1 <- rnorm(N,0,.75)		        
x2 <- rbinom(N,1,0.1)		       
X <- cbind(1,x1,x2)		        

## Simulate group testing data.
gt <- GT.data.Dorfman(theta.t,c.vec,Se.t,Sp.t,X)
Z <- gt$pool.resp
Y <- gt$indv.resp

## Fit the model with the starting value: param0 = (beta0,gamma0). One may find beta0 using the constant Se model as shown below. A reasonable choice for gamma0 is 0. 
J <- length(c.vec)                   
id <- cumsum(c(0,c.vec))              
beta0 <- optim(beta.t,Dorfman.const.loglik,Z=Z,Y=Y,X=X,cvec=c.vec,Sp=Sp.t,
Se=Se.t,J=J,id=id)$par
param0 <- c(beta0,0)

res <- Dil.EM.Dorfman(param0=param0,X=X,Z=Z,Y=Y,cvec=c.vec,Se=Se.t,
		      Sp=Sp.t,method="Nelder-Mead",hessian=TRUE,emtol=10^(-6),
		      emmaxit=500,GI=12000,burn=200,interval=c(-10000,10),
		      est.con=TRUE)

## Results from the dilution model.
res[[1]]
$par
[1] -3.1069268  2.0167373  0.9740853  0.9870223

$hessian
          [,1]      [,2]     [,3]      [,4]
[1,] 302.07759 193.81718 45.23482 -58.33001
[2,] 193.81718 228.29168 21.76993 -45.14087
[3,]  45.23482  21.76993 46.06786  -9.74376
[4,] -58.33001 -45.14087 -9.74376  16.88521

$p.value
[1] 0.1284383

Results from the constant Se model with Dorfman testing
in Zhang, Bilder, and Tebbs (2013)

res[[2]]
$par
[1] -3.1669524  1.9796215  0.9573082

$hessian
          [,1]      [,2]     [,3]
[1,] 318.52970 212.12085 48.49853
[2,] 212.12085 253.01413 25.02662
[3,]  48.49853  25.02662 48.56958
