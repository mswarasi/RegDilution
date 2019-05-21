
#########################################################################
#     R functions to find MLEs of the parameter 'theta' for master      #   
#     pool testing (MPT) corresponding to the manuscript titled,        #
#     "Group testing regression models with dilution submodels."        #
#                           Date: 12/21/2016.                           #
#########################################################################

            
#***********************************************************************#
# Log-likelihood function for the dilution model. 
# Specify the directory for "Dil_MPT_LogLik.dll".
#
MPT.Dil.loglik <- function(theta,Z,X,cvec,Se,Sp){
	dyn.load("Dil_MPT_LogLik.dll")
	beta <- theta[-length(theta)]
	gamma <- theta[length(theta)]
	id <- cumsum(c(0,cvec))
	J <- length(cvec)
	p <- g(X%*%beta)
	temp <- -9
	cmax <- as.integer(max(cvec))
	psz <- cvec[!duplicated(cvec)]
	se_mat <- matrix(-9,cmax,cmax)
	for(j in 1:length(psz)){
		se_mat[psz[j],1:psz[j]] <- h(gamma,1:psz[j],psz[j],Se)
	}
	res <- -.C("dil_mpt_loglik",temp,p,gamma,as.integer(Z),Sp,
	    as.integer(cvec),as.integer(J),as.integer(N),se_mat,cmax)[[1]]
	dyn.unload("Dil_MPT_LogLik.dll")
	return(res)
}

#***********************************************************************#
# Log-likelihood function of the constant Se model in 
# Vansteelandt, Goetghebeur, and Verstraeten (2000). 
# This is needed to perform the proposed likelihood 
# ratio test to detect dilution. 
# Specify the directory for "Van_LogLik.dll".
#               
MPT.Van.loglik <- function(beta,Z,X,cvec,Se,Sp){
	dyn.load("Van_LogLik.dll")
	J <- length(cvec)
	p <- g(X%*%beta)
	N <-length(p)
	temp <- -9
	res <- .C("van_loglik",temp,p,as.integer(Z),Sp,Se,as.integer(cvec),
			   as.integer(J),as.integer(N))[[1]]
	dyn.unload("Van_LogLik.dll")
	return(-res)
}


#***********************************************************************#
# Main function of the proposed model for master pool testing (MPT).
#
#
# Arguments:
# param0          = Initial value of the parameter 'theta'.
# Z               = J by 1 vector of master pooled responses.
# X               = Design matrix: N by length(beta) matrix.
# cvec            = J by 1 vector of pool sizes cj.
# Se              = Individual testing sensitivity.
# Sp              = Individual testing specificity.
# method          = Optimization method to pass to 'optim'. We found that 
#                   + 'Nelder-Mead' works the best for this model.
# hessian         = Calculates hessian matrix if TRUE.
# lower           = Lower bound for the "L-BFGS-B" method; see 'optim'.
# upper           = Upper bound for the "L-BFGS-B" method; see 'optim'.
# control         = A list of control parameters; see 'optim'.
# lower.con       = Lower bound for the constant Se model when
#                   + implementing the "L-BFGS-B" method; see 'optim'.
# upper.con       = Upper bound for the constant Se model when
#                   + implementing the "L-BFGS-B" method; see 'optim'.
# control.con     = A list of control parameters for the constant
#                   + Se model; see 'optim'.
# est.con         = Calculates results using the constant Se model
#                   + in Zhang, Bilder, and Tebbs (2013), if TRUE.
#
#
#
# Value: A list of two lists. The first list consists of the results
# 		 from the dilution model, and the second one consists of the 
# 		 results from the constant Se model when 'est.con' is true.
# 	         The output values are labeled identically as labeled for 'optim'. 
# 		 i.e, the list components are: par, hessian, convergence, etc.
#		 Additionally, the first list (dilution result) provides
# 		 p-value, labeled as 'p.value', for the proposed likelihood
# 		 ratio test to detect dilution.
#
#
# Details:
# (1) To improve computational efficiency, the code uses several
#     + FORTRAN 77 subroutines through the dynamic-link library (DLL).
#     + These subroutines only work in a 64-bit R package.
#
# (2) Specify DLLs' source directory.
#
# (3) Specify R's working directory if needed.
#
# (4) The program works with any link function 'g' and any submodel 'h'.
#
# (5) The optional arguments, such as, lower, upper, control, etc., 
#     + are identical to the arguments in optim; i.e., these arguments 
#     + can be specified similarly as in optim.
#
#
Dil.Obs.MPT <- function(param0,Z,X,cvec,Se,Sp,method="Nelder-Mead"
                        ,hessian=FALSE,lower=-Inf,upper=Inf,control=list(),
                        lower.con=-Inf,upper.con=Inf,control.con=list(),
                        est.con=FALSE){

#	Fit the dilution model. 
#					 
	dil <- optim(par=param0,fn=MPT.Dil.loglik,method=method,lower=lower,
                        upper=upper,control=control,hessian=hessian,Z=Z,X=X,
                        cvec=cvec,Se=Se,Sp=Sp)
				 
#	Value of the log-likelihood function for the dilution model.
#					 
	ll.dil <- (-1)*dil$value

#	Fit the constant Se model. 
#
	van <- optim(par=param0[-length(param0)],fn=MPT.Van.loglik,method=method,
                         lower=lower.con,upper=upper.con,control=control.con,
                         hessian=est.con,Z=Z,X=X,cvec=cvec,Se=Se,Sp=Sp)

#	Value of the log-likelihood function for the constant Se model.
#
	ll.van <- (-1)*van$value

#	Test statistic and p-value of the likelihood ratio test.
#	
	T.LR <- 2*ll.dil - 2*ll.van
	p.val <- 0.5*pchisq(T.LR,0,lower.tail=FALSE)
                         + 0.5*pchisq(T.LR,1,lower.tail=FALSE)

#	Return all results in a list.
#		
	res.dil <- c(dil,"p.value"=p.val)	
	res.con <- van
	if(est.con == FALSE) res.con <- NULL
	return(list(dilution=res.dil,constant=res.con))
}


###########################################################################
#                                                                         #                                                      
#    For illustration, simulate group testing data and fit the models.    #
#                                                                         #
###########################################################################

# The function 'GT.data.Dorfman' simulates group testing
# data using the configuration described in the paper.
#
# Arguments: 
# theta 	  = True parameter value; see examples below.
# cvec 		  = J by 1 vector of pool sizes cj; see examples below. 
# Se		  = Individual testing sensitivity.
# Sp 		  = Individual testing specificity.	
# X 		  = Design matrix: N by length(beta) matrix.
#
#
# Value: 
# Z               = J by 1 vector of master pooled responses.
#
GT.data.MPT <- function(theta,cvec,Se,Sp,X){
	beta <- theta[-length(theta)]
	gamma <- theta[length(theta)]
	J <- length(cvec)
	N <- sum(cvec)
	id <- cumsum(c(0,cvec))
	p <- g(X%*%beta)		# Individual disease probability.
	Ytil <- rbinom(N,1,p)	        # Individual true statuses.
	Z <- rep(-9,J)			# Vector of pooled testing results.
	for(j in 1:J){
		Ytemp <- Ytil[(id[j]+1):(id[j+1])]
		SY <- sum(Ytemp)+1
		p.tpk <- c( 1-Sp, h(gamma,1:cvec[j],cvec[j],Se) )
		p.pos <- p.tpk[SY]
		Z[j] <- rbinom(1,1,p.pos)
	}
	return(list("pool.resp"=Z))
}

#---------------------- Regression model g(.) ---------------------------#
#
# Specify the regression link function g. This program works 
# with any g, not limited to the logit link we used in the paper.
#
g <- function(t){
	res <- 1/(1+exp(-t))
	return(res)	
}

#---------------------- Dilution submodel h(.,.,.) -----------------------#
#
# Specify the dilution submodel 'h'. As in our simulation study,
# we use the submodel h derived from the logistic CDF. 
# 
# Note, the code works with any submodel h. The hessian matrix
# is calculated by optim using numerical derivative. 
# 
# Also note, exponential transformation has been implemented to 
# avoid constraint optimization of the dilution parameter lambda. 
# i.e., we estimate gamma, where lambda = exp(gamma), and 
# -Inf < gamma < Inf. 
#
h <- function(gamma,k,cj,Se){
	num <- exp(log(Se/(1-Se))-exp(gamma)*(cj-k)/cj)
	res <- num/(1+num)
	return(res)
}



#*************************************************************************#
# Specify the true parameter values for simulation.
#
lambda.t <- 3.8		     	          # Dilution parameter lambda.
gamma.t <- log(lambda.t) 	          # Transformed version of lambda.
beta.t <- c(-3,2,1)                       # Regression parameter.		
theta.t <- c(beta.t,gamma.t)              # Parameter to estimate.
Se.t <- 0.99                              # Individual testing sensitivity.
Sp.t <- 0.99                              # Individual testing specificity.

# Specify a vector of pool sizes, either constant or variable.
#
# c.vec <- rep(5,1000)                    # Configuration 1.
# c.vec <- rep(10,500)                    # Configuration 2.
c.vec <- c(rep(5,334),rep(10,333))        # Configuration 3.

N <- sum(c.vec)                           # Sample size. 

# Set a seed for reproducible results.
#
set.seed(123)
x1 <- rnorm(N,0,.75)		          # x1 from Normal(mu=0,sigma=0.75).
x2 <- rbinom(N,1,0.1)		          # x2 from Bernoulli(p = 0.1).
X <- cbind(1,x1,x2)                       # Design matrix.


# Set R's working directory.
#
setwd(dir = "C:\\programs")

# Simulate group testing data.
#                                          
gt <- GT.data.MPT(theta.t,c.vec,Se.t,Sp.t,X)
Z <- gt$pool.resp


# Fit the model with the starting value: param0 = (beta0,gamma0).
# One may find beta0 using the constant Se model as shown below. 
# A reasonable choice for gamma0 is 0. 
#
#
beta0 <- optim(beta.t,MPT.Van.loglik,Z=Z,X=X,cvec=c.vec,Se=Se.t,Sp=Sp.t)$par
param0 <- c(beta0,0)

res <- Dil.Obs.MPT(param0=param0,Z=Z,X=X,cvec=c.vec,Se=Se.t,
                   Sp=Sp.t,method="Nelder-Mead",hessian=TRUE,est.con=TRUE)

#
# Presented are estimation results from the data generated using the seed 123.
# Choose a different seed for different group testing data and results.	
#
								
# Results from the dilution model.
#
res[[1]]

# Results from the constant Se model in Vansteelandt, 
# Goetghebeur, and Verstraeten (2000).
#
res[[2]]

