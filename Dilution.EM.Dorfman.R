
#########################################################################
#   R functions to implement the EM algorithm for Dorfman testing (DT)	#
#   corresponding to the manuscript titled, "Group testing regression   #
#   models with dilution submodels."					#
#                         Date: 12/21/2016.			    	#	
#########################################################################


#***********************************************************************#
# Function to produce outcomes of flipping coins.
# This function was written by anonymous authors.
#
tosscoin <- function(times){
    temp <- list()
    for (i in 1:times){
        temp[[i]] <- c(0,1)
    }
    res <- expand.grid(temp, KEEP.OUT.ATTRS = FALSE)
return(as.matrix(res))
}

#***********************************************************************#
# E-step associated with estimating the regression parameter 'beta'.
# Specify directories for "eyijtilz1.dll" and "eyijtilz0.dll"
#
E.beta <- function(beta,X,Z,Y,cvec,gamma,Se,Sp,J,id){
	dyn.load("eyijtilz1.dll")
	dyn.load("eyijtilz0.dll")
	p <- g(X%*%beta)
	eyf <- rep(-9,N)	
	for(j in 1:J){
		pj <- p[(id[j]+1):(id[j+1])]
		yj <- Y[(id[j]+1):(id[j+1])]
		eyj <- rep(-9,cvec[j])
		Sej <- h(gamma,1:cvec[j],cvec[j],Se)		
		if(Z[j]==1){
			if(cvec[j]==1){
				eyf[(id[j]+1):(id[j+1])] <- Se*pj/((1-Sp)*(1-pj)+Se*pj)
			}
			if(cvec[j]>1){	
				ssp <- tosscoin(cvec[j]-1)
				rw <- dim(ssp)[1]
				cl <- dim(ssp)[2]			
				eyf[(id[j]+1):(id[j+1])] <- .C("eyijtilz1",eyj,pj,
				  as.integer(cvec[j]),as.integer(yj),as.integer(ssp),
				  as.integer(rw),as.integer(cl),Se,Sp,Sej)[[1]]
			}
		}
		if(Z[j]==0){
		    eyf[(id[j]+1):(id[j+1])] <- .C("eyijtilz0",eyj,pj,
						    as.integer(cvec[j]),Sp,Sej)[[1]]
		}
	}
	dyn.unload("eyijtilz1.dll")
	dyn.unload("eyijtilz0.dll")	
	return(eyf)
}

#***********************************************************************#
# E-step associated with estimating the dilution parameter 'lambda'.
# Specify directories for "esyzj.dll" and "eszj.dll".
#
E.gamma <- function(gamma,X,Z,Y,cvec,beta,Se,Sp,J,id){
	dyn.load("esyzj.dll")
	dyn.load("eszj.dll")
	p <- g(X%*%beta)
	expt <- matrix(-9,J,max(cvec))
	for(j in 1:J){
		ezj <- rep(-9,cvec[j])
		Yj <- Y[(id[j]+1):(id[j+1])]
		pj <- p[(id[j]+1):(id[j+1])]
		Sej <- h(gamma,1:cvec[j],cvec[j],Se)
		
		if(Z[j]==1){
			if(cvec[j]==1){
			     expt[j,1:cvec[j]] <- Se*pj/( (1-Sp)*(1-pj) + Se*pj )
			}
			if(cvec[j]>1){
				ssp <- tosscoin(cvec[j]-1)
				rw <- 2^(cvec[j]-1)
				cl <- cvec[j]-1
				expt[j,1:cvec[j]] <- .C("esyzj",ezj,as.integer(Yj),
				    pj,Sp,as.integer(cvec[j]),as.integer(ssp),
				    as.integer(rw),as.integer(cl),Se,Sej)[[1]]
			}
		}
		if(Z[j]==0){                 
			expt[j,1:cvec[j]] <- .C("eszj",ezj,pj,Sp,
			  			as.integer(cvec[j]),Sej)[[1]]
		}
	}
	dyn.unload("esyzj.dll")
	dyn.unload("eszj.dll")
	return(expt)
}

#***********************************************************************#
# Objective function to be maximized with respect to 'beta'.
#
Q.beta <- function(beta,eyij){
	p <- g(X%*%beta)
	res <- sum(eyij*log(p)+(1-eyij)*log(1-p))  
	return(-res)
}

#***********************************************************************#
# Objective function to be maximized with respect to 'gamma'.
#
Q.gamma <- function(gamma,Z,expt,cvec,Se,J){
	res <- rep(-9,J)
	for(j in 1:J){
		Sej <- h(gamma,1:cvec[j],cvec[j],Se)
		if(Z[j] > 0){
			res[j] <- sum(expt[j,1:cvec[j]]*log(Sej))
		}else{
			res[j] <- sum(expt[j,1:cvec[j]]*log(1-Sej))
		}
	}
	return(-sum(res))
}

#***********************************************************************#
# Observed log-likelihood function for the dilution model.
# This is needed to perform the likelihood ratio test.
# Specify directories for "obsrez1.dll" and "obsrez0.dll".
#
Dorfman.Dil.loglik <- function(theta,Z,Y,X,cvec,Sp,Se,J,id,plen){
	dyn.load("obsrez1.dll")
	dyn.load("obsrez0.dll")
	beta <- theta[-plen]
	gamma <- theta[plen]
	p <- g(X%*%beta)
	res <- rep(-9,J)
	out <- -9
	for(j in 1:J){
		pj <- p[(id[j]+1):(id[j+1])]
		yj <- Y[(id[j]+1):(id[j+1])]
		Sej <- h(gamma,1:cvec[j],cvec[j],Se)
		if(Z[j] > 0){
		    if(cvec[j]==1){
				pz1 <- (1-Sp)*(1-pj)+Se*pj
			}
		    if(cvec[j] > 1){			
				ssp <- tosscoin(cvec[j]-1)
				rw <- 2^(cvec[j]-1)
				cl <- cvec[j] - 1
				pz1 <- .C("obsrez1",out,as.integer(yj),pj,Sp,
				 as.integer(cvec[j]),as.integer(ssp),as.integer(rw),
				 as.integer(cl),Se,Sej)[[1]]
			}
			res[j] <- log(pz1)
		}else{
			pz0 <- .C("obsrez0",out,pj,Sp,as.integer(cvec[j]),Sej)[[1]]
			res[j] <- log(pz0)
		}
	}
	dyn.unload("obsrez1.dll")
	dyn.unload("obsrez0.dll")	
	return(-sum(res))
}

#***********************************************************************#
# Observed log-likelihood function of the constant Se 
# model in Zhang, Bilder, and Tebbs (2013). This is needed 
# to perform the likelihood ratio test.
# Specify directories for "upz1revands.dll" and "upz0revands.dll".
#
Dorfman.const.loglik <- function(beta,Z,Y,X,cvec,Sp,Se,J,id){
	dyn.load("upz1revands.dll")
	dyn.load("upz0revands.dll")
	p <- g(X%*%beta)
	res <- rep(-9,J)
	out <- -9
	for(j in 1:J){
		pj <- p[(id[j]+1):(id[j+1])]
		yj <- Y[(id[j]+1):(id[j+1])]	
		if(Z[j] > 0){
		    if(cvec[j]==1){
				pz1 <- (1-Sp)*(1-pj)+Se*pj
			}
			if(cvec[j]>1){
				ssp <- tosscoin(cvec[j]-1)
				row <- 2^(cvec[j]-1)
				col <- cvec[j] - 1
				pz1 <- .C("upz1revands",out,as.integer(yj),pj,Sp,
					  as.integer(cvec[j]),as.integer(ssp),
					  as.integer(row),as.integer(col),Se)[[1]]
			}
			res[j] <- log(pz1)
		}else{
			pz0 <- .C("upz0revands",out,pj,Sp,as.integer(cvec[j]),Se)[[1]]
			res[j] <- log(pz0)
		}
	}
	dyn.unload("upz1revands.dll")
	dyn.unload("upz0revands.dll")	
	return(-sum(res))
}

#***********************************************************************#
# Gibbs sampler: Generates individuals' true statuses given 
# Dorfman decoding data. The Gibbs sampler is needed to 
# approximate the covariance matrix using Louis's method.
# Specify the directory for "gibbs_ytil.dll".
#
gibbs.ytil <- function(dataf,p,gamma,Se,Sp,N,cvec,cmax,GI,burn){
	dyn.load("gibbs_ytil.dll")
	mat <- matrix(-9,N,GI)
	seed <- sample(2.1*10^9,1)
	psz <- cvec[!duplicated(cvec)]
	se_mat <- matrix(-9,cmax,cmax)
	for(j in 1:length(psz)){
		se_mat[psz[j],1:psz[j]] <- h(gamma,1:psz[j],psz[j],Se)
	}
	gibbs.ytil <- .C("gibbs_ytil",as.integer(mat),as.integer(dataf),
				p,Se,Sp,as.integer(N),as.integer(cmax),as.integer(GI),
				as.integer(burn),as.integer(seed),se_mat)[[1]]
	dyn.unload("gibbs_ytil.dll")
	Ymat <- matrix(gibbs.ytil,N,GI)
	return(Ymat)
}

#***********************************************************************#
# This function is needed to calculate the varience-covariance 
# matrix using Louis's method.
# Specify the directory for "cov_dl_gen.dll".
#
cov.dL <- function(param,dataf,Z,X,Se,Sp,cvec,N,J,plen,GI,burn){
	dyn.load("cov_dl_glgt.dll")
	beta <- param[-plen]
	gamma <- param[plen]
	p <- g(X%*%beta)
	cmax <- as.integer(max(cvec))
	res <- matrix(-9,N,GI)
	seed <- sample(10^10,1)
	Ymat <- gibbs.ytil(dataf,p,gamma,Se,Sp,N,cvec,cmax,GI,burn)
	rw <- cl <- plen
	cv <- matrix(-9,rw,cl)
	psz <- cvec[!duplicated(cvec)]
	dse_mat <- se_mat <- matrix(-9,cmax,cmax)
	for(j in 1:length(psz)){
		se_mat[psz[j],1:psz[j]] <- h(gamma,1:psz[j],psz[j],Se)
		dse_mat[psz[j],1:psz[j]] <- dh(gamma,1:psz[j],psz[j],Se)
	}
	res <- .C("cov_dl_glgt",cv,p,as.integer(Z),X,as.integer(N),
			as.integer(J),as.integer(GI),as.integer(cvec),
			as.integer(rw),as.integer(cl),as.integer(Ymat),
			as.integer(plen),se_mat,cmax,dse_mat)[[1]]
	dyn.unload("cov_dl_glgt.dll")
	return(res)
}

#***********************************************************************#
# This function is needed to calculate the varience-covariance 
# matrix using Louis's method.
#
DQsq.gamma <- function(gamma,ez,cvec,Se,J){
	res <- rep(-9,J)
	for(j in 1:J){
		cj <- cvec[j]
		H <- h(gamma,1:cj,cj,Se)
		DH <- dh(gamma,1:cj,cj,Se)
		D2H <- d2h(gamma,1:cj,cj,Se)
		temp <- Z[j]*(H*D2H-DH^2)/H^2-(1-Z[j])*((1-H)*D2H+DH^2)/(1-H)^2			
		res[j] <- sum(temp*ez[j,1:cj])
	}
	return(sum(res))
}

#***********************************************************************#
# This function is needed to calculate the varience-covariance 
# matrix using Louis's method.
#
cov.dQ <- function(param,Z,Y,X,cvec,Se,Sp,J,id,plen){
	blen <- plen-1
	beta <- param[-plen]
	gamma <- param[plen]
	p <- g(X%*%beta)	
	p.pcomp <- p*(1-p)
	dQ <- matrix(0,nrow=plen,ncol=plen)
	for(s in 1:blen){
		for(t in 1:blen){
			dQ[s,t] <- -sum(p.pcomp*X[ ,s]*X[ ,t])
		}
	}
	ez <- E.gamma(gamma,X,Z,Y,cvec,beta,Se,Sp,J,id)
	dQ[plen,plen] <- DQsq.gamma(gamma,ez,cvec,Se,J)
	return(dQ)
}

#***********************************************************************#
# This function is needed to run the Gibbs sampler
#
data.strc <- function(Z,Y,cvec,N,J,id){
	dtf <- matrix(0,nrow=N,ncol=4+max(cvec))
	for(j in 1:J){
		dtf[(id[j]+1):(id[j+1]),3] <- Z[j]
		dtf[(id[j]+1):(id[j+1]),4] <- cvec[j]	
		for(i in (id[j]+1):(id[j+1])){
			dtf[i,(4+1):(4+cvec[j])] <- (id[j]+1):id[j+1]
		}
	}
	dtf[ ,2] <- Y
	return(dtf)
}

#***********************************************************************#
# Main function to implement the EM algorithm for Dorfman testing.
#
#
# Arguments:
# param0 	  = Initial value of the parameter 'theta'.
# X 		  = Design matrix: N by length(beta) matrix.
# Z 		  = J by 1 vector of master pooled responses.
# Y 		  = N by 1 vector of individual retest results;
#		    + an individual's test response is 0 if not retested.
# cvec 		  = J by 1 vector of pool sizes cj.
# Se		  = Individual testing sensitivity.
# Sp		  = Individual testing specificity.
# method  	  = Optimization method to pass to 'optim'. We found that 
#		    + 'Nelder-Mead' works the best for this model.
# hessian 	  = Calculates hessian matrix if TRUE.
# emtol		  = Convergence criterion for the EM algorithm.
# emmaxit	  = Maximum number of EM iterates.
# GI 		  = Gibbs iterates after the burn-in period 'burn'. 
#		    + The Gibbs's sampler is needed to estimate the 
# 	            + covariance matrix using Louis's method.
# burn		  = Burn-in period in the Gibbs sampler. 
# lower		  = Lower bound for the "L-BFGS-B" method; see 'optim'.
# upper		  = Upper bound for the "L-BFGS-B" method; see 'optim'.
# control 	  = A list of control parameters; see 'optim'.
# interval	  = Interval to pass on to 'optimize' which 
#		    + estimates the reparametrized-version of the dilution
# 		    + parameter; that is, it estimates
# 		    + gamma, where lambda = exp(gamma), -Inf < gamma < Inf. 
# lower.con       = Lower bound for the constant Se model when
#		    + implementing the "L-BFGS-B" method; see 'optim'.
# upper.con       = Upper bound for the constant Se model when
#		    + implementing the "L-BFGS-B" method; see 'optim'.
# control.con     = A list of control parameters for the constant
#		    + Se model; see 'optim'.
# est.con	  = Calculates results using the constant Se model
# 		    + in Zhang, Bilder, and Tebbs (2013), if TRUE.
#
#
#
# Value: 	  (A list of two lists. The first one consists of 
#		   + results from the proposed model, and the second list
# 		   + consists of the results from the constat Se model.)
# par 		  = Maximum likelihood estimates. 
# hessian	  = Estimated hessian matrix. 
# p.value 	  = p-value of the likelihood ratio test to detect dilution. 
#
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
# (4) While the code works for any submodel 'h', which is specified 
#     + below, it only works for the logit regression link function 'g'.
#
# (5) The optional arguments, such as, lower, upper, control, etc., 
#     + are identical to the arguments in optim; i.e., these arguments 
#     + can be specified similarly as in optim.
#
#
Dil.EM.Dorfman <- function(param0,X,Z,Y,cvec,Se,Sp,method="Nelder-Mead"
			   ,hessian=FALSE,emtol=10^(-6),emmaxit=500,GI=12000,
			    burn=200,lower=-Inf,upper=Inf,control=list(), 
			    interval=c(-10000,10), lower.con=-Inf,upper.con=Inf,
			    control.con=list(), est.con=FALSE){
							
#       Define some global variables.
#								
	N <- length(Y)
	J <- length(cvec)
	id <- cumsum(c(0,cvec))
	plen <- length(param0)
	cmax <- max(cvec)
	dil.dt <- data.strc(Z=Z,Y=Y,cvec=cvec,N=N,J=J,id=id)
	dil.dt[ ,1] <- Y

#	E-steps.
#	
	eyf <- E.beta(beta=param0[-plen],X=X,Z=Z,Y=Y,cvec=cvec,gamma=param0[plen],
				   Se=Se,Sp=Sp,J=J,id=id)
	ezf <- E.gamma(gamma=param0[plen],X=X,Z=Z,Y=Y,cvec=cvec,beta=param0[-plen]
				   ,Se=Se,Sp=Sp,J=J,id=id)


#	M-steps: Estimate the regression and dilution parameter separately
#	
	b <- optim(par=param0[-plen],fn=Q.beta,hessian=FALSE,eyij=eyf)$par
	a <- optimize(f=Q.gamma,interval=interval,Z=Z,expt=ezf,cvec=cvec,
					Se=Se,J=J)$minimum
	param1 <- c(b,a)
	
	s <- 1
	while(max(abs(param0-param1)) > emtol){
		param0 <- param1

#	E-steps.
#		
		eyf <- E.beta(beta=param0[-plen],X=X,Z=Z,Y=Y,cvec=cvec,gamma=param0[plen]
						,Se=Se,Sp=Sp,J=J,id=id)
		ezf <- E.gamma(gamma=param0[plen],X=X,Z=Z,Y=Y,cvec=cvec,beta=param0[-plen]
						,Se=Se,Sp=Sp,J=J,id=id)

#	M-steps: Estimate the regression and dilution parameter separately
#
		b <- optim(par=param0[-plen],fn=Q.beta,hessian=FALSE,eyij=eyf)$par
		a <- optimize(f=Q.gamma,interval=interval,Z=Z,expt=ezf,cvec=cvec,
						 Se=Se,J=J)$minimum
		param1 <- c(b,a)
		if(s >= emmaxit) break
		s <- s+1
		print(c(s,param1))
	}

#	Calculate the hessian matrix using Louis's method.
#	
	if(hessian){
		dQ <- cov.dQ(param=param1,Z=Z,Y=Y,X=X,cvec=cvec,Se=Se,Sp=Sp,J=J,id=id,plen=plen)
		dL <- cov.dL(param=param1,dataf=dil.dt,Z=Z,X=X,Se=Se,Sp=Sp,cvec=cvec,
					  N=N,J=J,plen=plen,GI=GI,burn=burn)
		hess <- -dQ-dL
	}else{
		hess <- NULL 
	}

#	Value of the log-likelihood function for the dilution model.
#	
	ll.dil <- (-1)*Dorfman.Dil.loglik(theta=param1,Z=Z,Y=Y,X=X,cvec=cvec,Sp=Sp,
 					  Se=Se,J=J,id=id,plen=plen)

#	Fit the constant Se model in Zhang, Bilder, and Tebbs (2013)
#       when est.con is TRUE.
#
	res.con <- optim(par=param0[-plen],fn=Dorfman.const.loglik,method=method
			 ,hessian=hessian,lower=lower.con,upper=upper.con,control=control.con,
			 Z=Z,Y=Y,X=X,cvec=cvec,Sp=Sp,Se=Se,J=J,id=id)

#	Value of the log-likelihood for the constant Se model.
#			
	ll.con <- (-1)*res.con$value

#	Test statistic and p-value of the likelihood ratio test
#	to detect dilution.
#
	T.LR <- 2*ll.dil - 2*ll.con
	p.val <- 0.5*pchisq(T.LR,0,lower.tail=FALSE)+0.5*pchisq(T.LR,1,lower.tail=FALSE)

	res.dil <- list("par"=param1,"hessian"=hess,"p.value"=p.val)
	res.con2 <- list("par"=res.con$par,"hessian"=res.con$hessian)
	
#	Return all results in a list.
#	
	return(list(dilution=res.dil,constant=res.con2))
}


##########################################################################
#                                                                        #
#    For illustration, simulate group testing data and fit the models.   #
#                                                                        #
##########################################################################


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
# Z	 	  = J by 1 vector of master pooled responses.
# Y 		  = N by 1 vector of individual retest results;
#	            + an individual's test response is 0 if not retested.
#
GT.data.Dorfman <- function(theta,cvec,Se,Sp,X){
	beta <- theta[-length(theta)]
	gamma <- theta[length(theta)]
	J <- length(cvec)
	N <- sum(cvec)
	id <- cumsum(c(0,cvec))
	p <- g(X%*%beta)		# Individual disease probability.
	Ytil <- rbinom(N,1,p)	        # Individual true statuses.
	Y <- rep(0,N)			# Vector of individual retest results.
	Z <- rep(-9,J)			# Vector of pooled testing results.
	for(j in 1:J){
		Ytemp <- Ytil[(id[j]+1):(id[j+1])]
		SY <- sum(Ytemp)+1
		p.tpk <- c( 1-Sp, h(gamma,1:cvec[j],cvec[j],Se) )
		p.pos <- p.tpk[SY]
		Ztemp <- rbinom(1,1,p.pos)
		Z[j] <- Ztemp

		if(Ztemp > 0){
			prob <- ifelse(Ytemp>0,Se,1-Sp)
			Y[(id[j]+1):(id[j+1])] <- rbinom(length(prob),1,prob)
		}

	}
	return(list("pool.resp"=Z,"indv.resp"=Y))
}

#-------------------- Regression model g(.) -----------------------#
#
# Specify the regression link function 'g'. Note, the estimated 
# covariance matrix is correct when g is the logistic link function 
# as specified below.
#
g <- function(t){
	res <- 1/(1+exp(-t))
	return(res)	
}

#------------------ Dilution submodel h(.,.,.) --------------------#
#
# Specify the dilution submodel 'h' and its 1st and 2nd derivatives 
# denoted by dh and d2h. Note, the code works for any submodel h
# as specified below.
# 
# Note: Exponential transformation has been implemented to avoid 
# constraint optimization of the dilution parameter lambda. 
# i.e., we estimate gamma, where lambda = exp(gamma), and 
# -Inf < gamma < Inf. 
#
#
#
# This submodel is derived from the logistic CDF. It has been 
# used for the simulation study; see the paper. 
#
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


#*******************************************************************#
# Specify the true parameter values for simulation
#
lambda.t <- 3.8		     	  	# Dilution parameter lambda.
gamma.t <- log(lambda.t) 	  	# Transformed version of lambda.
beta.t <- c(-3,2,1)		 	# Regression parameter.	
theta.t <- c(beta.t,gamma.t)  	        # Parameters to estimate.
Se.t <- 0.99	 			# Individual testing sensitivity.
Sp.t <- 0.99		 		# Individual testing specificity.

# Specify a vector of pool sizes. 
#
c.vec <- c(rep(5,1000))		        # Configuration 1.
# c.vec <- rep(10,500)		        # Configuration 2.
# c.vec <- c(rep(5,334),rep(10,333))    # Configuration 3.

N <- sum(c.vec)                         # Sample size.

# Set a seed for reproducible results.
#
set.seed(123)
x1 <- rnorm(N,0,.75)		        # x1 from Normal(mu=0,sigma=0.75).	
x2 <- rbinom(N,1,0.1)		        # x2 from Bernoulli(p = 0.1).
X <- cbind(1,x1,x2)			# Design matrix.


# Set R's working directory.
#
setwd(dir = "C:\\programs")

# Simulate group testing data.
#
gt <- GT.data.Dorfman(theta.t,c.vec,Se.t,Sp.t,X)
Z <- gt$pool.resp
Y <- gt$indv.resp

#
# Note, Y must be a N by 1 vector of 0's and 1's. 
# If individuals of a pool are not retested, the corresponding
# individual test responses in Y are 0's.  
#


# Fit the model with the starting value: param0 = (beta0,gamma0).
# One may find beta0 using the constant Se model as shown below. 
# A reasonable choice for gamma0 is 0. 
#
#
J <- length(c.vec)                     # Number of pools.
id <- cumsum(c(0,c.vec))               # Cumulative pool sizes; see, head(id).
beta0 <- optim(beta.t,Dorfman.const.loglik,Z=Z,Y=Y,X=X,cvec=c.vec,Sp=Sp.t,
		Se=Se.t,J=J,id=id)$par
param0 <- c(beta0,0)

res <- Dil.EM.Dorfman(param0=param0,X=X,Z=Z,Y=Y,cvec=c.vec,Se=Se.t,
		      Sp=Sp.t,method="Nelder-Mead",hessian=TRUE,emtol=10^(-6),
		      emmaxit=500,GI=12000,burn=200,interval=c(-10000,10),
		      est.con=TRUE)

#
# Presented are estimation results from the data generated using the seed 123.
# Choose a different seed for different group testing data and results.	
#

# Results from the dilution model.
#
res[[1]]

# Results from the constant Se model with Dorfman testing
# in Zhang, Bilder, and Tebbs (2013)
#
res[[2]]

