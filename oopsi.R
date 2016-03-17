# Functions necessary to run the fast-oopsi algorithm in R
# Based on the Python implementation of the fast-oopsi algorithm
# by liubenyuan (https://github.com/liubenyuan/py-oopsi),
# which is a port from the matlab code by jovo (https://github.com/jovo/oopsi),
# the original developer of the algorithm (http://pqdtopen.proquest.com/doc/375440972.html?FMT=ABS)

# main function: fast(Ff = fluorescence time-series, 
#                     dt = time-step size, 
#                     iter_max = maximum # of iterations, 
#                     update = re-estimate parameters with every iteration)

# louisahsmith March 16, 2016

library(signal)
library(RSEIS)
library(Matrix)

oopsi_mad <- function(Ff){
  
	median(abs(Ff - median(Ff)))
  
}

oopsi_m <- function(gamma, Tt){
  
	ddata <- diag(-gamma, Tt-1)
	ddata <- cbind(rbind(rep(0, Tt-1), ddata), rep(0, Tt))
	diag(ddata) <- 1
	M <- Matrix::Matrix(ddata, sparse = T)
	
	return(M = M)
	
}

oopsi_init_par <- function(Ff, dt){
  
	epsilon <- 1e-16
	Tt <- length(Ff)
	Ff <- RSEIS::detrend(Ff)
	Ff <- (Ff - min(Ff))/(max(Ff) - min(Ff)) + epsilon
	a <- 1
	b <- median(Ff)
	lam <- 1
	gam <- 1 - dt/1
	sig <- oopsi_mad(Ff) * 1.4826
	P <- list(Tt = Tt, dt = dt, gamma = gam, alpha = a, 
		beta = b, sigma = sig, lambda = lam)
	
	return(list(Ff = Ff, P = P))
	
}

oopsi_est_map <- function(Ff, P){
  
	Tt <- P$Tt
	dt <- P$dt
	gam <- P$gamma
	a <- P$alpha
	b <- P$beta
	sig <- P$sigma
	lam <- P$lambda
	n <- 0.01 + rep(0, Tt)
	C <- signal::filter(filt = 1, a = c(1, -gam), x = n)
	llam <- (lam*dt) * rep(1, Tt)
	M <- oopsi_m(gam, Tt)
	grad_lnprior <- t(M)%*%llam
	H1 <-as(diag(1, Tt), "sparseMatrix") *  c(((a^2)/(sig^2)))
	z <- 1
	
	while(z > 1e-13){
	  
		D <- Ff - a*C - b
		lik <- 1/(2*(sig^2)) * (t(D)%*%D)
		post <- lik + t(llam)%*%n - z*sum(log(n), na.rm = T)
		s <- 1
		d <- 1
		
		while((sqrt(sum(d^2)) > 5e-2) & (s > 1e-3)){
		  
			glik <- -a/(sig^2)*(Ff - a*C - b)
			g <- c(glik + c(grad_lnprior@x)) - z*(t(M)%*%(1/n))
			H2 <- diag(1/(n^2), Tt)
			H <- H1 + z*(t(M)%*%H2%*%M)
			d <- solve(H, g)
			hit <- n/(M%*%d)
			hit <- hit[hit > 0]
			s <- ifelse(any(hit < 1), 0.99*min(hit), 1)
			post1 <- post + 1
			
			while(post1 > (post + 1e-7)){
			  
				C1 <- c(C) - s*c(d@x)
				n <- c((M%*%C1)@x)
				D <- Ff - a*C1 - b
				lik1 <- 1/(2*(sig^2))*(t(D)%*%D)
				post1 <- lik1 + t(llam)%*%n - z*sum(log(n), na.rm = T)
				s <- s/5
				
				if(s < 1e-20){
				  
					break
				  
				}
			}
			
			C <- C1
			post <- post1
			
		}
		
		z <- z/10
		
	}
	
	n[1:2] <- 1e-8
	n <- n/max(n)
	return(list(n = n, C = C, post = post))
	
}

oopsi_est_par <- function(n, C, Ff, P){
  
	Tt <- P$Tt
	dt <- P$dt
	gam <- P$gamma
	sig <- P$sigma
	lam <- P$lambda
	a <- 1
	b <- sum(Ff-C)/Tt
	D <- Ff - a*C - b
	mse <- t(D)%*%D
	sig <- sqrt(mse/Tt)
	lam <- Tt/(dt*sum(n))
	
	return(list(Tt = Tt, gamma = gam, alpha = a,
		beta = b, sigma = sig, lambda = lam, dt = dt))
	
}

fast <- function(Ff, dt = 0.02, iter_max = 1, update = T){

    listFP <- oopsi_init_par(Ff, dt)
    P <- listFP$P
    Ff <- listFP$Ff
    listnCpost <- oopsi_est_map(Ff, P)
    n <- listnCpost$n
    C <- listnCpost$C
    post <- listnCpost$post
    post_max <- post
    n_best <- n
    C_best <- C
    ml <- rep(1, iter_max)
    ml[1] <- post
    
    for(i in 1:iter_max){
      
    	if(update == T){
    	  
    		P <- oopsi_est_par(n, C, Ff, P)
    	}
      
    	listnCpost <- oopsi_est_map(Ff, P)
    	n <- listnCpost$n
    	C <- listnCpost$C
    	post <- listnCpost$post
    	
    	if(post > post_max){
    	  
    		n_best <- n
    		C_best <- C
    		post_max <- post
    	}
    	
    	ml[i+1] <- post
    	
    	if(((abs(ml[i+1] - ml[i])/ml[i]) < 1e-3) |
    		any(abs(ml[1:(i+1)] - ml[i+1]) < 1e-5)){
    	  
    		break
    	  
    		}
    }
    
    return(list(d = n_best, Cz = C_best))
}
