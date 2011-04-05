lasso <- function(y, ...) 
	UseMethod('lasso')
        
getLambdaNcoef <- function(y, x, lambda1, nCoef, track=FALSE, model='linear', standardize=FALSE){

  options("show.error.messages"=FALSE)
  fit <- try(penalized(y, x, lambda1=lambda1, standardize=standardize, steps=20, trace=track, model=model), silent=TRUE)
  options("show.error.messages"=TRUE)

  i <- 1
  while(i<=5 & class(fit)=='try-error'){
    lambda1 <- lambda1*2
    options("show.error.messages"=FALSE)
    fit <- try(penalized(y, x, lambda1=lambda1, standardize=standardize, steps=20, trace=track, model=model), silent=TRUE)
    options("show.error.messages"=TRUE)
  }
  
  if (class(fit)!='try-error'){
    nCoefs <- sapply(fit, function(x) length(coef(x)-1))
    lambdas <- sapply(fit, function(x) x@lambda1)
    lambda <- lambdas[which.min(abs(nCoefs - nCoef))]
  }else{
    cat('Selecting appropriate lambda for', floor(nCoef), 'non-zero coefficients failed. \n')
    lambda <- NULL
  }
  lambda
}

lasso.stability <- function(y, x=NULL, alpha=.5, subsampling=.5, nSubsampling=200, model='linear', pi_th=.6,  alpha.fwer=1, lambda1=NULL, steps=10, track=FALSE, standardize=FALSE, ...){
	
	require(penalized)
	if (is.null(x)) {
		data <- y
	}else {
		data <- NULL
		data$y <- y
		data$x <- x
	}
	
	px	<- ncol(data$x)
	nSample <- nrow(data$x)
	lambda <- ifelse(is.null(lambda1), sqrt(nSample + px)/10, lambda1)
	
	if (alpha==0){ 
		alpha <- runif(nSubsampling, .5, 1.5)
	}else{ 
		alpha <- rep(alpha, nSubsampling)
	}
	
	q_lambda <- sqrt(0.8*alpha.fwer*px)
	
	subset.sample <- sample(1:nSample, nSample * subsampling)
	subset.var <- sample(1:px, px*.5)
	
	x1 <- data$x
	x1[, subset.var] <- x1[, subset.var] * alpha[1]

        lambda <- getLambdaNcoef(data$y[subset.sample], x1[subset.sample,], lambda1=lambda, nCoef=q_lambda, model=model)
	if(is.null(lambda))
          stop('Estimation of lambda giving', q_lambda, 'non-zero coefficients failed.')
	
	mat	<- NULL
	n	<- NULL
	iter 	<- 1
	while ((iter <= nSubsampling)){
		subset.sample	<- sample(1:nSample, nSample * subsampling)
		subset.var	<- sample(1:px, px*.5)
		x1 <- data$x
		x1[, subset.var] <- x1[, subset.var] * alpha[1]
                options("show.error.messages"=FALSE)
		fit <- try(penalized(data$y[subset.sample], x1[subset.sample,], lambda1=lambda, standardize=standardize, steps=1, trace=track, model=model), silent=TRUE)
                options("show.error.messages"=TRUE)
		if (class(fit)!='try-error') {
			mat <- try(cbind(mat, fit@penalized))
			if (class(mat)!='try-error')
                          iter <- iter + 1
		}
	}

	beta.freq <- rowSums(mat!=0)/ncol(mat)
	beta <- (beta.freq >= pi_th)*1
	
	res <- list(beta.freq=beta.freq, beta=beta, mat=mat, residuals=data$y)
        class(res) <- 'lol'
        res
} 
	

lasso.multiSplit <- function(y, x=NULL, lambda1=NULL,  nSubsampling=200, model='linear', alpha=0.05, gamma.min=0.05, gamma.max=0.95, track=FALSE, ...){

  if (is.null(x)) {
    x <- y$x
    y <- y$y
  }

  px	<- ncol(x)
  nSample <- nrow(x)

  if (nSample!=length(y))
    stop('The number of samples must match between input x and y')

  if(is.null(lambda1))
          lambda1 <- sqrt(px+nSample)/10
	
	pmat <- matrix(1, ncol=nSubsampling, nrow=px)
	pmat.adj <- matrix(1, ncol=nSubsampling, nrow=px)
	i <- 0
  subsampling <- 0.5
        lambda <- getLambdaNcoef(y[1:(nSample/2)], x[1:(nSample/2),], lambda1=lambda1, nCoef=sqrt(.8*px))
                
 	while (i<nSubsampling){
		subset <- sample(1:nSample, nSample * subsampling)
                options("show.error.messages"=FALSE)
                fit <- try(penalized(y[subset], x[subset,], lambda1=lambda, steps=1, trace=track, model=model, ...), silent=TRUE)
                options("show.error.messages"=TRUE)

		if (class(fit)!='try-error') {
                  i <- i+1
                  coefs <- fit@penalized
                  nCoef <- sum(coefs!=0)
                  if ( nCoef > 0){
                    fit2 <- lm(y[-c(subset)]~x[-c(subset), coefs!=0])
                    if(!any(is.na(fit2$coefficients)))
                      pmat.adj[coefs!=0, i] <- summary(fit2)$coefficients[-1, 4] * nCoef
		}}
	}

	pmat.adj[pmat.adj>1] <- 1
	
	Q_gamma <- NULL
	for (gamma in seq(gamma.min, gamma.max, len=10))
		Q_gamma <- cbind(Q_gamma, apply(pmat.adj, 1, function(x) quantile(x/gamma, prob=gamma)))

	P <- (1 - log(0.05)) * apply(Q_gamma, 1, min)
	P[P>=alpha] <- 0
	beta <- rep(0, px)
	beta[P!=0] <- log(1/(P[P!=0]))
		
  res <- list(beta=beta, mat=Q_gamma, residuals=y, pmat=pmat.adj)
  class(res) <- 'lol'
  res
}


	
lasso.simultaneous <- function(y, x=NULL, model='linear', nSubsampling=200, alpha=.5, lambda1=NULL, track=FALSE, ...){

	if (is.null(x)) {
		x <- y$x
		y <- y$y
	}

	nSample <- nrow(x)
	px <- ncol(x)
	mat <- matrix(0, nrow=nSubsampling, ncol=px)
	n <- 0
	q_lambda <- sqrt(0.8*px)
	lambda <- ifelse(is.null(lambda1), sqrt(nSample + px)/10, lambda1)

        for (i in 1:nSubsampling){
		subset.var <- sample(1:px, px * .5 )
		sample.subset <- sample(1:nSample, nSample/2)
		x1 <- x
		x1[, subset.var] <- x1[, subset.var] * alpha

		data <- NULL
		data$x <- x1[sample.subset, ]
		data$y <- y[sample.subset]
                
                if(i==1)
                  lambda <- getLambdaNcoef(data$y, data$x, lambda1=lambda, nCoef=q_lambda)
                
                options("show.error.messages"=FALSE)
		B1 <- try(penalized(data$y, data$x, lambda1=lambda, steps=1, trace=track, model=model, standardize=FALSE, ...), silent=TRUE) 
                options("show.error.messages"=TRUE)
	
		data <- NULL
		data$x <- x1[-sample.subset, ]
		data$y <- y[-sample.subset]
                options("show.error.messages"=FALSE)
		B2 <- try(penalized(data$y, data$x, lambda1=lambda, steps=1, trace=track, model=model, standardize=FALSE,...), silent=TRUE) 
                options("show.error.messages"=TRUE)
		if( class(B1)!='try-error'  &  class(B2)!='try-error' ){
			mat[i, ] <- B1@penalized!=0 & B2@penalized!=0
			n <- n + 1
		}
	}
	beta <- colSums(mat!=0)/n
	
	res <- list(beta=beta, n=n, mat=mat)
        class(res) <- 'lol'
        res
}

lmMatrixFit <- function(y, x=NULL, mat, th=NULL){

	require(Matrix)

        if (is.null(x)) {
          x <- y$x
          y <- y$y
	}
	x <- scale(x)
        y <- scale(y)

	if(nrow(mat)!=ncol(y) | ncol(mat)!=ncol(x))
		stop('Input "mat" should be a ncol(y)-by-ncol(x) matrix!')
        
	nSample <- nrow(x)
	px <- ncol(x)
	py <- ncol(y)
	coefMat <- Matrix(0, nrow=py, ncol=px)
	resMat <- y	
	pvalMat <- Matrix(1, nrow=py, ncol=px)

        if((!is.null(th))&(!is.logical(mat))) 
           mat[abs(mat)<=th] <- 0
        nReg <- rowSums(mat!=0)

        if (is.logical(mat)){
          mat <- mat*1
        }else{
          for (i in which(nReg>=nSample)){
            b <- mat[i,]
            b[sort.list(abs(b), decreasing=TRUE)[(nSample):length(b)]] <- 0		
          }}
        

        for (i in which(nReg!=0)){
		b <- mat[i, ]
		res <- lm(y[,i] ~ x[, b!=0])
		pval <- summary(res)$coefficients[-1,4]
		b[b!=0]  <- res$coefficient[-1]
		resMat[, i] <- res$residuals
		coefMat[i, ] <- b
		pvalMat[i, b!=0] <- pval
	}
	
	res <- list(coefMat=coefMat, resMat=resMat, pvalMat=pvalMat)
        class(res) <- 'lolMatrix'
        res
}



matrixLasso <- function(y, x=NULL, method='cv', nameControl=FALSE, standardize=FALSE, track=0, lambda1=NULL, nFold=10, ...){

  require(Matrix)
   	cl <- match.call()
	allMethods <- c('cv', 'stability', 'multiSplit', 'simultaneous')
	if (is.numeric(method))
		method <- allMethods[method]
	
	if (is.null(x)) {
		data <- y
	}else{
		data <- NULL
		data$y <- y
		data$x <- x
	}	

  if(nrow(data$y)!=nrow(data$x))
    stop('The number of rows/samples should be the same for both data matrices!')
    
	class(data) <- method

	if (standardize) {
		data$y <- scale(data$y); 
		data$x <- scale(data$x)	
	}
	p <- dim(data$y)[2] 
	n <- dim(data$y)[1]
	q <- dim(data$x)[2]

	if (length(data$y) < 50 & dim(data$x)[2] < 50){ 
	
		if(identical(data$x[,1], data$y[,1])){
			X <- Matrix(0, p*n, p*(q-1))
			for (i in 1:p)
			X[((i-1)*n+1):(i*n), ((i-1)*(q-1)+1):(i*(q-1))] <- data$x[, -i]
		
		}else{
			X <- Matrix(0, p*n, p*q)
			for (i in 1:p)
			X[((i-1)*n+1):(i*n), ((i-1)*q+1):(i*q)] <- data$x
		}

		Y <- as.vector(data$y)
		Y <- Matrix(Y, length(Y), 1)
	
		data$y <- Y
		data$x <- X
		p <- dim(data$y)[2] 
		n <- dim(data$y)[1]
		q <- dim(data$x)[2]
	}else{
		X <- data$x
		Y <- data$y
	}

	resMat <- Matrix(0, n, p)
	coefMat<- Matrix(0, p, q)
	for (i in seq(p)){
		newX <- X
		if ((!is.null(colnames(data$y))) & nameControl) 
			if (colnames(data$y)[i] %in% colnames(data$x)) 
				newX <- data$x[, colnames(data$x) != colnames(data$y)[i]]
		
		data.v <- NULL
		data.v$y <- data$y[,i]
		if(class(data$y) == 'data.frame'){
			data.v$y <- unlist(data$y[,i])
		}
		data.v$x <- newX
		class(data.v) <- method
		fit <- try(lasso(data.v, standardize=standardize, track=ifelse(track==2, 1, 0), 
			lambda1=lambda1, nFold=nFold, ...), silent=TRUE)
		if (!('try-error'%in%class(fit))){	
			if ((!is.null(colnames(data$y))) & nameControl) {
				if (colnames(data$y)[i] %in% colnames(data$x)) {
				coefMat[colnames(data$x)!=(colnames(data$y)[i]), i] <- fit$beta
			}}else{
				coefMat[i,] <- fit$beta
				if (track>0)
					print(paste('Response', i, 'has', sum(fit$beta!=0), 'regulators.'))
			}
			resMat[,i] <- fit$residuals
		}else {
			resMat[,i] <- data.v$y 
		}
	}	
	
	colnames(coefMat) <- colnames(x)
	rownames(coefMat) <- colnames(y)
	
	res <- list(coefMat=coefMat, fit=fit, resMat=resMat)
        class(res) <- 'lolMatrix'
        res
}



lasso.cv <- function(y, x=NULL, lambda1=NULL, model='linear', steps=15, minsteps=5, log=TRUE, track=FALSE, standardize=FALSE, unpenalized=~0,  nFold=10, nMaxiter = Inf, ...){
	
	require(penalized)
	
	if (is.null(x)) {
		data <- y
	}else {
		data <- NULL
		data$y <- y
		data$x <- x
	}
	
	n <- length(data$y)
	q <- dim(data$x)[2]
	l.upper <- q/10
	l.lower <- 0
	conv <- FALSE
	iter <- 1
	time1 <- c(0,0,0)
	
	if(is.null(lambda1))
			lambda1 <- mean(c(l.upper, l.lower))

	while (iter<=5 & (!conv) & time1[3]<40){
			
		time1 <- proc.time()
                options("show.error.messages"=FALSE)
		fit <- try(profL1 (data$y, penalized=data$x, minlambda1=lambda1, fold=nFold, steps=steps, 
			minsteps=minsteps, log=log, save.predictions=FALSE, trace=track, 
			standardize=standardize, model=model, maxiter=nMaxiter), silent=TRUE) 
                options("show.error.messages"=TRUE)
		time2 <- proc.time()
		time1 <- time2-time1
		
		if (class(fit)=='try-error' ){
			
				l.lower <- lambda1 
				if (track) print('Estimated lambda too small, increase lambda\n')
			
		}else{ if (which.max(fit$cvl) == length(fit$cvl) ){
					
				l.upper <- lambda1
				if (track) print('Not converged, repeat optimization with smaller lambda\n')
				
		}else{ if( fit$lambda[1] < lambda1 ){ # initial lambda larger than maxlambda1,					
				l.upper <- fit$lambda[1]
				if (track) print('Input min lambda larger than max lambda, repeat\n')
			
		}else{ 
				conv <- TRUE 	
		}}}
				
		iter <- iter+1
		lambda1 <- mean(c(l.upper, l.lower))
	}
	
	residuals <- data$y
	beta <- rep(0, q)
	lambda <- 0
	if ((class(fit)!='try-error') & (length(fit)>1)){
		l <- fit$lambda[max(c(1, which.max(fit$cvl)))]
		fit <- try(optL1 (data$y, data$x, minlambda1=max(0.0001, l-2), maxlambda1=l+2, fold=10,  trace=track, standardize=standardize, model=model, maxiter=nMaxiter), silent=TRUE) 
		lambda <- ifelse (class(fit)!='try-error', fit$lambda, l)
		fit <- try(penalized(data$y, data$x, lambda1=lambda, trace=track, standardize=standardize, model=model, maxiter=nMaxiter), silent=TRUE)
		if (class(fit)!='try-error'){
			residuals <- fit@residuals
			beta <- fit@penalized
			lambda <- fit@lambda1
	}}
	
	res <- list(fit=fit, residuals=residuals, beta=beta, lambda=lambda, conv=conv)
        class(res) <- 'lol'
        res
}

plotGW <- function(data, pos, marks=NULL, fileType='png', file='plotGW', width=1000, height=500,  autoscale=FALSE, col=c('lightblue', 'lightgreen', 'darkblue', 'darkgreen'), legend=1:10, ylab='', pch=19, cex.axis=1.2 ,cex.lab=1.2, cex=.5, legend.pos='bottomright', mtext=NULL, mtext.side=2, mtext.at=NULL, mtext.line=3, ...){

  if (is.vector(data))
          data <- cbind(data, data)

  if(class(pos)=='matrix'){
	mark1 <- as.matrix(pos[,'chromosome_name']);
        rownames(mark1) <- 1:nrow(pos)
        mark1 <- unique(mark1, drop=F)
      }
  if(class(pos)=='numeric'){
	mark1 <- as.matrix(pos)
        rownames(mark1) <- 1:length(pos)
        mark1 <- unique(mark1, drop=F)
      }
  if(class(pos)=='data.frame'){
	mark1 <- as.matrix(pos[,'chromosome_name']);
        rownames(mark1) <- 1:nrow(pos)
        mark1 <- unique(mark1, drop=F)
      }
  
  if (autoscale & min(data) < 0 & max(data) > 0) {
	s <- - max(data)/min(data)
	data[data<0] <- data[data<0] * s
      }
  if (fileType=='pdf'){
    if (width>20 | height>20){
      width <- 10
      height <- 8
    }
    pdf(paste(file, fileType, sep='.'), width=width, height=height)
  }else{
    png(paste(file, fileType, sep='.'), width=width, height=height)
  }
  
  plot(data[,1], xaxt='n',  ylab=ylab, col=col[1], ylim=c(min(data), max(data)), xlab='Chromosome', pch=pch, cex=cex,...)
  for (i in 1:ncol(data))
    points(data[,i],  col=col[i], pch=pch, cex=cex, ...)
    
  axis(1, at=as.numeric(rownames(mark1)), labels=as.vector(mark1), las=1, cex.axis=cex.axis)
  abline(v=as.numeric(rownames(mark1)), lty=2, col='grey')
  points(rep(0, nrow(data)), col='black', pch=pch, cex=cex, ...)
  if (!is.null(marks))
    points(x=marks, y=rep(0, length(marks)), pch=4, cex=2.5, col=2)
  legend(legend.pos, legend=legend, pch=pch, inset=.01, col=col, cex=cex*2)
  if (is.null(mtext.at))
    mtext.at <- c(min(data)/2, max(data)/2)
  mtext(mtext, side=mtext.side, at=mtext.at, line=mtext.line, cex=cex.axis)
  dev.off()
}


print.lol <- function(x,...){
  cat(" Non-zero coefficients:",  sum(x$beta!=0))
  cat("\n from a total of", length(x$beta), "predictors\n")  
}

print.lolMatrix <- function(x,...){
  cat(" Non-zero coefficients in total:", sum(x$coefMat!=0))
  cat("\n from a total of", ncol(x$coefMat), "predictors")
  cat("\n and", nrow(x$coefMat), "responses\n")
}
