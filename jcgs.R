computeCX <- function(X, mode = "Bernoulli"){
         P <- X%*%t(X)
         if(mode == "Bernoulli"){
         D <- rowSums(P - P^2)
         }
         if(mode == "Poisson"){
           D <- rowSums(P)
         }
         D[D < 0] <- 0
         tmp1 <- eigen(t(X)%*%X)
         tmp2 <- tmp1$vectors %*% diag(1/(tmp1$values^2)) %*% t(tmp1$vectors)
         CX <- sum(diag(tmp2 %*% t(X*sqrt(D)) %*% (X*sqrt(D))))
         return(sqrt(CX))
}

stfp <- function(A, dim){
    require("irlba")

    A.svd <- irlba(A, nu = dim, nv = dim)
    A.svd.values <- A.svd$d[1:dim]
    A.svd.vectors <- A.svd$v[,1:dim]
    if(dim == 1)
        A.coords <- sqrt(A.svd.values) * A.svd.vectors
      else
        A.coords <- A.svd.vectors %*% diag(sqrt(A.svd.values))
    
   return(A.coords)
}

rg.sample <- function(P,mode = "Bernoulli"){
  n <-  nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  if(mode == "Bernoulli"){
  A[col(A) > row(A)] <- runif(n*(n-1)/2)
  A <- (A + t(A))
  A <- (A < P) + 0 ;
  diag(A) <- 0
  }
  if(mode == "Poisson"){
    lambdas <- P[col(P) > row(P)]
    lambdas <- lambdas*(lambdas > 0)
    A[col(A) > row(A)] <- rpois(n*(n-1)/2, lambdas)
    A <- (A + t(A))
  }
  return(A)
}

procrustes <- function(X,Y, type = "I"){
    if(type == "C"){
        X <- X/norm(X, type = "F")
        Y <- Y/norm(Y, type = "F")
    }
    if(type == "D"){
        tX <- rowSums(X^2)
        tX[tX <= 1e-15] <- 1
        tY <- rowSums(Y^2)
        tY[tY <= 1e-15] <- 1
        X <- X/sqrt(tX)
        Y <- Y/sqrt(tY)
    }

    tmp <- t(X) %*% Y
    tmp.svd <- svd(tmp)
    W <- tmp.svd$u %*% t(tmp.svd$v)
    return(list(error = norm(X%*%W - Y, type = "F"), W = W))
}

## Compute the orthogonal procrustes statistics given the embedding
T.semipar <- function(X,Y, type = "I"){
    if(type == "I"){
        return(procrustes(X,Y, type)$error/(computeCX(X) + computeCX(Y)))
    }
    else if(type == "C"){
        norm.X <- norm(X, type = "F")
        norm.Y <- norm(Y, type = "F")
        return(procrustes(X, Y, type = "C")$error/(2*(computeCX(X)/norm(X, type = "F") +
                                                          computeCX(Y)/norm(Y, type = "F"))))
    }
    else {
        return(procrustes(X, Y, type = "D")$error/(2*(1/min(sqrt(rowSums(X^2))) +
                                                    1/min(sqrt(rowSums(Y^2))))))
    }
}

## Compute the orthogonal procrustes statistics along with the
## critical values and the p-values given the embedding.
## X and Y are the embeddings
## B is the number of bootstrap samples
## alpha is the significance level of the null hypothesis
## type is the type of test;
## type = "I" corresponds to test of equality;
## type = "C" corresponds to test of equality up to scaling
## type = "D" corresponds to test of equality up to diagonal transformation
T.semipar.bootstrap <- function(X,Y,B,alpha, type, mode = "Bernoulli"){
    TX.vec <- numeric(B)
    TY.vec <- numeric(B)
    stat <- procrustes(X,Y, type)$error

    for(b in 1:B){
            A <- rg.sample(X %*% t(X), mode)
            Xhat0 <- stfp(A,ncol(X))
            A <- rg.sample(X %*% t(X), mode)
            Xhat1 <- stfp(A,ncol(X))

            A <- rg.sample(Y %*% t(Y), mode)
            Yhat0 <- stfp(A,ncol(Y))
            A <- rg.sample(Y %*% t(Y), mode)
            Yhat1 <- stfp(A,ncol(Y))
        
        TX.vec[b] <- procrustes(Xhat0, Xhat1, type)$error
        TY.vec[b] <- procrustes(Yhat0, Yhat1, type)$error
   }
    return(list(cv1 = quantile(TX.vec, probs = c(1 - alpha)), 
                cv2 = quantile(TY.vec, probs = c(1 - alpha)), 
                p.val = max((sum(TX.vec > stat) + 0.5)/length(TX.vec), (sum(TY.vec > stat) + 0.5)/length(TY.vec)), 
                stat = stat))
}

## Run the experiment for section 4.2
## herm_Graph.Rd is a RData file storing the chemical synapes graph Ac
## and the electrical gap junction graph Ag.
c.elegans <- function(){
    load("herm_Graph.Rd")
    Ac <- (Ac + t(Ac))/2
    Ag <- (Ag + t(Ag))/2

    Ac[Ac > 0] <- 1
    Ag[Ag > 0] <- 1

    sum(Ac)
    sum(Ag)

    Xhat.Ac <- stfp(Ac, dim = 6)
    Xhat.Ag <- stfp(Ag, dim = 6)

    result <- T.semipar.bootstrap(Xhat.Ac, Xhat.Ag, B = 1000, alpha = 0.05, type = "C")
    return(result)
}

## This experiment is best done by parallelization
## To obtain Table 4.1 in the manuscript, run this experiment 
## for various choices of epsilon and n.
## To obtain Table 4.2 in the manuscript, run this experiment
## for various choices of epsilon and n. Also set type = "C".
## in the call to T.semipar.bootstrap
bootstrap.experiment1 <- function(epsilon, n, nmc){
  
  library("foreach")
  library("doParallel")
  registerDoParallel(1) ## Replace by number of threads
  library("igraph")
  
  B0 <- matrix(c(0.5,0.2, 0.2, 0.5), nrow = 2)
  rho <- c(0.4,0.6)
  
  B1 <- matrix(c(0.5 + epsilon,0.2,0.2,0.5 + epsilon), nrow = 2)
  B <- 200
  alpha <- 0.05
  
  result <- foreach(i = 1:nmc, .combine = 'cbind') %dopar% {
    tau <- as.vector(rmultinom(1, n, prob = rho))
    A0 <- as_adj(sbm.game(n, B0, tau))
    A1 <- as_adj(sbm.game(n, B0, tau))
    A2 <- as_adj(sbm.game(n, B1, tau))
    
    Xhat0 <- stfp(A0,2)
    Xhat1 <- stfp(A1,2)
    Xhat2 <- stfp(A2,2)
    
    null.i <- T.semipar.bootstrap(Xhat0, Xhat1, B, alpha, type = "I")
    alt.i <-  T.semipar.bootstrap(Xhat0, Xhat2, B, alpha, type = "I")
    resulti <- c(null.i, alt.i)
  }
  
  return(result)
}

## This experiment generates the figure for section 4.3
runKKI <- function(){
  
  ## The data file for glist is to be obtained from the openconnectome project at
  ## http://openconnecto.me/data/public/MR/m2g_v1_1_0/KKI2009/outputs/smallgraphs/mat/
  ## The above link contains download link for the 42 small KKI brain graphs.
  ## Once they are downloaded, save them into the working directory of the current source file.
  
  require("R.matlab")
  glist <- list()
  files.list <- dir(pattern = "\\.mat")
  i <- 1
  for(f in files.list){
    glist[[i]] <- readMat(f)$graph
    i <- i + 1
  }
  
  Xhat.list <- list()
  (m <- length(glist))
  dmax <- 10
  
  for(i in 1:m) {
    gi <- as.matrix(glist[[i]][])
    gi <- (gi + t(gi))/2
    gi.log <- log(gi + 1)
    Xhat.list[[i]] <- stfp(gi.log, dmax)
  }
  
  Smat <- Pmat <- matrix(0, m, m)
  d <- 4
  type <- "C"
 
  require(parallel)
  require(doMC)
  registerDoMC(cores=detectCores()-1)
  getDoParWorkers()
  
   for(i in 1:m) {
      cat("d = ", d, ", i = ", i, "\n")
      out <- foreach (j=1:m) %dopar% {
        #        for(j in 1:m) {
        out <- T.semipar.bootstrap((Xhat.list[[i]])[,1:d], (Xhat.list[[j]])[,1:d], B=200, 
                                   alpha=0.05, type=type, mode = "Poisson")
       
      }
      Pmat[i,] <- sapply(out,"[[",3)
      Smat[i,] <- sapply(out,"[[",4)
    }
    save(Pmat,Smat,file=paste0("Smat-Pmat-type",type,"-dhat",d,".Rbin"))
  
  Matrix::image(Matrix(Smat),useAbs=FALSE)
  Matrix::image(Matrix(Pmat),useAbs=FALSE)
  
  require(popbio)
  pair <- read.csv("kki42_subjectinformation.csv")
  (subj <- pair[,1])
  (rord <- pair[seq(1,41,by=2),1])
  (cord <- pair[seq(2,42,by=2),1])
  subjp <- cbind(rord,cord)
  subjp <- subjp[order(subjp[,1]),]
  ord <- as.vector(t(subjp))
  
  d <- 4
  type <- "C"
  #    par(mfrow=c(2,2))
  #    for (d in 2:5) {
    print(load(paste0("Smat-Pmat-type",type,"-dhat",d,".Rbin")))
    rownames(Pmat) <- colnames(Pmat) <- 1:42
    Pmat0 <- Pmat[ord,ord]
    #    Smat <- Smat[rord,cord]
    image2(Pmat0,text.cex=0.3,label.cex=0.5,round=2,box.offset=0,labels=c(2,3))
    #        title(paste("dhat =",d))
    dev.print(pdf,paste0("pmat-kki70-type",type,"-dhat",d,".pdf"))
  
}

res <- bootstrap.experiment1(epsilon = .05, n=200, nmc=10)

