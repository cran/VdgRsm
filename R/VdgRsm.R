#library(permute)  # permute
#library(AlgDesign)  # use quad(.)
#library(DoE.wrapper)

####################   spv function   #####################

spv <- function(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
                des.names = c("Design 1","Design 2","Design 3"),
                scale = TRUE, add.pts = TRUE, big.RAM = FALSE){  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.grid = 15
  n.var <- ncol(design.matrix)
  
  if(n.var > 7){stop("The maximum number of design factors allowed is 7 \n")}
  if((n.var == 6 || n.var == 7) && big.RAM == FALSE){
    print("Simulated points will be not added unless option big.RAM = TRUE")
    add.pts = FALSE
  }
  
  
  if(scale == TRUE){ 
    temp.radii <- apply(design.matrix, MARGIN = 1, norm2)
    scale.factor <- sqrt(n.var)/max(temp.radii)
    design.matrix <- scale.factor*design.matrix
  }
  
  design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
  N.generated <- 15*n.var*(2^(n.var))
  model.X <- model.matrix( ~quad(.) , design.matrix)
  nrow <- nrow(model.X)
  ncol <- ncol(model.X)
  M.inv<- nrow*solve(t(model.X)%*%model.X)
  
  # Design 2
  if(!is.null(design.matrix.2)){
    n.var.2 <- ncol(design.matrix.2)
    if(scale == TRUE){ 
      temp.radii.2 <- apply(design.matrix.2, MARGIN = 1, norm2)
      scale.factor.2 <- sqrt(n.var.2)/max(temp.radii.2)
      design.matrix.2 <- scale.factor.2*design.matrix.2
    }
    design.matrix.2<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.2))),byrow=FALSE,ncol=n.var.2))  
    model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
    nrow.2 <- nrow(model.X.2)
    ncol.2 <- ncol(model.X.2)
    if(ncol != ncol.2){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
  } # Design 2
  # Design 3
  if(!is.null(design.matrix.3)){
    
    n.var.3 <- ncol(design.matrix.3)
    if(scale == TRUE){ 
      temp.radii.3 <- apply(design.matrix.3, MARGIN = 1, norm2)
      scale.factor.3 <- sqrt(n.var.3)/max(temp.radii.3)
      design.matrix.3 <- scale.factor.3*design.matrix.3
    }
    design.matrix.3<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.3))),byrow=FALSE,ncol=n.var.3))  
    model.X.3 <- model.matrix( ~quad(.) , design.matrix.3)
    nrow.3 <- nrow(model.X.3)
    ncol.3 <- ncol(model.X.3)
    if(ncol != ncol.3){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.3<- nrow.3*solve(t(model.X.3)%*%model.X.3)   
  } # Design 3  
  
  
  Matrix.V <- matrix(numeric(0), ncol = 4, nrow = n.grid)
  colnames(Matrix.V) <- c("Radius","  max.V","  min.V"," average.V")
  
  if(!is.null(design.matrix.2)){
    Matrix.V.2 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V.2) <- c("Radius","  max.V","  min.V"," average.V")  
  }# Design 2
  
  if(!is.null(design.matrix.3)){
    Matrix.V.3 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V.3) <- c("Radius","  max.V","  min.V"," average.V")  
  }# Design 3
  
  
  # Generate points on a sphere
  rand.sphere <- matrix(c(runif(n=(n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))),byrow= TRUE, ncol=n.var)
  radius <- seq(from = 0, to = sqrt(n.var), length.out = n.grid)
  norm2.rand <- apply(rand.sphere, MARGIN = 1, norm2)
  norm2.rand.rep <- matrix(rep(norm2.rand, each = n.var),byrow= TRUE, ncol=n.var)
  
  set.seed(1234567)
  gen.array <- array(rand.sphere/norm2.rand.rep,dim = c(N.generated,n.var,n.grid))
  Array.with.R <- rep(radius,each = (N.generated*n.var))*gen.array
  
  
  for(ii in 1:n.grid){
    pred.model  <- model.matrix( ~quad(.), as.data.frame(Array.with.R[,,ii]))
    Result.V    <- pred.model%*%M.inv%*%t(pred.model)
    Var.pred    <- diag(Result.V)
    minimum <- min(Var.pred)
    maximum <- max(Var.pred)
    average <- mean(Var.pred)
    R.rho <- radius[ii]
    Matrix.V[ii,] <- c(R.rho,maximum,minimum,average)  
    
    # Design 2
    if(!is.null(design.matrix.2)){
      Result.V.2    <- pred.model%*%M.inv.2%*%t( pred.model)
      Var.pred.2    <- diag(Result.V.2)
      minimum.2 <- min(Var.pred.2)
      maximum.2 <- max(Var.pred.2)
      average.2 <- mean(Var.pred.2)
      Matrix.V.2[ii,] <- c(R.rho,maximum.2,minimum.2,average.2)     
    } 
    # Design 3
    if(!is.null(design.matrix.3)){
      Result.V.3    <- pred.model%*%M.inv.3%*%t( pred.model)
      Var.pred.3    <- diag(Result.V.3)
      minimum.3 <- min(Var.pred.3)
      maximum.3 <- max(Var.pred.3)
      average.3 <- mean(Var.pred.3)
      Matrix.V.3[ii,] <- c(R.rho,maximum.3,minimum.3,average.3)     
    } 
  }
  
  # Find lim.max
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.max <- max(c(Matrix.V[,2], Matrix.V.2[,2], Matrix.V.3[,2]))
  }else{
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Matrix.V[,2], Matrix.V.2[,2]))
    }else{ 
      lim.max <- max(Matrix.V[,2])
    }
  }
  # Find lim.min
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.min <- min(c(Matrix.V[,3], Matrix.V.2[,3], Matrix.V.3[,3]))
  }else{
    if(!is.null(design.matrix.2) ){
      lim.min <- min(c(Matrix.V[,3], Matrix.V.2[,3]))
    }else{
      lim.min <- min(Matrix.V[,3])
    }
  }

  #windows(width = 5, height = 5)
  #par(mfrow=c(1,2),
  #    mai = c(1, 1, 1, 0.5),
  #    omi = c(0.2, 0.1, 0.1, 0.1),
  #    mgp = c(2, 1, 0),
  #    xpd = FALSE)
  par(mfrow=c(1,1),
      mai = c(0.75, 0.75, 0.75, 0.375),
      omi = c(0.075, 0.0375, 0.0375, 0.0375),
      mgp = c(2, 1, 0),
      xpd = FALSE)
  plot(Matrix.V[,1],Matrix.V[,2], type = "l", lwd = 2, lty = 1, 
       ylim = c(lim.min-4, lim.max + 15), col = "#2E2E2E",
       xlab = "Radius", 
       ylab = "Scaled Variance",
       panel.first = grid())
  lines(Matrix.V[,1],Matrix.V[,3], lwd = 2, col = "#2E2E2E", lty = 6)
  lines(Matrix.V[,1],Matrix.V[,4], lwd = 2, col = "#2E2E2E", lty = 3)
  
  
  legend(x = 0, y = lim.max + 16,  legend = c("Max","Avg","Min"),lty=c(1,3,6),
         lwd=c(2,2,2),col=c(1,1,1),
         inset = 0.05, bg="transparent", bty = "n")
  
  # Design 2
  if(!is.null(design.matrix.2)){
    lines(Matrix.V.2[,1],Matrix.V.2[,2], lwd = 2, col = "#EE4000", lty = 1)
    lines(Matrix.V.2[,1],Matrix.V.2[,3], lwd = 2, col = "#EE4000", lty = 6)
    lines(Matrix.V.2[,1],Matrix.V.2[,4], lwd = 2, col = "#EE4000", lty = 3)
  }
  # Design 3
  if(!is.null(design.matrix.3)){
    lines(Matrix.V.3[,1],Matrix.V.3[,2], lwd = 2, col = "#3A5FCD", lty = 1)
    lines(Matrix.V.3[,1],Matrix.V.3[,3], lwd = 2, col = "#3A5FCD", lty = 6)
    lines(Matrix.V.3[,1],Matrix.V.3[,4], lwd = 2, col = "#3A5FCD", lty = 3)
  }
  
  # Put legend Design names
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    legend(x = 0.6, y = lim.max + 16,  
           legend = c(des.names[1],des.names[2],des.names[3]),lty=c(1,1,1),
           lwd=c(2,2,2),col=c("#2E2E2E","#EE4000","#3A5FCD"),
           inset = 0.05, bg="transparent", bty = "n")
  }else{
    if(!is.null(design.matrix.2)){
      legend(x = 0.6, y = lim.max + 16,  
             legend = c(des.names[1],des.names[2]),lty=c(1,1),
             lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      legend(x = 0.6, y = lim.max + 16,  
             legend = c(des.names[1]),lty=c(1),
             lwd=c(2),col=c("#2E2E2E"),
             inset = 0.05, bg="transparent", bty = "n")
      
    }}
  
  # abline
  abline(h = ncol, col = "#8B8682",  lwd = 1)
  text(sqrt(n.var)-0.15, lim.min -3 , paste("p =",ncol),col = "#8B8682")
  
  if(add.pts == TRUE){
    # Filling points within max and min
    
    if(n.var <= 4){
      rand.sphere.2times <- matrix(c(runif(n=(4*n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))),byrow= TRUE, ncol=n.var)
      norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
      norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var)   
      radius.random      <- runif(n = 4*N.generated, min = 0, max = sqrt(n.var))
      Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
      pred.model.random  <- model.matrix( ~quad(.), Surface.pts.random)
    }
    
    if(n.var > 4){
      rand.sphere.2times <- matrix(c(runif(n=(2*n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))),byrow= TRUE, ncol=n.var)
      norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
      norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var)   
      radius.random      <- runif(n = 2*N.generated, min = 0, max = sqrt(n.var))
      Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
      pred.model.random  <- model.matrix( ~quad(.), Surface.pts.random)
    }
    
    Result.V.random    <- pred.model.random%*%M.inv%*%t(pred.model.random)
    Var.pred.random    <- diag(Result.V.random)
    points(radius.random,Var.pred.random, pch = ".", col = "#2E2E2E")
    
    # Filling points within max and min for Design #2
    if(!is.null(design.matrix.2)){
      Result.V.random.2    <- pred.model.random%*%M.inv.2%*%t(pred.model.random)
      Var.pred.random.2    <- diag(Result.V.random.2)
      points(radius.random,Var.pred.random.2, pch = ".", col = "#EE4000")
    }
    # Filling points within max and min for Design #3
    if(!is.null(design.matrix.3)){
      Result.V.random.3    <- pred.model.random%*%M.inv.3%*%t(pred.model.random)
      Var.pred.random.3    <- diag(Result.V.random.3)
      points(radius.random,Var.pred.random.3, pch = ".", col = "#3A5FCD")
    }
  }
  
  par(mfrow=c(1,1))
  
  # Return Values
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    return(list(design.1 = Matrix.V, design.2 = Matrix.V.2, 
                design.3 = Matrix.V.3))
  }else{
    if(!is.null(design.matrix.2)){
      return(list(design.1 = Matrix.V, design.2 = Matrix.V.2))
    }else{
      return(Matrix.V)      
    }
    
  }
}


######################### fds.sphere function  ########################

fds.sphere <- function(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
                       des.names = c("Design 1","Design 2","Design 3"),
                       scale = TRUE){
  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.var <- ncol(design.matrix)
  
  if(n.var > 7){stop("The maximum number of design factors allowed is 7 \n")}
  
  
  if(scale == TRUE){ 
    temp.radii <- apply(design.matrix, MARGIN = 1, norm2)
    scale.factor <- sqrt(n.var)/max(temp.radii)
    design.matrix <- scale.factor*design.matrix
  }
  
  design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
  N.generated <- 15*n.var*(2^(n.var))
  model.X <- model.matrix( ~quad(.) , design.matrix)
  nrow <- nrow(model.X)
  ncol <- ncol(model.X)
  M.inv<- nrow*solve(t(model.X)%*%model.X)
  
  # Design 2
  if(!is.null(design.matrix.2)){
    n.var.2 <- ncol(design.matrix.2)
    if(scale == TRUE){ 
      temp.radii.2 <- apply(design.matrix.2, MARGIN = 1, norm2)
      scale.factor.2 <- sqrt(n.var.2)/max(temp.radii.2)
      design.matrix.2 <- scale.factor.2*design.matrix.2
    }
    design.matrix.2<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.2))),byrow=FALSE,ncol=n.var.2))  
    model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
    nrow.2 <- nrow(model.X.2)
    ncol.2 <- ncol(model.X.2)
    if(ncol != ncol.2){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
  } # Design 2
  # Design 3
  if(!is.null(design.matrix.3)){
    
    n.var.3 <- ncol(design.matrix.3)
    if(scale == TRUE){ 
      temp.radii.3 <- apply(design.matrix.3, MARGIN = 1, norm2)
      scale.factor.3 <- sqrt(n.var.3)/max(temp.radii.3)
      design.matrix.3 <- scale.factor.3*design.matrix.3
    }
    design.matrix.3<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.3))),byrow=FALSE,ncol=n.var.3))  
    model.X.3 <- model.matrix( ~quad(.) , design.matrix.3)
    nrow.3 <- nrow(model.X.3)
    ncol.3 <- ncol(model.X.3)
    if(ncol != ncol.3){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.3<- nrow.3*solve(t(model.X.3)%*%model.X.3)   
  } # Design 3  
  
  # Filling points within max and min
  
  if(n.var <= 4){
    rand.sphere.2times <- matrix(c(runif(n=(4*n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))),byrow= TRUE, ncol=n.var)
    norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
    norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var)   
    radius.random      <- runif(n = 4*N.generated, min = 0, max = sqrt(n.var))
    Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
    pred.model.random  <- model.matrix( ~quad(.), Surface.pts.random)
  }
  
  if(n.var > 4){
    rand.sphere.2times <- matrix(c(runif(n=(1*n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))),byrow= TRUE, ncol=n.var)
    norm2.rand.2times <- apply(rand.sphere.2times, MARGIN = 1, norm2)
    norm2.rand.rep.2times <- matrix(rep(norm2.rand.2times, each = n.var),byrow= TRUE, ncol=n.var)   
    radius.random      <- runif(n = 1*N.generated, min = 0, max = sqrt(n.var))
    Surface.pts.random <- radius.random*as.data.frame(rand.sphere.2times/norm2.rand.rep.2times)
    pred.model.random  <- model.matrix( ~quad(.), Surface.pts.random)
  }
  
  Result.V.random    <- pred.model.random%*%M.inv%*%t(pred.model.random)
  Var.pred.random    <- diag(Result.V.random)
  
  # Design #2
  if(!is.null(design.matrix.2)){
    Result.V.random.2    <- pred.model.random%*%M.inv.2%*%t(pred.model.random)
    Var.pred.random.2    <- diag(Result.V.random.2)
  }
  #Design #3
  if(!is.null(design.matrix.3)){
    Result.V.random.3    <- pred.model.random%*%M.inv.3%*%t(pred.model.random)
    Var.pred.random.3    <- diag(Result.V.random.3)
  }
  
  # Find lim.max
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.max <- max(c(Var.pred.random, Var.pred.random.2, Var.pred.random.3))
  }else{
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Var.pred.random, Var.pred.random.2))
    }else{ 
      lim.max <- max(Var.pred.random)
    }
  }
  # Find lim.min
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.min <- min(c(Var.pred.random, Var.pred.random.2, Var.pred.random.3))
  }else{
    if(!is.null(design.matrix.2)){
      lim.min <- min(c(Var.pred.random, Var.pred.random.2))
    }else{ 
      lim.min <- min(Var.pred.random)
    }
  }
  
  #windows(width = 5, height = 5)
  #par(mfrow=c(1,2),
  #    mai = c(1, 1, 1, 0.5),
  #    omi = c(0.2, 0.1, 0.1, 0.1),
  #    mgp = c(2, 1, 0),
  #    xpd = FALSE)
  par(mfrow=c(1,1),
      mai = c(0.75, 0.75, 0.75, 0.375),
      omi = c(0.075, 0.0375, 0.0375, 0.0375),
      mgp = c(2, 1, 0),
      xpd = FALSE)
  
  # Make FDS plots
  plot(seq(0, 1, 0.01),quantile(Var.pred.random, probs = seq(0, 1, 0.01), type = 4),
       type = "l",lwd = 2, col = "#2E2E2E",
       xlab = "Fraction of design space",ylab = "Scaled Variance",
       ylim = c(lim.min-4, lim.max + 15), panel.first = grid())
  
  # Make FDS plots for Design 2
  if(!is.null(design.matrix.2)){
    lines(seq(0, 1, 0.01),quantile(Var.pred.random.2, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#EE4000")
  }
  # Make FDS plots for Design 3
  if(!is.null(design.matrix.3)){
    lines(seq(0, 1, 0.01),quantile(Var.pred.random.3, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#3A5FCD")
  }
  
  # Put legend in VDG plots
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    legend(x = 0, y = lim.max + 16,  
           legend = c(des.names[1],des.names[2],des.names[3]),lty=c(1,1,1),
           lwd=c(2,2,2),col=c("#2E2E2E","#EE4000","#3A5FCD"),
           inset = 0.05, bg="transparent", bty = "n")
  }else{
    if(!is.null(design.matrix.2)){
      legend(x = 0, y = lim.max + 16,  
             legend = c(des.names[1],des.names[2]),lty=c(1,1),
             lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      legend(x = 0, y = lim.max + 16,  
             legend = c(des.names[1]),lty=c(1),
             lwd=c(2),col=c("#2E2E2E"),
             inset = 0.05, bg="transparent", bty = "n")
    }}
  # abline
  abline(h = ncol, col = "#8B8682",  lwd = 1)
  text(sqrt(n.var)-0.15, lim.min -3 , paste("p =",ncol),col = "#8B8682")
}

#########################  cpv function ####################

cpv <- function(design.matrix, add.pts = FALSE){
  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.var <- ncol(design.matrix)
  if(n.var > 6){stop("The maximum number of design factors allowed is 6 \n")}
  critical.R <- numeric(0)
  
  for(ii in 0:n.var){
    critical.R <- append(critical.R,sqrt(ii))
  }
  LR <- length(critical.R)
  
  if(n.var > 4){
    for(jj in 1:(LR-1)){
      critical.R <- 
        append(critical.R,
               seq(from = critical.R[jj], to = critical.R[jj+1], length.out = 3))
    }
  }else{
    for(jj in 1:(LR-1)){
      critical.R <- 
        append(critical.R,
               seq(from = critical.R[jj], to = critical.R[jj+1], length.out = 8))
    }
  }
  
  critical.R<- unique(sort(critical.R))
  
  if(n.var <= 4){
    n.grid <- length(critical.R)
    N.generated <- 7000
    model.X <- model.matrix( ~quad(.) , design.matrix)
    nrow <- nrow(model.X)
    ncol <- ncol(model.X)
    M.inv<- nrow*solve(t(model.X)%*%model.X)
    
    Matrix.V <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V) <- c("Radius","  max.V","  min.V"," average.V")
    # Generate points on a sphere
    set.seed(1234567)
    rand.sphere <- matrix(c(runif(n=(n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))),byrow= TRUE, ncol=n.var)
    norm2.rand <- apply(rand.sphere, MARGIN = 1, norm2)
    norm2.rand.rep <- matrix(rep(norm2.rand, each = n.var),byrow= TRUE, ncol=n.var)
    generated.pts<- as.data.frame(rand.sphere/norm2.rand.rep)
    
    
    count <- 0
    for(R in critical.R){
      count <- count + 1
      Surface.pts <- R*generated.pts
      Conbind.col<- cbind(Surface.pts,rowSums(abs(Surface.pts[,])>1))
      Surface.pts <- Conbind.col[which(Conbind.col[,(n.var+1)] == 0),]
      Surface.pts <- Surface.pts[,-(n.var+1)]    
      
      if( dim(Surface.pts)[1] == 0){ break }
      pred.model  <- model.matrix( ~quad(.), Surface.pts)
      Result.V    <- pred.model%*%M.inv%*%t( pred.model)
      Var.pred    <- diag(Result.V)
      
      minimum <- min(Var.pred)
      maximum <- max(Var.pred)
      average <- mean(Var.pred)
      Matrix.V[count,] <- c(R,maximum,minimum,average) 
      
    } 
    
    Matrix.V<- matrix(Matrix.V[is.na(Matrix.V) == FALSE],byrow = FALSE, ncol = 4)
    colnames(Matrix.V) <- c("Radius","  max.V","  min.V"," average.V")
    
    R <- critical.R[n.grid]
    Surface.pts<- as.data.frame(matrix(sample(c(-1,1),size = n.var*50,replace = TRUE),byrow = TRUE, ncol = n.var))
    
    pred.model  <- model.matrix( ~quad(.), Surface.pts)
    Result.V    <- pred.model%*%M.inv%*%t( pred.model)
    Var.pred    <- diag(Result.V)
    minimum <- min(Var.pred)
    maximum <- max(Var.pred)
    average <- mean(Var.pred)
    Matrix.V<- rbind(Matrix.V,c(R,maximum,minimum,average) )
    
    lim.max <- max(Matrix.V[,2],na.rm = TRUE)
    lim.min <- min(Matrix.V[,3],na.rm = TRUE)
    
    # Find lim.max
    lim.max <- max(Matrix.V[,2])
    # Find lim.min
    lim.min <- min(Matrix.V[,3])
    
    
    #windows(width = 5, height = 5)
    par(mfrow=c(1,1),
        mai = c(0.75, 0.75, 0.75, 0.375),
        omi = c(0.075, 0.0375, 0.0375, 0.0375),
        mgp = c(2, 1, 0),
        xpd = FALSE)
    
    plot(Matrix.V[,1],Matrix.V[,2], type = "l", lwd = 2, lty = 4, 
         ylim = c(lim.min-4, lim.max + 13), col = "#2E2E2E",
         xlab = "Radius", 
         ylab = "Scaled Variance",
         panel.first = grid())
    lines(Matrix.V[,1],Matrix.V[,3], lwd = 2, col = "#2E2E2E", lty = 2)
    legend(x = 0, y = lim.max + 14,  legend = c("Max","Min"),lty=c(4,2),
           lwd=c(2,2),col=c(1,1),
           inset = 0.05, bg="transparent", bty = "n")
    abline(h = ncol, col = "#8B8682",  lwd = 1)
    text(sqrt(n.var)-0.15, lim.min -3 , paste("p =",ncol),col = "#8B8682")
    
    
    # Filling points
    if(add.pts == TRUE){
      Cuboidal.pts.random  <- matrix(runif(n = N.generated*n.var, min = -1, max = 1),
                                     byrow = TRUE, ncol = n.var)      
      Cuboidal.pts.random  <- as.data.frame(Cuboidal.pts.random)
      pred.model.random  <- model.matrix( ~quad(.), Cuboidal.pts.random)
      radius.random<- apply(Cuboidal.pts.random, MARGIN = 1, norm2)
      Result.V.random    <- pred.model.random%*%M.inv%*%t(pred.model.random)
      Var.pred.random    <- diag(Result.V.random)
      points(radius.random,Var.pred.random, pch = ".", col = "#2E2E2E")
    }
    
    # Return Values
    return(Matrix.V)   
  }
  if(n.var > 4 && n.var <= 6){
    n.grid <- length(critical.R)
    N.generated <- 15000
    model.X <- model.matrix( ~quad(.) , design.matrix)
    nrow <- nrow(model.X)
    ncol <- ncol(model.X)
    M.inv<- nrow*solve(t(model.X)%*%model.X)
    
    Matrix.V <- matrix(NA, ncol = n.grid, nrow = N.generated)
    
    # Generate points on a sphere
    set.seed(1234567)
    rand.sphere <- matrix(c(runif(n=(n.var*N.generated),min=-sqrt(n.var),max=sqrt(n.var))),byrow= TRUE, ncol=n.var)
    norm2.rand <- apply(rand.sphere, MARGIN = 1, norm2)
    norm2.rand.rep <- matrix(rep(norm2.rand, each = n.var),byrow= TRUE, ncol=n.var)
    generated.pts<- as.data.frame(rand.sphere/norm2.rand.rep)
    
    
    count <- 0
    for(R in critical.R[1:(n.grid - 1)]){
      count <- count + 1
      Surface.pts <- R*generated.pts
      Conbind.col<- cbind(Surface.pts,rowSums(abs(Surface.pts[,])>1))
      Surface.pts <- Conbind.col[which(Conbind.col[,(n.var+1)] == 0),]
      Surface.pts <- Surface.pts[,-(n.var+1)]    
      
      #if( dim(Surface.pts)[1] == 0){ break }
      pred.model  <- model.matrix( ~quad(.), Surface.pts)
      Result.V    <- pred.model%*%M.inv%*%t( pred.model)
      Answer <- diag(Result.V)
      Matrix.V[1:length(Answer),count]    <- Answer
    }
    
    R <- critical.R[n.grid]
    Surface.pts<- as.data.frame(matrix(sample(c(-1,1),size = n.var*50,replace = TRUE),byrow = TRUE, ncol = n.var))
    
    pred.model  <- model.matrix( ~quad(.), Surface.pts)
    Result.V    <- pred.model%*%M.inv%*%t( pred.model)
    Var.pred    <- diag(Result.V)
    Answer <- diag(Result.V)
    Matrix.V[1:length(Answer),n.grid]    <- Answer
    #windows(width = 5, height = 5)
    par(mfrow=c(1,1),
        mai = c(0.75, 0.75, 0.75, 0.375),
        omi = c(0.075, 0.0375, 0.0375, 0.0375),
        mgp = c(2, 1, 0),
        xpd = FALSE)
    boxplot(Matrix.V, use.cols = TRUE, outline = TRUE, range = 9.5, col = "lightgray",
            whisklty = 1, staplelty = 0, ylab = "Scaled Variace",
            names = round(critical.R,2), xlab = "Radius", lwd = 2)
  }
  
}

########################### fds.cube function  ########################

fds.cube <- function(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
                     des.names = c("Design 1","Design 2","Design 3")){
  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.var <- ncol(design.matrix)
  
  if(n.var > 6){stop("The maximum number of design factors allowed is 6 \n")}
  
  
  design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
  model.X <- model.matrix( ~quad(.) , design.matrix)
  nrow <- nrow(model.X)
  ncol <- ncol(model.X)
  M.inv<- nrow*solve(t(model.X)%*%model.X)
  
  # Design 2
  if(!is.null(design.matrix.2)){
    n.var.2 <- ncol(design.matrix.2)
    
    design.matrix.2<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.2))),byrow=FALSE,ncol=n.var.2))  
    model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
    nrow.2 <- nrow(model.X.2)
    ncol.2 <- ncol(model.X.2)
    if(ncol != ncol.2){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
  } # Design 2
  # Design 3
  if(!is.null(design.matrix.3)){
    
    n.var.3 <- ncol(design.matrix.3)
    
    design.matrix.3<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.3))),byrow=FALSE,ncol=n.var.3))  
    model.X.3 <- model.matrix( ~quad(.) , design.matrix.3)
    nrow.3 <- nrow(model.X.3)
    ncol.3 <- ncol(model.X.3)
    if(ncol != ncol.3){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.3<- nrow.3*solve(t(model.X.3)%*%model.X.3)   
  } # Design 3 
  
  
  
  
  if(n.var <= 4){
    Cuboidal.pts.random  <- matrix(runif(n = 6000*n.var, min = -1, max = 1),
                                   byrow = TRUE, ncol = n.var)      
    Cuboidal.pts.random  <- as.data.frame(Cuboidal.pts.random)
    pred.model.random  <- model.matrix( ~quad(.), Cuboidal.pts.random)
    radius.random<- apply(Cuboidal.pts.random, MARGIN = 1, norm2)
  }
  
  if(n.var > 4){
    Cuboidal.pts.random  <- matrix(runif(n = 11000*n.var, min = -1, max = 1),
                                   byrow = TRUE, ncol = n.var)      
    Cuboidal.pts.random  <- as.data.frame(Cuboidal.pts.random)
    pred.model.random  <- model.matrix( ~quad(.), Cuboidal.pts.random)
    radius.random<- apply(Cuboidal.pts.random, MARGIN = 1, norm2)
  }
  
  Result.V.random    <- pred.model.random%*%M.inv%*%t(pred.model.random)
  Var.pred.random    <- diag(Result.V.random)
  
  # Design #2
  if(!is.null(design.matrix.2)){
    Result.V.random.2    <- pred.model.random%*%M.inv.2%*%t(pred.model.random)
    Var.pred.random.2    <- diag(Result.V.random.2)
  }
  #Design #3
  if(!is.null(design.matrix.3)){
    Result.V.random.3    <- pred.model.random%*%M.inv.3%*%t(pred.model.random)
    Var.pred.random.3    <- diag(Result.V.random.3)
  }
  
  # Find lim.max
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.max <- max(c(Var.pred.random, Var.pred.random.2, Var.pred.random.3))
  }else{
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Var.pred.random, Var.pred.random.2))
    }else{ 
      lim.max <- max(Var.pred.random)
    }
  }
  # Find lim.min
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.min <- min(c(Var.pred.random, Var.pred.random.2, Var.pred.random.3))
  }else{
    if(!is.null(design.matrix.2)){
      lim.min <- min(c(Var.pred.random, Var.pred.random.2))
    }else{ 
      lim.min <- min(Var.pred.random)
    }
  }
  
  #windows(width = 5, height = 5)
  par(mfrow=c(1,1),
      mai = c(0.75, 0.75, 0.75, 0.375),
      omi = c(0.075, 0.0375, 0.0375, 0.0375),
      mgp = c(2, 1, 0),
      xpd = FALSE)
  
  # Make FDS plots
  plot(seq(0, 1, 0.01),quantile(Var.pred.random, probs = seq(0, 1, 0.01), type = 4),
       type = "l",lwd = 2, col = "#2E2E2E",
       xlab = "Fraction of design space",ylab = "Scaled Variance",
       ylim = c(lim.min-4, lim.max + 15), panel.first = grid())
  
  # Make FDS plots for Design 2
  if(!is.null(design.matrix.2)){
    lines(seq(0, 1, 0.01),quantile(Var.pred.random.2, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#EE4000")
  }
  # Make FDS plots for Design 3
  if(!is.null(design.matrix.3)){
    lines(seq(0, 1, 0.01),quantile(Var.pred.random.3, probs = seq(0, 1, 0.01), type = 4), lwd = 2, col = "#3A5FCD")
  }
  
  # Put legend in VDG plots
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    legend(x = 0, y = lim.max + 16,  
           legend = c(des.names[1],des.names[2],des.names[3]),lty=c(1,1,1),
           lwd=c(2,2,2),col=c("#2E2E2E","#EE4000","#3A5FCD"),
           inset = 0.05, bg="transparent", bty = "n")
  }else{
    if(!is.null(design.matrix.2)){
      legend(x = 0, y = lim.max + 16,  
             legend = c(des.names[1],des.names[2]),lty=c(1,1),
             lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      legend(x = 0, y = lim.max + 16,  
             legend = c(des.names[1]),lty=c(1),
             lwd=c(2),col=c("#2E2E2E"),
             inset = 0.05, bg="transparent", bty = "n")
    }}
  # abline
  abline(h = ncol, col = "#8B8682",  lwd = 1)
  text(0.90, lim.min -3 , paste("p =",ncol),col = "#8B8682")
}

##########################  hyper.vdg #######################

hyperarcs.vdg<- function(design.matrix, design.matrix.2 = NULL, design.matrix.3 = NULL, 
                     des.names = c("Design 1","Design 2","Design 3")){
  
  norm2 <- function(x){return(sqrt(sum(x^2)))}
  shuffle.fun <- function(row.vec){
    num.var <- length(row.vec)
    row.vec[shuffle(num.var)]
  }
  
  n.var <- ncol(design.matrix)
  if(n.var > 6){stop("The maximum number of design factors allowed is 6 \n")}
  n.grid <- 15
  design.matrix<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix))),byrow=FALSE,ncol=n.var))  
  
  model.X <- model.matrix( ~quad(.) , design.matrix)
  nrow <- nrow(model.X)
  ncol <- ncol(model.X)
  M.inv<- nrow*solve(t(model.X)%*%model.X)
  
  # Design 2
  if(!is.null(design.matrix.2)){
    n.var.2 <- ncol(design.matrix.2)
    
    design.matrix.2<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.2))),byrow=FALSE,ncol=n.var.2))  
    model.X.2 <- model.matrix( ~quad(.) , design.matrix.2)
    nrow.2 <- nrow(model.X.2)
    ncol.2 <- ncol(model.X.2)
    if(ncol != ncol.2){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.2<- nrow.2*solve(t(model.X.2)%*%model.X.2)   
  } # Design 2
  # Design 3
  if(!is.null(design.matrix.3)){
    
    n.var.3 <- ncol(design.matrix.3)
    design.matrix.3<- as.data.frame(matrix(as.numeric(paste(unlist(design.matrix.3))),byrow=FALSE,ncol=n.var.3))  
    model.X.3 <- model.matrix( ~quad(.) , design.matrix.3)
    nrow.3 <- nrow(model.X.3)
    ncol.3 <- ncol(model.X.3)
    if(ncol != ncol.3){
      stop("Designs need to have the same number of design factors")
    }
    M.inv.3<- nrow.3*solve(t(model.X.3)%*%model.X.3)   
  } # Design 3  
  
  Matrix.V.A <- matrix(numeric(0), ncol = 4, nrow = n.grid )
  colnames(Matrix.V.A) <- c("Hypercube Radius","  max.V","  min.V"," average.V")
  
  if(!is.null(design.matrix.2)){
    Matrix.V.A.2 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V.A.2) <- c("Hypercube Radius","  max.V","  min.V"," average.V")
  }# Design 2
  
  if(!is.null(design.matrix.3)){
    Matrix.V.A.3 <- matrix(numeric(0), ncol = 4, nrow = n.grid)
    colnames(Matrix.V.A.3) <- c("Hypercube Radius","  max.V","  min.V"," average.V")
  }# Design 3
  
  radius <- seq(from = 0, to = 1, length.out = n.grid)
  N.generated <- 25*n.var*2^(n.var)
  
  count<- 0
  for(R in radius){
    count <- count + 1
    set.seed(1234567)
    mat.shuffle <- as.matrix(shuffleSet(n = n.var, nset = N.generated,check = FALSE))
    
    full.basket <- matrix(c(sample(c(-R,R), size = N.generated*(n.var-2), replace = TRUE)),
                          ncol = (n.var-2), byrow = TRUE)
    mat.before <- cbind(full.basket,runif(n = N.generated, min = -R, max = R)
                        ,runif(n = N.generated, min =  -R, max = R) )    
    Cuboidal.pts <- as.data.frame(t(apply(mat.before, MARGIN = 1, shuffle.fun)))
    
    pred.model  <- model.matrix( ~quad(.), Cuboidal.pts)
    Result.V    <- pred.model%*%M.inv%*%t( pred.model)
    Var.pred    <- diag(Result.V)
    minimum <- min(Var.pred)
    maximum <- max(Var.pred)
    average <- mean(Var.pred)
    Matrix.V.A[count,] <- c(R,maximum,minimum,average) 
    
    ###################################
    # Design 2
    if(!is.null(design.matrix.2)){
      Result.V.2    <- pred.model%*%M.inv.2%*%t( pred.model)
      Var.pred.2    <- diag(Result.V.2)
      minimum.2 <- min(Var.pred.2)
      maximum.2 <- max(Var.pred.2)
      average.2 <- mean(Var.pred.2)
      Matrix.V.A.2[count,] <- c(R,maximum.2,minimum.2,average.2)     
    } 
    # Design 3
    if(!is.null(design.matrix.3)){
      Result.V.3    <- pred.model%*%M.inv.3%*%t( pred.model)
      Var.pred.3    <- diag(Result.V.3)
      minimum.3 <- min(Var.pred.3)
      maximum.3 <- max(Var.pred.3)
      average.3 <- mean(Var.pred.3)
      Matrix.V.A.3[count,] <- c(R,maximum.3,minimum.3,average.3)      
    } 
  }
  # Find lim.max
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.max <- max(c(Matrix.V.A[,2], Matrix.V.A.2[,2], Matrix.V.A.3[,2]))
  }else{
    if(!is.null(design.matrix.2)){
      lim.max <- max(c(Matrix.V.A[,2], Matrix.V.A.2[,2]))
    }else{ 
      lim.max <- max(Matrix.V.A[,2])
    }
  }
  # Find lim.min
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    lim.min <- min(c(Matrix.V.A[,3], Matrix.V.A.2[,3], Matrix.V.A.3[,3]))
  }else{
    if(!is.null(design.matrix.2) ){
      lim.min <- min(c(Matrix.V.A[,3], Matrix.V.A.2[,3]))
    }else{
      lim.min <- min(Matrix.V.A[,3])
    }
  }
  
  #windows(width = 5, height = 5)
  par(mfrow=c(1,1),
      mai = c(0.75, 0.75, 0.75, 0.375),
      omi = c(0.075, 0.0375, 0.0375, 0.0375),
      mgp = c(2, 1, 0),
      xpd = FALSE)
  plot(Matrix.V.A[,1],Matrix.V.A[,2], type = "l", lwd = 2, lty = 1, 
       ylim = c(lim.min-4, lim.max + 15), col = "#2E2E2E",
       xlab = "Hypercube Radius", 
       ylab = "Scaled Variance",
       panel.first = grid())
  lines(Matrix.V.A[,1],Matrix.V.A[,3], lwd = 2, col = "#2E2E2E", lty = 6)
  lines(Matrix.V.A[,1],Matrix.V.A[,4], lwd = 2, col = "#2E2E2E", lty = 3)
  
  
  legend(x = 0, y = lim.max + 16,  legend = c("Max","Avg","Min"),lty=c(1,3,6),
         lwd=c(2,2,2),col=c(1,1,1),
         inset = 0.05, bg="transparent", bty = "n")
  
  # Design 2
  if(!is.null(design.matrix.2)){
    lines(Matrix.V.A.2[,1],Matrix.V.A.2[,2], lwd = 2, col = "#EE4000", lty = 1)
    lines(Matrix.V.A.2[,1],Matrix.V.A.2[,3], lwd = 2, col = "#EE4000", lty = 6)
    lines(Matrix.V.A.2[,1],Matrix.V.A.2[,4], lwd = 2, col = "#EE4000", lty = 3)
  }
  # Design 3
  if(!is.null(design.matrix.3)){
    lines(Matrix.V.A.3[,1],Matrix.V.A.3[,2], lwd = 2, col = "#3A5FCD", lty = 1)
    lines(Matrix.V.A.3[,1],Matrix.V.A.3[,3], lwd = 2, col = "#3A5FCD", lty = 6)
    lines(Matrix.V.A.3[,1],Matrix.V.A.3[,4], lwd = 2, col = "#3A5FCD", lty = 3)
  }
  # Put legend Design names
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    legend(x = 0.4, y = lim.max + 16,  
           legend = c(des.names[1],des.names[2],des.names[3]),lty=c(1,1,1),
           lwd=c(2,2,2),col=c("#2E2E2E","#EE4000","#3A5FCD"),
           inset = 0.05, bg="transparent", bty = "n")
  }else{
    if(!is.null(design.matrix.2)){
      legend(x = 0.4, y = lim.max + 16,  
             legend = c(des.names[1],des.names[2]),lty=c(1,1),
             lwd=c(2,2),col=c("#2E2E2E","#EE4000"),
             inset = 0.05, bg="transparent", bty = "n")
    }else{
      legend(x = 0.4, y = lim.max + 16,  
             legend = c(des.names[1]),lty=c(1),
             lwd=c(2),col=c("#2E2E2E"),
             inset = 0.05, bg="transparent", bty = "n")
      
    }}
  
  # Return Values
  if(!is.null(design.matrix.2) && !is.null(design.matrix.3)){
    return(list(design.1 = Matrix.V.A, design.2 = Matrix.V.A.2, 
                design.3 = Matrix.V.A.3))
  }else{
    if(!is.null(design.matrix.2)){
      return(list(design.1 = Matrix.V.A, design.2 = Matrix.V.A.2))
    }else{
      return(Matrix.V.A)      
    }
    
  }
}

