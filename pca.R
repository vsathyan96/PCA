library(stats)
#Read file from csv file
data <- read.csv("heart.csv", head=T)
target <- data[,14]
m <- length(data[1,])-1
n <- length(data[,1])
#Exclude target variable
data <- data[,1:m]
datanum <-as.matrix(data)

#Extract mean, std dev for reconstruction
mu <- apply(datanum, 2, mean)
stdev <- apply(datanum,2,sd)

#Preprocess data [normalize/feature scaling]
datanum <- scale(datanum,center = TRUE,scale = TRUE)

#Compute PCA
retVar = 0.90 # 90% Variance retained
pca <- function(X,retVar) {
  Sigma <- (t(X)%*%X)*(1/m)
  decomp <- svd(Sigma,nu=n,nv=m,LINPACK = FALSE)
  D <- decomp$d
  U <- decomp$u
  V <- decomp$v
  k=0
  z=data.frame(matrix(vector(),303,11))
  #Choose k (number of principal components)
  s <- sum(D)
  for (i in 1:m) {
    val <- sum(D[1:i])/s
    if(val>=retVar) {
      k=i
      break
    }
  }
  
  #Reduce dimension
  Ureduce <- U[,1:k]
  z <- t(t(Ureduce)%*%t(X))
  mylist <- list("z" = z,"Ureduce"=Ureduce)
  return(mylist)
}

pcaData <- pca(datanum,retVar)

#Reconstruct from compressed version
Xapprox <- t(pcaData$Ureduce%*%t(pcaData$z))
colnames(Xapprox) <- colnames(datanum)
for (i in 1:m) {
  Xapprox[,i] <- (Xapprox[,i]*stdev[i])+mu[i]
}

