
# Part I
dist <- function (x,y){ # sum of | x_i - y_i | from 1 to d
  sum <- 0 #Initial sum
  for(i in 1: length(x)){
    sum <- sum + abs(x[i] - y[i]) # loop from 1 to dimension d 
  }
  return (sum)
  #return the distance
}


dist_mat <- function(X , y){ #distance matrix of X_row and vector y
  n = nrow(X)
  v = numeric(n)
  for (i in 1: n){
    v[i] = dist(X[i, ], y)
  }
  
  return (v) # v = {dist(X_1:, y),dist(X_2:, y),...,dist(X_n:, y)}.
}



dist_mat_fast <- function(X, y){ #distance matrix of X_row and vector y
  y_matrix <- matrix(y, nrow = nrow(X), ncol = length(y), byrow = T)
  diff_matrix <- abs(X - y_matrix)
  v <- rowSums(diff_matrix)
  
  return (v) # v = {dist(X_1:, y),dist(X_2:, y),...,dist(X_n:, y)}.
}



# Part II
w_fun <- function(z, h){ # vector z and positive scalar h
  
  w <- exp(-(z^2)/ (h^2))
  
  return (w) # w_i = exp(−(z_i)^2 /h_2), for i = 1,...,n.
}

w_mean <- function(y, X, h){
  n <- nrow(X) #number of row of matrix X
  z <- dist_mat(X, y) # distance matrix of X_row and vector y
  w <- w_fun(z, h) # w_fun function
  
  sum_wX <- 0
  for(i in 1: n){
    sum_wX <- sum_wX + (w[i] * X[i, ]) 
  }
  sum_w <- sum(w)
  y1 <- sum_wX / sum_w
  
  return (y1) # y1 = sum of (w_i * X_i) / (w_l) for i = 1,..,n  and l = 1,...,n
}             # where w_i = exp( -dist(y, X_i)^2 / h^2 )


w_mean_vett <- function(Y, X ,h){
  Y1_matrix <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  
  for (i in 1:nrow(Y)){
    Y1_matrix[i, ] <- w_mean(Y[i, ], X, h)
  }
  
  return(Y1_matrix)  # Y1_i = w_mean(Y_i, X, h) 
}                    # where i = 1,...,n  Y is n x d matrix


# Part III

#Algorithm that will shift the rows of X
phase_1 <- function(X, h, maxit, delta){
  X1 <- X # Set matrix X to be X1
  n = nrow(X) #number of row of matrix X
  
  for (i in 1: maxit){
    X1_old = X1 #Set X1 to be X1_old
    X1 <- w_mean_vett(X1_old, X1_old, h) #Shifted X1
    row_distance <- numeric (n) 
    
    for (j in 1:n){
      row_distance <- dist(X1[j, ], X1_old[j, ]) #distance between shifted X1 and X1_old
    } # data point have reached the max number of iterations
    if(all(row_distance < delta) ){
      break
    } # data point shifted until the shifted data point have converged
  } 
  return(X1) # returning the shifted data points of X1.
} 


#Phase 2
#Shifted data points are merged into few clusters
#The corresponding original data points are attributed to the corresponding cluster
phase_2 <- function(Xt, eps){
  n <- nrow(Xt) #number of row of matrix Xt
  C <- list () #initialize C as the empty set
  cl <- numeric(n)
  
  for(i in 1: n){
    if (length(C) > 0){
      for( j in 1: length(C)){
        if( dist( Xt[i], C[[j]] ) < eps ){  # if distance of X_i and C_j is less then eps 
          cl[i] <- j # assume X_i belongs to the j-th cluster
          break
        }
      }
    }
    if (cl[i] == 0){ 
      C[[length(C) + 1]] <- Xt[i, ] # #add Xt_i to C
      cl[i] <- length(C) #set cl_i to the number of elements of C
    }
  }
  return(list(C = C, cl = cl)) #Return the list
}



# Function that used to plot the data using a scatter plot
plot_clustering <- function(X, C, cl){
  if (length(C) > 0) {
    # Convert list of cluster centers to a matrix 
    C_matrix <- matrix(0, nrow = length(C), ncol = length(C[[1]])) 
   
     for( i in 1:length(C) ){
      C_matrix[i, ] <- C[[i]]
    }
    
    plot(X, col = cl, pch = 20, xlab = "X1", ylab = "X2", main = "Clustering Result")
    points(C_matrix, col = 1:length(C), pch = 8, cex = 1, lwd = 2)
    } 
  else {
    plot(X, col = cl, pch = 20, xlab = "X1", ylab = "X2", main = "Clustering Result")
  }
  
}

#Test to ensure the function work on matrix and vector
set.seed(416)
X <- matrix(rnorm(10*3), 10, 3)
y <- rnorm(3)
dist_mat(X, y) #distance of matrix X and vector y
dist_mat_fast(X, y) #distance of matrix X and vector y


#Code to generate some simulated data
set.seed(123)
n <- 150
X <- rbind(
  cbind(rnorm(n / 3, mean = 0, sd = 0.3), rnorm(n / 3, mean = 0, sd = 0.3)),
  cbind(rnorm(n / 3, mean = 1, sd = 0.3), rnorm(n / 3, mean = 1, sd = 0.3)),
  cbind(rnorm(n / 3, mean = 2, sd = 0.3), rnorm(n / 3, mean = 0, sd = 0.3))
)


#Testing function using the following code
#Should produce 4 plots, showing the effect of h
par(mfrow = c(2, 2)) 
kk <- 1
for(ii in c(0.1, 1, 2, 10)){
  plot(X)
  mns <- w_mean_vett(Y = X, X = X, h = ii)
  points(mns, col = 2, pch = 2)
}


#Perform the clustering point and plot the result
Xt <- phase_1(X = X, h = 0.5, maxit = 300, delta = 1e-3)
clustering <- phase_2(Xt = Xt, eps = 1e-2)
plot_clustering(X = X, C = clustering$C, cl = clustering$cl)


# Part IV
data("iris") #Load the iris data set
variables <- colnames(iris)[1:4] #First four variables
iris_pairs <- list()


#Create a list of 6 new data sets from the iris data set.
index <- 1 
for (i in 1:(length(variables) - 1)) {
  for (j in (i + 1):length(variables)) {
    iris_pairs[[index]] <- iris[ , c(variables[i], variables[j])]
    index <- index + 1
  }  #Picking every pairs of variables from first four variables
}


#Visualize the assignments obtained from the clustering algorithm
par(mfrow = c(2,3)) 
for(i in 1: length(iris_pairs)){
  current_data <- as.matrix(iris_pairs[[i]]) #convert a list to a matrix
  
  Xt <- phase_1(X = current_data, h = 1, maxit = 300, delta = 1e-3)
  clustering <- phase_2(Xt = Xt, eps = 1e-2)
  plot_clustering(X = current_data, C = clustering$C, cl = clustering$cl)
  #Plot the results
}
