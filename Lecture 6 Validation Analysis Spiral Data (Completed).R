rm(list = ls())
# A Function for creating continuous colour palette. I.e., colour gradient bins and
# associated colours:
library(colorspace)
color.gradient = function(x, colors=c('magenta','white','lightblue'), colsteps=50)
{
  colpal = colorRampPalette(colors)
  return( colpal(colsteps)[ findInterval(x, seq(min(x),max(x), length=colsteps)) ] )
}

# Read in the data and populate a design-matrix + response matrix:
dat = read.table('Spiral_Data_2022.txt',h = T)
X   = as.matrix(dat[,-3])
Y   = as.matrix((dat[,3]==1)*1) # Transform to 0 and 1. 

# Plot the data to see the pattern:
cols = c('red','blue')
plot(X[,2]~X[,1],col = cols[Y+1], pch = 16)


# Specify some appropriate activation functions:
# Hidden layers
sig1 = function(z)
{
  tanh(z) # trust me, bro. (See Lecture 7)
}

# Output layer
sig2 = function(z)
{
  1/(1+exp(-z))  # This one is determined by the data (response). 
}

# Write a function that evaluates a neural network. 
# Takes arguments:
# X     - Design matrix
# Y     - Response matrix
# theta - A vector of parameters
# m     - Number of nodes on the hidden layer
# nu    - A regularisation parameter 

neural_net = function(X,Y,theta,m,nu)
{
   # Infer dimensions:
   N = dim(X)[1]
   p = dim(X)[2]
   q = dim(Y)[2]
   
   d = c(p,m,q)

   # Populate weight and bias matrices:
   index = 1:(d[1]*d[2]) # This is d_0 x d_1 
   W1    = matrix(theta[index],d[1],d[2])     # W1 is a d_0 x d_1 matrix
   index = max(index)+1:(d[2]*d[3]) # This is NEXT d_1 x d_2 elements
   W2    = matrix(theta[index],d[2],d[3])    # W1 is a d_1 x d_2 matrix
   index = max(index)+1:d[2]
   b1    = matrix(theta[index],d[2],1)  # b1 is a d_1 x 1 column vector
   index = max(index)+1:d[3]
   b2    = matrix(theta[index],d[3],1)  # b2 is a d_2 x 1 column vector
   
   ones = matrix(1,N,1) # This is for repeating bias vectors N times
   # Evaluate the updating equation in matrix form
  
     A0 = t(X)
     A1 = sig1(t(W1)%*%A0+b1%*%t(ones))
     A2 = sig2(t(W2)%*%A1+b2%*%t(ones))
  
     Y_hat =  t(A2)
     error = rep(0,N)
     wh0   = which(Y==0)
     wh1   = which(Y==1)
     error[wh0]  = log(1-Y_hat[wh0])
     error[wh1]  = log(Y_hat[wh1])
     
     #error = (Y*log(Y_hat)+(1-Y)*log(1-Y_hat)) # C_i values
     
   # Evaluate an appropriate objective function and return some predictions:
   E1 = -sum(error) # C = sum(C_i) 
   E2 = E1/N+nu/N*(sum(W1^2)+sum(W2^2))
   # Return a list of relevant objects:
   return(list(A2 = A2,A1 = A1, E1 = E1, E2 = E2))
}

# Create an objective function and evaluate it at a random coordinate in the 
# parameter space:
nu  = 0
m   = 10

p     = dim(X)[2] 
q     = dim(Y)[2]
npars = p*m+m*q+m+q
theta_rand = runif(npars,-1,1)


nu = 0
obj_pen = function(pars)
{
  res_model = neural_net(X,Y,pars,m,nu)
  return(res_model$E2) # E2 is the penalised objective!!!
}
obj_pen(theta_rand)

# Fit the neural network using a standard optimizer in R:
res_opt = nlm(obj_pen,theta_rand, iterlim = 1000)
res_opt

# Draw a response curve over the 2D input space to see
# what pattern the neural network predicts
M   = 5
x1  = seq(-3,+3,length = M)
x2  = seq(-3,+3,length = M)
xx1 = rep(x1,M)
xx2 = rep(x2,each = M)

abline(h= x2,v = x1, lwd = 0.3)

points(xx2~xx1,pch = 16, cex = 0.5)


XX = cbind(xx1,xx2)
YY = cbind(rep(1,M^2))
res_fitted = neural_net(XX,YY,res_opt$estimate,m,nu)

plot(XX[,2]~XX[,1],col = color.gradient(res_fitted$A2),pch = 16,cex = 0.5)
points(X[,2]~X[,1], col = cols[Y+1], pch = 16)


#===========================================================
# Validation Analysis
#===========================================================
set.seed(2)
N       = dim(X)[1]
set     = sample(1:N,0.8*N,replace = FALSE)
X_train = X[set,,drop = FALSE] #
Y_train = Y[set,,drop = FALSE]
X_val   = X[-set,,drop = FALSE]
Y_val   = Y[-set,,drop = FALSE]

obj_pen = function(pars)
{
  res_model = neural_net(X_train,Y_train,pars,m,nu)  # FIT ON THE TRAINING DATA!!!!
  return(res_model$E2) # E2 is the penalised objective!!!
}

n_nus = 15
nu_seq= exp(seq(-5,1,length = n_nus))
val_error = rep(0,n_nus)
for(i in 1:n_nus)
{
	nu = nu_seq[i] # Pick nu
	res_opt = nlm(obj_pen,theta_rand, iterlim = 1000) # Fit for that nu
	
	res_val      = neural_net(X_val,Y_val,res_opt$estimate,m,0)
  val_error[i] = res_val$E1
  
  print(paste0('Validation run:',i,' nu = ', nu))
  
  res_fitted = neural_net(XX,YY,res_opt$estimate,m,nu)

  plot(XX[,2]~XX[,1],col = color.gradient(res_fitted$A2),pch = 16,cex = 0.5, main = paste0('nu = ',round(nu,4)))
  points(X[,2]~X[,1], col = cols[Y+1], pch = 16)

  
}

plot(val_error~nu_seq, type = 'b',pch = 16)
plot(val_error~log(nu_seq), type = 'b',pch = 16)







