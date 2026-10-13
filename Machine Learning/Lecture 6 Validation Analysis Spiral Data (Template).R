rm(list = ls())
# A Function for creating continuous colour palette. I.e., colour gradient bins and
# associated colours:
library(colorspace) 
color.gradient = function(x, colors=c('magenta','white','lightblue'), colsteps=50) # quatitative color palette according to values of x
{
  colpal = colorRampPalette(colors)
  return( colpal(colsteps)[ findInterval(x, seq(min(x),max(x), length=colsteps)) ] )
}

# Read in the data and populate a design-matrix + response matrix:
dat = read.table('Spiral_Data_2022.txt',h = T)
X   = as.matrix(dat[,-3])
Y   = as.matrix((date[,3]==1)*1)

# Plot the data to see the pattern:
cols = c('red','blue')
plot(X[,2]~X[,1],col=cols[Y+1], pch = 16)


# Specify some appropriate activation functions:
# Hidden layers
sig1 = function(z)
{
  tanh(z) # trust for now (lecture 7 will justify this)
  
}



# Output layer
sig2 = function(z)
{
  1/(1+exp(-z)) # This one is determined by the repose data, we want to mimick repose
}

# Write a function that evaluates a neural network. 
# Takes arguments:
# X     - Design matrix
# Y     - Response matrix
# theta - A vector of parameters
# m     - Number of nodes on the hidden layer
# nu    - A regularisation parameter 

neural_net = function(X,Y,theta,m,nu) #theta are our initil parameter 
{
   # Infer dimensions:
   N = dim(X)[1] #data
   p = dim(X)[2]#reponse
   q = dim(Y)[2]#output dim, num categories
   
   d = c(p,m,q) # p=input dim, m-=hidden dim, q=output dim

   # Populate weight and bias matrices:
   index = 1:(d[1]*d[2]) # This is d_0 *d_1
   W1    = matrix(theta[index],d[1],d[2]) # d0*d1 matrix
   index = max(index)+1:(d[2]*d[3]) # next elements d1*d2
   W2    = matrix(theta[index],d[2],d[3]) #d2*d3 matrix
   index = max(index)+1:d[2]
   b1    = matrix(theta[index],d[2],1) #b1 id a d_1 *1 column vectore
   index = matrix(theta[index],d[3],1) 
   b2    = matrix(theta[index],d[3],1)#b1 id a d_1 *1 column vectore
   
   ones = matrix(1,N,1)
   # Evaluate the updating equation in matrix form
   #for(i in 1:N)
   #{
     A0 = t(X) # input matrix,features
     A1 = sig1(t(W1)%*%A0+b1%*%t(ones)) # hidden layer
     A2 = sig2(t(W2)%*%A1+b2%*%t(ones)) #output layer
   #}
     Y_hat = t(A2) # transpose output
     
     error=rep(0,N)
     wh0=which(Y==0)
     wh1=which(Y==1)
     
     error[wh0] = # -Y*log(Y_hat)+(1-Y)*log(1-Y_hat) # Cross Ent Objective Values at i
     error[wh1] = -Y*log(Y_hat)+(1-Y)*log(1-Y_hat) 
   # Evaluate an appropriate objective function and return some predictions:
   E1 = sum(error)
   E2 = E1+nu*(sum(W1^2))+sum(W2^2) # penalized/regularized objective we want to minimize
   # Return a list of relevant objects:
   return(list(A2 = A2,A1 = A1, E1 = E1, E2 = E2))
}

# Create an objective function and evaluate it at a random coordinate in the 
# parameter space:
nu  = 0.025 # reg parametr
m   = 10

p     = dim(X)[2]
q     = dim(Y)[2]
npars = p*m+m*q+m+q # num parameters across the network
theta_rand = runif(npars,-1,1) # random set of parameters


nu = 0
obj_pen = function(pars) # objective function to minimize
{
  res_model = nueral_net(X,Y,pars,m,nu)
  return(res_modeL$E1)
}
obj_pen(theta_rand) # cost value of the nn we have just coded- we want to minimize this

# Fit the neural network using a standard optimizer in R:
res_opt = nlm(obj_pen,theta_rand, iterlim = 1000) # use nlm function to optimize the above objective
res_opt

# Draw a response curve over the 2D input space to see
# what pattern the neural network predicts
M  = 100
x1 = seq(-3,3,length=M)
x2 = seq(-3,3,length=M)
xx1 = rep(x1,M) # repeat x1 M times
xx2 = rep(x1,each=M) # Fixes these

abline(h=x2,v=x1,lwd=0.3)

points(xx2,~xx1, cex=0.5, pch = 16)

XX = cbind(xx1,xx2)
YY = cbind(rep(1,M^2))
res_fitted = nueral_net(XX,YY,res_opt$estimate,m,nu)

plot(XX[,2]~XX[,1],col=color.gradient(res_fitted$A2),pch=16)
points(X[,2]~X[,1], col = cols[Y+1], pch = 16)


#===========================================================
# Validation Analysis
#===========================================================

set.seed(2)
N   = dim(X)[1]
set =sample(1:N,0.8*N,replace=F)
# Split up the data
X_train =X[set,,drop=F] # drop ensures that the dimensions remain the same
Y_train =Y[set,,drop=F]
X_val =X[-set,,drop=F]
Y_val =Y[-set,,drop=F]


obj_pen = function(pars) # objective function to minimize
{
  res_model = nueral_net(X_train,Y_train,pars,m,nu) # Fit on training data
  return(res_modeL$E1)
}


n_nu=10 # values for nus

nu_seq=exp(seq(-5,0,length=n_nus)) # using exp to start off small and ramp them up quickly
val_error=rep(0,nu_nus)

for(i in 1:n_nus){
  
  nu=nu_seq[i]
  res_val=nlm(obj_pen,theta_rand,iterlim=1000)
  res_opt=nueral_net(X_val,Y_val,res_val$estimate,m,nu)
  
  
}
