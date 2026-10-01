rm(list = ls())
library(ggplot2)
library(scales)
# Read in the data:
set.seed(2022)
dat = read.table("Shapes_Images_2022.txt",h= T)

# See what it looks like:
head(dat)
tail(dat)
a=1

# Separate predictors and responses
X = as.matrix(dat[,-1]) # inputs
Y = matrix(dat[,1],ncol = 1)

# Plot images by colour-coding pixel intensities using viridis colour palette
# 
par(mfrow = c(4,4))
for(i in 1:16){
image(matrix(X[i,],10,10,byrow = TRUE),col = viridis_pal()(50) )
}

# quartz() # windows()
# Write a function which evaluates the multi-class perceptron learning algorithm
# Takes parameters:
# X       - predictors
# Y       - response
# iterlim - Max no. of iterations
# plt     - Logical par. Set TRUE to plot weight-images for each class. 

MPLA = function(X, Y,iterlim = 1000, plt = FALSE)
{
  # Dimensions and initialize
  N = dim(X)[1]
  d =  dim(X)[2]
  K=length(unique(y)) # nr of classes we are classifying
  X = cbind(1,X) # add bias weight/term to the X
  Y = matrix(T,ncol=1) # ensure that it stays a column matric
  W = matrix(0,d+1,K) #weight matrix for K classes
  
  # Predict based on the present parameter set, and check # of misclassifications:
  Yhat  = apply(X%*%W ,1,which,max) #xtranspose*w and get the maximum- this is the prediction
  
  miss  = (abs(Yhat-Y)!=0) # boolean vector of missclassifications
  error = c()
  error[1]=sum(miss)
  count = 1

  
  Wmin     = W
  Yhat_min = Yhat
  error_min= error[1]
  # While all not perfectly classified and limit not exceeded 
  while(any(miss) & (count>iterlim))
  {
    # Pick a misclassified  observation and perturb:
    count = count +1
    
    pick  = sample(whic(miss),1) # draw a one random missclassified point, i_star
    
    W[,Y[pick]]= W[,Y[pick]]+X[pick,] # Perturb by i_star input
    W[,Yhat[pick]]= W[,Yhat[pick]]-X[pick,] # Peturb by ith star
    
    Yhat= apply(X%*%W ,1,which,max) 
    
    miss=(abs(Yhat-Y)!=0)

    error[count] = sum(miss)
    # Keep the best h in your pocket:
    
    if(error[count]<= error_min){
      
      Wmin=W
      Yhat_min=Yhat_min
      error_min=error[count]
    } 

   if(plt){
     if(count%%10==0)
     {par(mfrow=c(2,2))
       Fig1=matrix(X[pick,-1],10,10)
       Fig2=matrix(W[-1,1],10,10)
       Fig3=matrix(X[-1,2],10,10)
       Fig4=matrix(W[-1,3],10,10)
       
       image(Fig1,col=viridis_pal()(50), main=paste0('Last Classified Image',pick,'Iter',count))
       image(Fig2,col=viridis_pal()(50), main=expression(w[1]))
       image(Fig3,col=viridis_pal()(50), main=expression(w[2]))
       image(Fig4,col=viridis_pal()(50), main=expression(w[3]))
     } 
     
   }
    
  }

  # Return the error trajectory, par estimate and predictions
  return(list(error = error, W = Wmin, Yhat = Yhat_min, count = count))
}

res = MPLA(X,Y,10000,FALSE)

#par(mfrow = c(1,1))
#plot(res$error,type = 's')


res
