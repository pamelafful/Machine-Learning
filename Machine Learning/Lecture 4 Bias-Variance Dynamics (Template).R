# Example 2.8 in slides

rm(list=ls()) # listing all of the objects / flushes the workspace
# Set up a target function/pattern:
k=1
f = function(x){sin(k*pi*x)}

# Write a function that generates a data set of size N 
# with pattern f. Return a list of predictors (x) and responses (y).
dgp = function(N,f,sig)
{ x=runif(N,-1,1)
e=rnorm(N,0,sig)
y=f(x)+e # add noise to the system
return (list(y=y,x=x))

}
res=dgp(2,f,0);res

1
# Plot a single realization of the data and overlay the target function:
plot(,xlim = c(-1,1), ylim = c(-1,1))
xx_lat =seq(-1,1,1/100) 
lines(f(xx_lat)~xx_lat)


# Function for bias-variance

bias_var=function(N,order,sig,M=5000,dx=1/100){
  xx_lat=seq(-1,1,dx) #dx steps of the lattice
  Nx=length(xx_lat)
  
  g_bar=xx_lat*0 # empty vector of zeros, Average model
  G_D=matrix(0,M,Nx) # Particular models fitted at each of M sims
  
  plot(f(xx_lat)~xx_lat,type='l',ylim=c(-3,3))
  
  test_error_D=0
  test_error=0
  
  for(i in 1:M){ # number of M data simulations
    
    dat_in=dgp(N,f,sig)
    
    if(order==0){
      res_model=lm(y~1,data=dat_in)
    }
    
    if(order>0){
      res_model=lm(y~poly(x,order),data=dat_in)
    }
    
     # Fit H0
    gd= predict(res_model,data.frame(x=xx_lat)) #get predictions of the above
    lines(gd~xx_lat,col=4,lwd=0.2) # for each model realization
    G_D[i,]=gd # store predictions for this realization of the data
    g_bar=g_bar+gd #Running mean calculation
    
    # Sample estimates
    dat_out=dgp(N,f,sig)
    pred_out=predict(res_model,data.frame(dat_out$x))
    test_error_D=mean((dat_out$y-pred_out)**2)
    test_error=test_error+test_error_D # Running mean
    
    
  }
  
  # Calculate the bias using left point Riemann 
  phi=1/2 # from distribution of X! U(-1,1)
 g_bar=g_bar/M # average prediction at xx_lattice 
 lines(g_bar~xx_lat,col=4,lwd=0.2,col='green') # average model
 
 bias2=sum((gbar-f(xx_lat))**2)[-Nx]*phi*dx # Look up Riemann Sum
 
 
 # Model Variance (Two Integration steps)
 ones=matrix(1,M,1)
 
 var_at_x= colSums((G_d-ones%*%g_bar)**2)/M # repaeting gbar along the rows; Monte Carlo step
 lines(c(g_bar+1*sqrt(var_at_x))~xx_lat,col='red',lwd=2)
 lines(c(g_bar-1*sqrt(var_at_x))~xx_lat,col='red',lwd=2)
 
 var=sum( var_at_x[-Nx]*phi*dx) # Rieman again
 
 both=bias2+var # OOS Error
 
 ret=list(bias=bias2,var=var,both=both,test_error=test_error)
 
 return(ret)
 
  
 
}

# Run the above
bias_var(5,0,0,M=1000,dx=1/100) # look through arguments

# Even id underlying pattern is complex, we are better of using a simpler model- du to the noise floor (or something)
# Mean of our model's OOS error is an estimate of our BV/Expected value of your OOS error
  


# Calculate the (squared) bias for an appropriate H. 
# Use M = 1000 draws of the data and a step-size of dx = 1/100 for 
# the left-(end)point approximation to the relevant integral.


# Calculate the variance for an appropriate H. 
# This integration involves two steps.

# Calculate the expected Out-of-Sample error (over all D).
# Use sample estimates for M simulations. 
# Call it Test_Error

# Now, wrap all of this in a function, bias_var, that returns  
# bias2, var, bias2 + var, and Test_Error, for particular N,
# order, M, and dx and answer questions for the configurations
# in the appended Example 2.8.
