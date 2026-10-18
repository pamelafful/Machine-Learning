set.seed(2023)
f_x=function(x,sig=0){
  epsi=rnorm(1,0,sig)
  
  if(x<0){
    target=abs(x+1)-0.5
    return (target+epsi)
    
  }
  if(x>=0){
    target=abs(x-1)-0.5
    
    return(target+sig)
    
  }
}

# Simulte 5 time
n_sim=seq(1,5,1)
N=5

# plot(NULL, axes = FALSE, xlim=c(-2,2),ylim=c(-0.5,0.5),xlab='xi',ylab='f(xi)',main='Target Function Simulation')
# axis(1, pos = 0)
# axis(2, pos = 0)

par(mfrow = c(3,2 ))
colours=c('green','blue','lightblue','orange','red')
for(i in seq_along(n_sim)){
  xi=sort(runif(N,-2,2))
  y=sapply(xi,f_x)
  plot( y~xi,lwd=2,col=colours[1],type='l', main=paste('Target Function Simulation:',i))
  points(y~xi,lwd=2,col=colours[1],pch=1)
  
  # Plot Hypotheses
  h1=lm(y~xi)
  h1_predict=unname(predict(h1,data.frame(xi)))
  lines(h1_predict~xi,lwd=2,col=colours[3])
  points(h1_predict~xi,lwd=2,col=colours[3],pch=1)
  # 
  h2=lm(y~sin(pi*xi)+cos(pi*xi)-1) # no coeff regression
  h2_predict=unname(predict(h2,data.frame(xi)))
  lines(h2_predict~xi,lwd=1,col=colours[4])
  points(h2_predict~xi,lwd=1,col=colours[4],pch=1)
  # 
  
  legend("bottomright", 
         legend = c("f(x)", "H1", "H2"), # Text labels
         col = c(colours[1],colours[3],colours[4]),                           
         pch = 1, 
         cex=0.95,
         title = "Sample Plots",
         bty='n')  
  
}


# Funtion Arg
# f_tru: True DGP
# hyp: Hypothesis function
# sig: nr of std deviation from mean
# N: size of data
# nsim= nr of simulations/draws

oos_estimate=function(f_true,hyp,sig,N, nsim=1000){
  
  #RHS BV trade-off
  bias_store=c()
  prediction_sim_store=matrix(NA,nrow=N,ncol=nsim) # store predictions for each simulation in this matrix,col= prediction at simulation i
  pred_var_score=c()
  oos_err_store=c()
  for(i in 1:nsim){
    xi=sort(runif(N,-2,2)) # sample N samples
    y=sapply(xi,f_x) # apply the f_true
    
    h_predict=unname(predict(hyp,data.frame(xi))) # get prediction of the hypothesis on the sample
    
    # Calculate the bias at draw i
    bias_i= (y-h_predict)**2
    bias_store=c(bias_store,bias_i) # store the bias at i
    
    # Store predictions
    prediction_sim_store[,i]=h_predict # at each column, insert the prediction
    
    # get running mean of the prediction
    # use sapply for the first simulation draw
    mean_prediction=ifelse(i>1,apply(prediction_sim_store[,1:i],1,mean),sapply(prediction_sim_store[,1:i],mean)) 
    
    # how far each prediction is from the mean prediction
    pred_var_i=(h_predict-mean_prediction)**2
    pred_var_score=c(pred_var_score,pred_var_i)
    
    
    #LHS
    # create OOS datapoint
    xi_oos=sort(runif(N,-2,2)) # draw at random again
    y_oos=sapply(xi_oos,f_x)
    
    xi_oos_df=data.frame(xi=xi_oos)
    
    h_predict_oos=unname(predict(hyp,data=xi_oos_df))
    
    # calculate the error SE^2
    oos_err_i=(h_predict_oos-y_oos)**2
    oos_err_store=c(oos_err_store,oos_err_i)
    
    
    
  }
  
  #  get the averages
  model_bias2=round(mean(bias_store),4)
  model_var=round(mean(pred_var_score),4)
  both=round((model_bias2+model_var+sig),4)
  oos_err_est=round(mean(oos_err_store),4)
  estimate_diff=round(abs(both-oos_err_est),4)
  both_with_noise=both+sig
  
  
  return(list(model_bias2=model_bias2,model_var=model_var,both=both,oos_err_est=oos_err_est,estimate_diff=estimate_diff))
}


# Test the above
f_true=f_x
sig=0
N=5



h1_out=oos_estimate(f_true,h1,sig,N, nsim=1000)

h2_out=oos_estimate(f_true,h2,sig,N, nsim=1000)


comaprison_table=rbind(data.frame(h1_out),data.frame(h2_out)) %>% 
  mutate(model=c('h1(x)','h2(x)')) %>% 
  column_to_rownames(var = "model")

kable(comaprison_table, col.names = c('Bias','Variance','Both:Bias+Var','OOS Error','Estimate Differance'), align = "cc")

#BV Dynamics

set.seed(321)

sigmas=seq(0,0.5,0.1)


f_true=f_x
N=5

# apply function for all sigmas
h1_sigma_function_values=sapply(sigmas,oos_estimate,f_true=f_x,hyp=h1,N=5,nsim=1000)
h2_sigma_function_values=sapply(sigmas,oos_estimate,f_true=f_x,hyp=h2,N=5,nsim=1000)

# create bias plot
plot(NULL,xlim=c(0,0.5),ylim=c(0,0.095),xlab='sigma',ylab='bias2',main=' Model Bias vs Noise Parameter (sigma)',)
lines(sigmas,unlist(h1_sigma_function_values[1,]),col='blue',lwd=1.5)
points(sigmas,unlist(h1_sigma_function_values[1,]),col='blue',pch = 1)
lines(sigmas,unlist(h2_sigma_function_values[1,]),col='red',lwd=1.5)
points(sigmas,unlist(h2_sigma_function_values[1,]),col='red',pch = 1)


h1_oos_bias_mean=round(mean(unlist(h1_sigma_function_values[1,])),4)
h2_oos_bias_mean=round(mean(unlist(h2_sigma_function_values[1,])),4)


text3=paste('h1(x) Avg Bias:',h1_oos_bias_mean)
text4=paste('h2(x) Avg Bias:',h2_oos_bias_mean)

plot_text=paste(text3,text4,sep = "\n")

text(x =0.08, y = 0.04, labels =plot_text ,
     srt = 0, adj = c(0.5, -0.5), col = "blue",cex=1)


legend(x=0.4,y=0.04, 
       legend = c("h1(x)", "h2(x)"), # Text labels
       col = c('blue','red'),                           
       pch = 1, 
       cex=0.85,
       title = "Hypothesis",
       bty='n')  


plot(NULL,xlim=c(0,0.5),ylim=c(0,0.075),xlab='sigma',ylab='Variance',main=' Model Variance vs Noise Parameter (sigma)',)
lines(sigmas,unlist(h1_sigma_function_values[2,]),col='blue',lwd=1.5)
points(sigmas,unlist(h1_sigma_function_values[2,]),col='blue',pch = 1)
lines(sigmas,unlist(h2_sigma_function_values[2,]),col='red',lwd=1.5)
points(sigmas,unlist(h2_sigma_function_values[2,]),col='red',pch = 1)



h1_oos_variance_mean=round(mean(unlist(h1_sigma_function_values[2,])),4)
h2_oos_variance_mean=round(mean(unlist(h2_sigma_function_values[2,])),4)


text3=paste('h1(x) Avg Var:',h1_oos_variance_mean)
text4=paste('h2(x) Avg Var:',h2_oos_variance_mean)

plot_text=paste(text3,text4,sep = "\n")

text(x =0.08, y = 0.04, labels =plot_text ,
     srt = 0, adj = c(0.5, -0.5), col = "blue",cex=1)



legend(x=0.4,y=0.04,  
       legend = c("h1(x)", "h2(x)"), # Text labels
       col = c('blue','red'),                           
       pch = 1, 
       cex=0.85,
       title = "Hypothesis",
       bty='n')  


plot(NULL,xlim=c(0,0.5),ylim=c(0.07,0.138),xlab='sigma',ylab='E[E_OOS_Err]',main=' Model Expected OOS Error vs Noise Parameter (sigma)',)
lines(sigmas,unlist(h1_sigma_function_values[4,]),col='blue',lwd=1.5)
points(sigmas,unlist(h1_sigma_function_values[4,]),col='blue',pch = 1)
lines(sigmas,unlist(h2_sigma_function_values[4,]),col='red',lwd=1.5)
points(sigmas,unlist(h2_sigma_function_values[4,]),col='red',pch = 1)

legend(x=0.4,y=0.11, 
       legend = c("h1(x)", "h2(x)"), # Text labels
       col = c('blue','red'),                           
       pch = 1, 
       cex=0.75,
       title = "Hypothesis",
       bty='n')  




h1_oos_estimate_mean=round(mean(unlist(h1_sigma_function_values[4,])),4)
h2_oos_estimate_mean=round(mean(unlist(h2_sigma_function_values[4,])),4)




text3=paste('h1(x) Avg OOS Error Est:',h1_oos_estimate_mean)
text4=paste('h2(x) Avg OOS Error Est:',h2_oos_estimate_mean)

plot_text=paste(text3,text4,sep = "\n")

text(x =0.099, y = 0.105, labels =plot_text ,
     srt = 0, adj = c(0.5, -0.5), col = "blue",cex=1)


# Averages across

h1_both=h1_oos_variance_mean+h1_oos_bias_mean
h2_both=h2_oos_variance_mean+h2_oos_bias_mean
h1_both;h2_both
