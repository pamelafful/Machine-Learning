lq=function(x,q){
  coeff=2**q
  summands_store=c()
  for(k in 0:q) { # These are the ks we are summing over
    
    summands_k=(x**k)*(choose(q,k)*choose(((q+k-1)/2),q)  )
    summands_store=c(summands_store,summands_k)
  }
  
  f_x=coeff*sum(summands_store)
  
  return(f_x)
  
}


# plot the above
# USE sapply to vectorize the calculation across a supplied vector
x_list=list(x=seq(-1,1,0.1))

x=seq(-1,1,0.1)
Qs=seq(0,5,1)

q_cols=c('blue','red','magenta','lightblue','purple','green')



par(mfrow = c(2, 3),oma = c(0, 0, 3, 0))

for(i in 1:length(Qs)){ #c(min(sapply(x,lq,q=5)),max(sapply(x,lq,q=5)))
  q_i=Qs[i]
  plot(x,sapply(x,lq,q=q_i),col=q_cols[i],lwd=1.2,xlab='x',ylab='Lq(x)',main= paste('Lq(x) q=',q_i))
  
  
}

mtext("Legendre Polynomials of order q", side = 3, outer = TRUE, cex = 1.5, font = 2)

set.seed(30)

alpha_qs=runif(3,-1,1)
beta_qs=runif(3,-1,1)

target_alpha=function(x,alpha_q,Qf){
  store_q=c()
  for(q_i in 0:Qf){
    
    store_q=c(store_q,alpha_q*(x**q_i))
  }
  
  ret=sum(store_q)
  return(ret)
}


target_beta=function(x,beta_q,Qf){
  store_q=c()
  for(q_i in 0:Qf){
    Lq=lq(x,q_i)
    store_q=c(store_q,beta_q*Lq)
  }
  
  ret=sum(store_q)
  return(ret)
}


# Plot the above target_alpha

plot(NULL,xlim=c(-1,1),ylim=c(-5,5),xlab='x',ylab='f(x)',main='Equation 2 Random Target Functions')

q_cols=c('blue','red','magenta','lightblue','purple','green')

for(i in 1:length(alpha_qs)){
  alpha=alpha_qs[i]
  
  lines(x,sapply(x,target_alpha,alpha_q=alpha,Qf=5),col=q_cols[i],lwd=1.2)
}


legend("bottomleft", 
       legend = c( paste('alpha:',round(alpha_qs[1],2)),paste('alpha:',round(alpha_qs[2],2)),paste('alpha:',round(alpha_qs[3],2))), # Text labels
       col = q_cols[1:3],                           
       pch = 1, 
       cex=0.85,
       title = "Equation 2",
       bty='n') 

# Plot the above target_beta


plot(NULL,xlim=c(-1,1),ylim=c(-2,2),xlab='x',ylab='f(x)',main='Equation 3 Random Target Functions')

for(i in 1:length(beta_qs)){
  beta_i=beta_qs[i]
  lines(x,sapply(x,target_beta,beta_q=beta_i,Qf=5),col=q_cols[i],lwd=1.2)
}


legend("topleft", 
       legend = c( paste('beta:',round(beta_qs[1],2)),paste('beta:',round(beta_qs[2],2)),paste('beta:',round(beta_qs[3],2))), # Text labels
       col = q_cols[1:3],                           
       pch = 1, 
       cex=0.85,
       title = "Equation 3",
       bty='n')



