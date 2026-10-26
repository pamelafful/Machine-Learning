rm(list = ls())
#========================================================================
#                 1. Set up a game universe (Encoding)
#========================================================================

set.seed(2023)
J       = 30
o1      = runif(J,-1,+1)    # x-coordinates of objects
o2      = runif(J,-0.5,+1)    # y-coordinates of objects
objects = cbind(o1,o2)
xt      = cbind(runif(1,-0.8,0.8),-1)     # 1 random starting position
delt    = 0.05             # How big are the steps you can take?
rad     = 0.05             # How close to an object before you crash?

# Just a function to plot the game:
plot_game = function(xt,objects,rad,cols = 1)
{
  plot(objects[,2]~objects[,1], type ='n', col = 'darkgreen', ylim = c(-1,1), xlim = c(-1,1),xlab ='South',ylab = 'West')
  mtext('North',3);mtext('East',4)
  points(xt[,2]~xt[,1],pch = 9,col = cols)
  symbols(objects[,1],objects[,2],circles = rep(rad,dim(objects)[1]),ylim = c(-1,1),xlim = c(-1,1),inches = FALSE,add = TRUE,bg = 'darkgreen')
  abline(v = c(-1,1), lty = 3,col = 'black',lwd = 2)
  abline(h = c(-1,1), lty = 3,col = c('black','black'), lwd = 2)
}

plot_game(xt,objects,rad,'black')

#========================================================================
#                  2. Evaluating the game state
#========================================================================

game_status=function(xt,objects, rad)
{
  status    = 0
  J         = dim(objects)[1]
  ones      = matrix(1,J,1)
  dists     = sqrt(rowSums((objects-ones%*%xt)^2))
  min_dist  = min(dists)
  status    = +1*(xt[2] >= +1)-1*((xt[2] < -1)|(xt[1]< -1)|(xt[1]> +1)|(min_dist<rad))
  ret = list(status = status, closest = objects[which.min(dists),,drop = FALSE])
  return(ret)
}

# Just randomly move piece for now:
#control = function(xt,theta)
#{
#  return(cbind(runif(1,-1,1),1))
#}

play = function(x0,delt,objects,rad,theta,plt = FALSE,trace = FALSE)
{
  k           =   1 # Count how many steps
  xt          =   x0 # Set the initial coordinate(s) for the drone.
  trajectory  = matrix(NA,2,100)
  trajectory[,1] = xt

  # Check the game status:
  res_status  = game_status(xt,objects,rad)

  if(plt)
  {
  	plot_game(xt,objects,rad,c('red','black','green')[res_status$status+2])
  }
  while((res_status$status==0)&(k<100))
  {
    k = k+1
    
    # Now, let the model update the position of the pieces:
    ct = control(xt,theta)
    xt = xt + ct*delt

    # Check the game status after the positions are updates:
    res_status  = game_status(xt,objects,rad)

    if(trace)
    {
    	  trajectory[,k] = xt
    }
    if(plt)
    {
      plot_game(xt,objects,rad,c('red','black','green')[res_status$status+2])
    	title(main = paste0('Step = ',k))
    }
  }
  if(trace)
  {
  	lines(trajectory[2,]~trajectory[1,],col = 'lightgrey')
  }
  return(list(k = k, status = res_status$status,xt= xt,trajectory = trajectory))
}

#play(xt,delt,objects,rad,theta,plt = TRUE,trace = TRUE)


#========================================================================
#                  3. Control (Giving a Model Agency.)
#========================================================================

model = function(xt,theta,nodes)
{
  # Infer dimensions:
   N = dim(xt)[1]
   p = dim(xt)[2]
   q = 2
   dims = c(p,nodes,q)
  
   # Populate weight and bias matrices:
   index = 1:(dims[1]*dims[2])
   W1    = matrix(theta[index],dims[1],dims[2])
   index = max(index)+1:(dims[2]*dims[3])
   W2    = matrix(theta[index],dims[2],dims[3])
   index = max(index)+1:(dims[3]*dims[4])
   W3    = matrix(theta[index],dims[3],dims[4])

   index = max(index)+1:(dims[2])
   b1    = matrix(theta[index],dims[2],1)
   index = max(index)+1:(dims[3])
   b2    = matrix(theta[index],dims[3],1)
   index = max(index)+1:(dims[4])
   b3    = matrix(theta[index],dims[4],1)

   ones    = matrix(1,N,1)
   a0      = t(xt)
   
   # Evaluate the updating equation in matrix form
   a1 = tanh(t(W1)%*%a0+b1%*%t(ones))
   a2 = tanh(t(W2)%*%a1+b2%*%t(ones))
   a3 = tanh(t(W3)%*%a2+b3%*%t(ones))

   # Return a list of relevant objects:
   return(list(a3 = t(a3)))
}

p     = 2
q     = 2
nodes = 2
npars = p*nodes+nodes*nodes+nodes*q+nodes+nodes+q
npars
theta_rand = runif(npars,-1,1)
res_net = model(xt[,,drop=FALSE],runif(npars,-1,1),rep(nodes,2))

control = function(xt,theta)
{
  res_model = model(xt[,,drop=FALSE],theta,rep(nodes,2))
  return(res_model$a3)
}


#========================================================================
#                  4. Objectives and Fitness
#========================================================================

T_games = 1

play_T_games = function(theta)
{
	score = 0
	for(i in 1:T_games)
	{
   xt       = cbind(runif(1,-0.8,0.8),-1)
   res_play = play(xt,delt,objects,rad,theta,plt = FALSE,trace = FALSE)
   score    = score+(res_play$status==+1)*1
	}
  return(score/T_games)
}
play_T_games(runif(npars,-1,1))


#========================================================================
#                  5. Evolutionary Learning
#========================================================================

library('GA')
GA <- ga(type = 'real-valued',
	       fitness = play_T_games,
	       lower = rep(-10,npars),
	       upper = rep(+10,npars),
	       popSize = 200, # Here's the issue
	       maxiter = 100,
	       keepBest = TRUE) # Set popSize = 100 to get the bad run from class. 
plot(GA)

theta_hat = GA@solution[1,]
theta_hat
 
 
# Does it work?


win = c()
for(i in 1:10)
{
  xt_try    = cbind(runif(1,-0.8,0.8),-1)
  plot_game(xt_try,objects,rad,'black')
  res_final = play(xt_try,delt,objects,rad,theta_hat,plt = FALSE,trace = TRUE)
  lines(res_final$trajectory[2,]~res_final$trajectory[1,])
  win[i] = (res_final$status==1)
}
mean(win)

