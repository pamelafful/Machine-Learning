 rm(list= ls(all = TRUE))
 # Fake some data:
 N  = 200

 x1  = runif(N,-1,1)
 x2  = runif(N,-1,1)
 X   = cbind(x1,x2)
 #y   = c(-1,1)[(0.2*x2+0.5*x1>=0.1)+1] # Encode responses as -1,1
 y   = c(-1,1)[(x2^2+x1^2>=0.7)+1]
 plot(x2~x1,pch = c(4,16)[(y+1)/2+1])

 library('quadprog')
 # Create data matrix:
 X  = cbind(x1,x2) 
 # Some space for the matrix in our QP problem.
 DD = matrix(0,N,N) # Dmat 

 # Define a polynomial kernel
 gm = 1 #coeff
 cf = 1 # constant
 dg = 2 # order of polynomial
 KK = function(X1,X2)
 {
     return((cf+gm*t(X1)%*%X2)^dg)
 }
 
 # Assign inner products to D
 for(i in 1:N)
 {
     for(j in 1:N)
     {
         DD[i,j] = y[i]*y[j]*KK(X[i,],X[j,]) # Second part of Lagragian  formulation yn*ym*zn*zm
     }
 }
 
 
 # D inverse is not a friendly chop... So we do what any respectable
 # applied mathematician does... We cheat (a little) by adding a small component
 # to the diagonal:
 eps = 5e-6
 DD  = DD +eps*diag(N)
   
 # Now for the constraints:
 # We have quality constraints and inequality constraints:
 # Eq  :    y'a = 0 t(A)b>bo
 # Ineq:    a >= 0 alpah=bvec
 Amat =  cbind(y,diag(N)) # y will be on first row of t(Amat) A, to ensure sum(apha*y)=0
 bvec =  matrix(0,N+1,1) #RHS
 d    =  matrix(1,N,1) # coeffiencents of all alpha is 1, d in formulation 

 # Solve for lagrange mults a:
 res = solve.QP(DD, d, Amat, bvec, meq=1, factorized=FALSE)
 a   = res$solution # alphas
 plot(a,type= 'h', main = expression(alpha[i]),xlab = expression(i), lwd = 2) # plot non-zero aphas/ support vectors
 
 # Recover the model:

 pad.a     = round(a,3)
 wh        = which.max(a) # find max alpha indicies (why?)
 
 # Find the ws:
 ww        = t(a*y)%*%X # work back the wieghts/ X be substituted by the kernel val
 
 T1 = rep(0,N) 
 for(i in 1:N)
 {
     T1[i] = sum(a*y*KK(X[i,],t(X))) # matriz to get predictions/plane
 } # 
 
 
 # How to find the intercept?
 #intercept = 1/y[wh]-ww%*%X[wh,]
 intercept = 1/y[wh]-sum(a*y*KK(X[wh,],t(X))) # get intercept back

 # Fitted values:
 #yhat      = sign(X%*%t(ww)+intercept[1])
 yhat = sign(T1+intercept[1])
 
 par(mfrow = c(2,2))
 plot(x2~x1,pch = c(1,16)[(y+1)/2+1])
 plot(x2~x1, pch = 16, col = c('grey','blue')[(yhat+1)/2+1])
 
 wh.text = which(pad.a!=0)
 points(x2~x1, pch = 1, col = c(NA,'red')[(pad.a>0)+1],cex=2)
 text(x2[wh.text]~x1[wh.text],labels = paste0('n = ',wh.text))
 
 
 # Now compare to a package solution:
 library('e1071')
 #?svm
 dat    = data.frame(y = as.factor(y), x2=x2, x1 =x1)
 model  = svm(y~., data = dat, cost = 1000, kernel = 'polynomial', coef0 = cf, degree = dg, gamma = gm,scale =F)
 
 # Our solution
 plot(x2~x1,pch = c(1,16)[(y+1)/2+1])
 points(x2~x1, pch = 1,col = c(NA,'red')[(pad.a>0)+1],cex=2)
 
 # Package solution
 points(model$SV[,1]~model$SV[,2],pch = '+', col = 'blue',cex = 2)
 plot(model$coefs~model$index,type = 'h',xlim = c(0,N))

