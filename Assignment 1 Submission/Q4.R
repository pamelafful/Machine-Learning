dat = read.table('Collider_Data_2022.txt', h = TRUE,stringsAsFactors =T) 

# Plot

plot(NULL, xlim=c(min(dat['X1']),max(dat['X1'])),ylim=c(min(dat['X2']),max(dat['X2'])), xlab='X1',ylab='X2',main='Classification Scatter Plot')

points(dat[dat$Yi1==1,]$X1,dat[dat$Yi1==1,]$X2,col='blue',pch=10)

points(dat[dat$Yi2==1,]$X1,dat[dat$Yi2==1,]$X2,col='red',pch=10)

points(dat[dat$Yi3==1,]$X1,dat[dat$Yi3==1,]$X2,col='purple',pch=10)

legend(x='topleft',
       legend = c("Yi1", 'Yi2','Yi3'),
       col=c('blue','red','purple'),
       pch=c(10,10,10),
       title = "Groups",
       cex=0.85,
       bty = "o")


#soft max

soft_max = function(z) {  # assume z is a matrix
  
  z_shifted = sweep(z, 2, apply(z, 2, max), FUN = "-") # subtract maximum from each column to avoid default R behaviour when adding numbers of varying magnitudes
  
  
  num = exp(z_shifted)
  
  # Sum exponentials within each column
  sums = colSums(num)
  
  # Normalize each column
  result = sweep(num, 2, sums, FUN = "/")
  
  return(result)
}


#NN

X_feat=cbind(dat$X1,dat$X2)
Y=cbind(dat$Yi1,dat$Yi2,dat$Yi3)

choose_label=function(vector_prob){
  
  labels_vec=ifelse(vector_prob==max(vector_prob),1,0)
  return(labels_vec)
  
}


feedforward_nn=function(X_input,y_input,theta,m,nu,activation_fun_2=soft_max,activation_fun_1=tanh){
  
  # theta is the vector of initial weights
  #it is of length p*m+m*q+m+q (W1+W2+b1+b2)
  # in this code q=d2
  
  y=as.matrix(y_input)
  
  # Initialize relevant layer matrices
  
  
  N=dim(X_input)[1] # nr obs
  p=dim(X_input)[2] # nr of cols=nr of nodes in input layer for obs i
  #1. nr of nodes in input layer=2 features across N observations
  
  # Going into layer 1
  d0=p # nodes in input
  d1=m # nodes in hidden
  d2=dim(y)[2]
  
  X=t(as.matrix(X_input)) # obs as cols
  
  index_chunk1=p*m #W1
  index_chunk2=index_chunk1+1 #W2
  index_chunk3=((index_chunk2+(m*d2))-1)+1 #B1
  index_chunk4=((index_chunk3+m)-1)+1 #B2
  
  
  # weight matrix from input layer to hidden layer
  W1=matrix(theta[1:(index_chunk1)],nrow=d0,ncol=d1) 
  
  #weight matrix from hidden layer to output layer
  W2=matrix(theta[(index_chunk2):((index_chunk2+(m*d2))-1)],nrow=d1,ncol=d2) # 
  
  # bias matrix going into hidden layer
  b1=matrix(theta[(index_chunk3):((index_chunk3+m)-1)],nrow=m,ncol=1) 
  
  b1_ones_matrix=matrix(1,nrow=1,ncol=N)
  b1_mat=b1%*%b1_ones_matrix
  
  
  # bias matrix going  into output/layer 2
  b2=matrix(theta[index_chunk4:((index_chunk4+d2)-1)],nrow=d2,ncol=1)
  b2_ones_matrix=matrix(1,nrow=1,ncol=N) # Repeat this N times
  b2_mat=b2%*%b2_ones_matrix
  
  # Do one forward calc and get prediction
  
  # first layer computation
  first_layer_comp=t(W1)%*%X+b1_mat
  
  # apply activation
  first_layer_comp=activation_fun_1(first_layer_comp)
  
  # second layer computation
  
  second_layer_comp=activation_fun_2(t(W2)%*%first_layer_comp+b2_mat)
  second_layer_comp_no_activation=t(W2)%*%first_layer_comp+b2_mat
  
  
  #predictions across N observations 
  # prediction = max probability
  y_preds=second_layer_comp # 3 by N
  y_pred_label=apply(y_preds,1,choose_label) # one hot encoding
  
  
  # control for very small values in the matrix
  log_compute= function(mat, threshold = 1e-6) {
    ifelse(mat < threshold, log1p(mat), log(mat))
  }
  
  
  CE_mat=t(y)*log_compute(y_preds)
  CE_error=-mean(apply(CE_mat,1,sum)) # normal non_penalized oobj
  
  CE_error_L2=CE_error+(nu/N)*(sum(W1^2)+sum(W2^2)) # regularized obj
  
  
  
  return(list(first_layer_comp = first_layer_comp,second_layer_comp = second_layer_comp, CE_error = CE_error, CE_error_L2 = CE_error_L2,second_layer_comp_no_activation=second_layer_comp_no_activation))
  
  
}


# Validation


set.seed(2022)

y_input=Y
X_input=X_feat
# 
m   = 10 # number nodes in hidden layer
p     = dim(X_input)[2]
q     = dim(as.matrix(Y))[2]

npars = p*m+m*q+m+q
theta_rand = runif(npars,-1,1)
N       = dim(X_feat)[1]
set     = sample(1:N,0.5*N,replace = FALSE) #50 50 split
X_train = X_feat[set,,drop = FALSE] #
Y_train = Y[set,,drop = FALSE]
X_val   = X_feat[-set,,drop = FALSE]
Y_val   = Y[-set,,drop = FALSE]

obj_pen = function(pars)
{
  res_model = feedforward_nn(X_train,Y_train,pars,m,nu,activation_fun_2=soft_max,activation_fun_1=tanh)  
  return(res_model$CE_error_L2) # Reg Obj
}

n_nus = 15
nu_seq= exp(seq(-2,2,length = n_nus)) # sequence of nus to validate against
val_error = rep(0,n_nus)
for(i in 1:n_nus)
{
  nu = nu_seq[i] # Pick nu
  res_opt = nlm(obj_pen,theta_rand, iterlim = 1000) # Fit for that nu
  
  #Evaluate on Validation Set with optimal parameters and particular nu
  
  res_val= feedforward_nn(X_val,Y_val,res_opt$estimate,m,nu,activation_fun_2=soft_max,activation_fun_1=tanh)
  
  # Get reg validation Error
  val_error[i] = res_val$CE_error_L2
  
  #print(paste0('Validation run:',i,' nu = ', nu))
}


#Plot

min_val_error=min(val_error)



best_nu_index=which(val_error==min(val_error))
best_nu=nu_seq[best_nu_index];


# Compare the error of the nu the right of the best

next_val_error=val_error[(best_nu_index+1):length(val_error)]

next_seq=nu_seq[(best_nu_index)+1:length(nu_seq)]

next_best_nu=next_seq[which.min(next_val_error)] ;

next_min_val_error=min(next_val_error)




plot(val_error~log(nu_seq), type = 'b',pch = 16,xlim=c(min(log(nu_seq)),max(log(nu_seq))+1), main='Validation Error vs Regularization Parameter (nu)', xlab='lambda (log)',ylab='validation error')

abline(v =log(best_nu), col = "purple", lty = "dashed", lwd = 2)

abline(v =log(next_best_nu), col = "blue", lty = "dashed", lwd = 2)

text(x = log(best_nu)+0.33, y = 30, labels = paste("best nu:",round(best_nu,3), 'Val Error:',round(min_val_error,2),sep = "\n"),
     srt = 0, adj = c(0.5, -0.5), col = "purple",cex=0.68)


text(x = log(next_best_nu)+0.55, y = 30, labels = paste("next best nu:",round(next_best_nu,3), 'Val Error:',round(next_min_val_error,2),sep = "\n"),
     srt = 0, adj = c(0.5, -0.5), col = "blue",cex=0.68)
# 
# 


# Reps Curve

library(colorspace)

# Fit probability color mapping for each class onto the feature space



# Helper functions

# Assign color pallet based on indexing
ass_pallet=function(ind_val,col_steps=50){
  if(ind_val==1){
    
    # Class 1, high to low prob color mapping
    pallet=colorRampPalette(c("lightblue", 'white',"darkblue")) 
    
  }
  
  if(ind_val==2){
    
    # Class 2
    pallet=colorRampPalette(c("mistyrose",'white', "darkred"))
    
  }
  if(ind_val==3){
    
    # Class 3
    pallet=colorRampPalette(c("mediumpurple1",'white', "purple4"))
    
  }
  
  return(pallet(col_steps))
}


# Assign a color to prediction based on its class and strength of probability value

col_pred=function(x,col_steps=50){
  # x=prediction matrix from ff nn
  
  max_ind=apply(x,2,which.max) #get index of max probability
  max_prob=apply(x,2,max) # get value of max probability
  col_range=seq(0,1,length.out = col_steps) 
  col_pallet=sapply(max_ind,ass_pallet) #prob colors for each class, matrix
  color_indices= findInterval(max_prob, col_range) # for each obs, find its indx in a range of 0 to 1
  
  col_pallet_N=ncol(col_pallet)
  selected_cols= col_pallet[cbind(color_indices, seq_len(col_pallet_N))] # get appropriate index in the pallette assignment matrix
  
  return(selected_cols)
}

pal1 = ass_pallet(1)
pal2 = ass_pallet(2)
pal3 = ass_pallet(3)

obj_pen = function(pars)
{
  res_model =feedforward_nn(X_input,y_input,pars,m,next_best_nu,activation_fun_2=soft_max,activation_fun_1=tanh)
  return(res_model$CE_error_L2) # Reg Obj
}

# Optimize weights using next best nu
res_opt_best_nu = nlm(obj_pen,theta_rand, iterlim = 1000)



M   = 160
x1  = seq(-5,+5,length = M)
x2  = seq(-5,+5,length = M)
xx1 = rep(x1,M)
xx2 = rep(x2,each = M)

layout(matrix(1:4, nrow = 1),
       widths = c(5, 0.7, 0.7, 0.7))

par(mar = c(4, 4, 3, 1))

plot(NULL, xlim=c(-5,5),ylim=c(-5,5),xlab='x1',ylab='x2', main='Feed Forward NN Response Curve')

abline(h= x2,v = x1, lwd = 0.3,col='black')




XX = cbind(xx1,xx2)
YY1 = rep(1,M^2)
YY2=rep(1,M^2)
YY3=rep(1,M^2)
YY=cbind(YY1,YY2,YY3)


# Again use next best nu
res_fitted=feedforward_nn(XX,YY,res_opt_best_nu$estimate,m,next_best_nu,activation_fun_2=soft_max,activation_fun_1=tanh) 

points(xx2~xx1,pch = 16, col=col_pred(res_fitted$second_layer_comp) ,cex = 2)

points(dat[dat$Yi1==1,]$X1,dat[dat$Yi1==1,]$X2,col='blue',pch=10)

points(dat[dat$Yi2==1,]$X1,dat[dat$Yi2==1,]$X2,col='red',pch=10)

points(dat[dat$Yi3==1,]$X1,dat[dat$Yi3==1,]$X2,col='purple',pch=10)

add_colour_bar = function(pal, class_name) {
  
  par(mar = c(4, 0.5, 3, 2))
  
  probs <- seq(0, 1, length.out = length(pal))
  
  # Draw continuous-looking gradient
  z = matrix(rep(probs, each = 2), nrow = 2)
  
  image(x = c(0, 1),
        y = probs,
        z = z,
        col = pal,
        zlim = c(0, 1),
        axes = FALSE,
        xlab = "",
        ylab = "",
        main = class_name)
  
  axis(4, at = seq(0, 1, by = 0.2),
       las = 1, cex.axis = 0.75)
  
  box()
}

#colour mapping
add_colour_bar(pal1, "Class 1")
add_colour_bar(pal2, "Class 2")
add_colour_bar(pal3, "Class 3")

