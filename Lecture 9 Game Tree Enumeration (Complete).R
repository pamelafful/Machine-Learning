rm(list = ls())

#===============================================================================
# Populate the state matrix
#===============================================================================
S = matrix(0,9,8)
# Row Sums:
S[1:3,1] = 1
S[4:6,2] = 1
S[7:9,3] = 1
# Col Sums:
S[c(1,4,7),4] = 1
S[c(2,5,8),5] = 1
S[c(3,6,9),6] = 1
# Diag Sums
S[c(1,5,9),7] = 1
S[c(3,5,7),8] = 1
S

#===============================================================================
# Evaluate the Game State rho(m,S)
#===============================================================================
rho = function(m,S)
{
  m = matrix(m,ncol = 1)
  player1  = any(t(m)%*%S == -3)
  player2  = any(t(m)%*%S == +3)
  m_1      = (m == -1)
  m_2      = (m == +1)
  tie      = sum((t(m_1)%*%S>0)*(t(m_2)%*%S>0)) == 8
  winner   = c(-1,0,+1)[c(player1,tie,player2)]
  terminal = player1|tie|player2
  ret      = list(terminal = terminal,winner = winner)
  return(ret)
}

m = as.matrix(c(1,1,-1,0,0,0,0,0,-1))
matrix(m,3,3,byrow = TRUE)
rho(m,S)

m = as.matrix(c(1,1,-1,0,0,-1,0,0,-1)) 
matrix(m,3,3,byrow = TRUE)
rho(m,S)

#===============================================================================
# game_tree(m,k) returns a vector of -1, 0, and +1s delineating all terminal 
# game states.
#===============================================================================
game_tree= function(m,k)
{
	g = c()
	game_state = rho(m,S)
	if(game_state$terminal)
	{
		 g = c(g,game_state$winner)
		 return(g)
	}else{
		index = which(m == 0)
		for(i in 1:length(index))
		{
			x           = m
			x[index[i]] = k
			g           = c(g,game_tree(x,-1*k))
		}
		return(g)
	}
}
#m = as.matrix(c(0,0,0,0,0,0,0,0,0))
m = as.matrix(c(1,1,0,0,0,0,0,0,-1))

res   = game_tree(m,-1)
n_g   = 	length(res)
Xwins = sum(res==-1)
Draws = sum(res== 0)
Owins = sum(res==+1)
c(n_g,Xwins,Draws,Owins)

#===============================================================================
# Write a function that draws the game tree.
#===============================================================================
draw_game_tree= function(m,k,x0 =0,d = 0,lims =c(-1,+1))
{
	if(d == 0) 
	{
		 plot(1,1,type ='n',xlim = lims, ylim = c(0,sum(m==0)+2))
	}
	g = c()
	game_state = rho(m,S)
	if(game_state$terminal)
	{
		 text(x0,d,labels = c('X','.','O')[game_state$winner+2],cex = 0.75)
		 g = c(g,game_state$winner)
		 return(g)
	}else{
		index = which(m == 0)
		x_seq = seq(lims[1],lims[2],length = length(index)+1)
		x_mid = x_seq[-1]-0.5*diff(x_seq)[1]
		for(i in 1:length(index))
		{
			segments(x0,d,x_mid[i],d+1, col = 'grey')
			text(x_mid[i],d+1,labels = index[i],pos = 1,cex = 0.7)
			x           = m
			x[index[i]] = k
			g           = c(g,draw_game_tree(x,-1*k,x_mid[i],d+1,lims = c(x_seq[i],x_seq[i+1])))
		}
		return(g)
	}
}

res  = draw_game_tree(m,-1)
