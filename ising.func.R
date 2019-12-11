# Computation
#############
# basic Ising computation step
computeIsing <- function(num.iter, J, beta) {
  for(i in 1:num.iter) {
    # choose random spin
    x=round(runif(1,min=1,max=N))
    y=round(runif(1,min=1,max=N))
    # compute energy change to flip:
    nb = spin[(x %% N)+1,y] + spin[((x-2) %% N)+1,y] + 
      spin[x,(y %% N)+1] + spin[x,((y-2) %% N)+1]
    dE = 2*J*spin[x,y]*nb
    if (dE<0) { 
      spin[x,y] <<- -spin[x,y] 
    } else {
      # flip coin
      if (runif(1) < exp(-dE*beta)) {
        spin[x,y] <<-  -spin[x,y]
      }
    }
  }
}


computeIsingRandExp <- function(num.iter, J, beta) {
  xa=round(runif(num.iter,min=1,max=N))
  ya=round(runif(num.iter,min=1,max=N))
  # precompute all 5 possible values
  expVal = exp(-2*J*beta*(-4:4))  
  for(i in 1:num.iter) {
    # choose random spin
    x=xa[i]
    y=ya[i]
    # compute energy change to flip:
    nb = spin[(x %% N)+1,y] + spin[((x-2) %% N)+1,y] + 
      spin[x,(y %% N)+1] + spin[x,((y-2) %% N)+1]
    dE = spin[x,y]*nb
    if (dE<0) { 
      spin[x,y] <<- -spin[x,y] 
    } else if (runif(1) < expVal[(dE+5)]) {
      spin[x,y] <<-  -spin[x,y]
    }
  }
}




computeIsingRand <- function(num.iter, J, beta) {
  xa=round(runif(num.iter,min=1,max=N))
  ya=round(runif(num.iter,min=1,max=N))
  for(i in 1:num.iter) {
    # choose random spin
    x=xa[i]
    y=ya[i]
    # compute energy change to flip:
    nb = spin[(x %% N)+1,y] + spin[((x-2) %% N)+1,y] + 
      spin[x,(y %% N)+1] + spin[x,((y-2) %% N)+1]
    dE = 2*J*spin[x,y]*nb
    if (dE<0) { 
      spin[x,y] <<- -spin[x,y] 
    } else {
      # flip coin
      if (runif(1) < exp(-dE*beta)) {
        spin[x,y] <<-  -spin[x,y]
      }
    }
  }
}

computeIsing1DRand <- function(num.iter, J, beta) {
  xa=round(runif(num.iter,min=1,max=N))
  ya=round(runif(num.iter,min=1,max=N))
  for(i in 1:num.iter) {
    # choose random spin
    x=xa[i]
    y=ya[i]
    # compute energy change to flip:
    nb = spin[(x %% N)+1,y] + spin[((x-2) %% N)+1,y] 
    dE = 2*J*spin[x,y]*nb
    if (dE<0) { 
      spin[x,y] <<- -spin[x,y] 
    } else {
      # flip coin
      if (runif(1) < exp(-dE*beta)) {
        spin[x,y] <<-  -spin[x,y]
      }
    }
  }
}


computeIsing2DTriangularRand <- function(num.iter, J, beta) {
  xa=round(runif(num.iter,min=1,max=N))
  ya=round(runif(num.iter,min=1,max=N))
  for(i in 1:num.iter) {
    # choose random spin
    x=xa[i]
    y=ya[i]
    # compute energy change to flip:
    # top, top-right, right
    nb = spin[x,(y %% N)+1] + spin[(x %% N)+1,(y %% N)+1] + spin[(x %% N)+1,y]
    dE = 2*J*spin[x,y]*nb
    if (dE<0) { 
      spin[x,y] <<- -spin[x,y] 
    } else {
      # flip coin
      if (runif(1) < exp(-dE*beta)) {
        spin[x,y] <<-  -spin[x,y]
      }
    }
  }
}



totalEnergy <- function(N,J) {
  totE = 0
  for(x in 1:N) {
    for(y in 1:N) {
      nb = spin[(x %% N)+1,y] + spin[((x-2) %% N)+1,y] + 
        spin[x,(y %% N)+1] + spin[x,((y-2) %% N)+1]
      totE = totE + J*spin[x,y]*nb
    }
  }
  totE / 4
}


rasterGraph <- function(spinMatrix) {
  df = expand.grid(x = 1:N, y = 1:N)
  df$spin = as.vector(spinMatrix)
  ggplot(df, aes(x,y, fill=spin)) + geom_raster() + 
    scale_x_continuous(expand=c(0,0))+
    scale_fill_gradientn(colors=c("red","blue"))+
    scale_y_continuous(expand=c(0,0))+
    theme_bw() + theme(legend.position='none')
}




