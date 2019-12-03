#########################################
# Compute the temperature dependence of the
# Magnetization for the Ising 2D model and
# test the speed in R
#
# see https://arxiv.org/pdf/0803.0217.pdf
# see http://micro.stanford.edu/~caiwei/me334/Chap12_Ising_Model_v04.pdf 
# see https://github.com/basilwong/monte-carlo-2D-ising
#
# (c) 2019 Thomas Gredig
#########################################


# Ising2D Model Parameters
##########################
library(ggplot2)
library(raster)

N = 10  # array size
J = 1   # interaction strength
beta = 3  # inverse temperature
conv = 20000   # convergence factor
reInit = FALSE # re-initialize random matrix
path.FIGS = 'images'
path.DATA = 'data'
file.runTime = file.path(path.DATA,'runTimes.csv')

# Computation
#############
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



rasterGraph <- function(spinMatrix) {
  df = expand.grid(x = 1:N, y = 1:N)
  df$spin = as.vector(spinMatrix)
  ggplot(df, aes(x,y, fill=spin)) + geom_raster() + 
    scale_x_continuous(expand=c(0,0))+
    scale_fill_gradientn(colors=c("blue","red"))+
    scale_y_continuous(expand=c(0,0))+
    theme_bw() + theme(legend.position='none')
}


# Array Initialization
######################
spin = matrix(data=sign(runif(N*N)-0.5), nrow=N)

# Sample Output
###############
print(rasterGraph(spin))
ggsave(file.path(path.FIGS,paste0('Ising2D-',N,'x',N,'-Random.png')), width=4,height=4,dpi=220)
computeIsing(conv*N*N/100, J, beta)
print(rasterGraph(spin))
ggsave(file.path(path.FIGS,paste0('Ising2D-',N,'x',N,'-Domains.png')), width=4,height=4,dpi=220)

# Computation Intesive Run: M vs T
##################################
d.runTimeAll = read.csv(file.runTime, stringsAsFactors = FALSE)
d.runTime = data.frame(N,conv,date=Sys.Date(), start.time=as.numeric(Sys.time()),end.time=0,diff.s=0)
Mavg = c()
M2avg = c()
TSeq = seq(0.5,4, by=0.05)
bSeq = 1/TSeq
for(b in bSeq) {
  print(b)
  if(reInit) { spin = matrix(data=sign(runif(N*N)-0.5), nrow=N) }
  computeIsing(conv*N*N, J, b)
  Mavg = c(Mavg, sum(spin))
  M2avg = c(M2avg, sum(spin*spin))
}
d.runTime$end.time = as.numeric(Sys.time())
d.runTime$diff.s = d.runTime$end.time-d.runTime$start.time
d.runTime$reInit = reInit
d.runTimeAll = rbind(d.runTimeAll,d.runTime)
write.csv(d.runTimeAll,file=file.runTime,row.names = FALSE)

# Graphing of Data
##################
N=50
conv = 800
d = read.csv(file.path(path.DATA,paste0('Ising2D-',N,'x',N,'-c',conv,'.csv')))
d = data.frame(
  beta = bSeq,
  T.J = 1/bSeq/J,
  J,
  N,
  conv,
  M = Mavg,
  M2 = M2avg,
  chi = ((M2avg*M2avg/(N*N))-(Mavg*Mavg/(N*N)))
)
ggplot(d, aes(T.J, abs(M))) +
  geom_point(col='red', size=2) + 
  ggtitle(paste('N=',N,'x',N,' conv=',conv, ' reInit=',reInit)) + 
  xlab('T/J') +
  ylab('|M|') +
  theme_bw()
ggsave(file.path(path.FIGS,paste0('Ising2D-',N,'x',N,'-c',conv,'.png')), width=4, height=3, dpi=220)
write.csv(d,file.path(path.DATA,paste0('Ising2D-',N,'x',N,'-c',conv,'.csv')), row.names=FALSE)

d$dM = c(0,abs(diff(d$M)))
ggplot(d, aes(T.J, dM)) +
  geom_point(col='red', size=2) + 
  ggtitle(paste('N=',N,'x',N,' conv=',conv, ' reInit=',reInit)) + 
  xlab('T/J') +
  ylab(expression(paste(chi))) +
  theme_bw()
