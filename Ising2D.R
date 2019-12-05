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

# Parameters
############
N = 80           # array size
J = 1           # interaction strength
conv.eq = 1000   # convergence to equilibrium
conv = 1000      # measurements
reInit = FALSE  # re-initialize for new temperature
TSeq = J*seq(1.2,3.8, by=0.02)  # temperature range

path.FIGS = 'images'
path.DATA = 'data'
file.runTime = file.path(path.DATA,'runTimes.csv')
if(file.exists(file.runTime)) { 
  d.runTimeAll = read.csv(file.runTime, stringsAsFactors = FALSE) 
} else {
  d.runTimeAll = data.frame()
}

source('ising.func.R')

# Array Initialization
######################
spin = matrix(data=sign(runif(N*N)-0.5), nrow=N)

# Sample Output
###############
rasterGraph(spin)
ggsave(file.path(path.FIGS,paste0('Ising2D-',N,'x',N,'-Random.png')), width=4,height=4,dpi=220)
computeIsingRandExp(conv.eq*N*N, J, 1/2.4)
totalEnergy(N,J)
rasterGraph(spin)
ggsave(file.path(path.FIGS,paste0('Ising2D-',N,'x',N,'-Domains.png')), width=4,height=4,dpi=220)

# Computation Intesive Part: M vs T
##################################
d.runTime = data.frame(N,J,conv.eq,conv, reInit, date=Sys.Date(), 
                       start.time=as.numeric(Sys.time()),
                       end.time=0,diff.s=0)
bSeq = 1/TSeq
result = data.frame()
for(b in bSeq) {
  print(paste("Temp: ",1/b))
  if (reInit) { spin = matrix(data=sign(runif(N*N)-0.5), nrow=N) }
  computeIsingRandExp(conv.eq*N*N, J, b)
  Mavg = 0
  M2avg = 0
  for(i in 1:conv) {
    computeIsingRandExp(N*N, J, b)
    Ms = sum(spin)
    Mavg = Mavg + Ms
    M2avg = M2avg + Ms*Ms
  }
  result = rbind(result, data.frame(b, conv, Mavg, M2avg,
                                    T.J = 1/b/J,
                                    chi=(M2avg/(conv*N*N) - Mavg*Mavg/(conv*conv*N*N))*b))
}

# save timing
d.runTime$end.time = as.numeric(Sys.time())
d.runTime$diff.s = d.runTime$end.time-d.runTime$start.time
d.runTimeAll = rbind(d.runTimeAll,d.runTime)
write.csv(d.runTimeAll,file=file.runTime,row.names = FALSE)

# Graphing of Data
##################
ggplot(result, aes(T.J, abs(Mavg)/(conv*N*N))) +
  geom_point(col='red', size=2) + 
  ggtitle(paste('N=',N,'x',N,' conv=',conv, ' reInit=',reInit)) + 
  xlab('T/J') +
  ylab('|M|') +
  geom_vline(xintercept=2.27, stroke=5,size=2, alpha=0.5) + 
  theme_bw()
ggsave(file.path(path.FIGS,paste0('Ising2D-',N,'x',N,'-c',conv,'.png')), width=4, height=3, dpi=220)
write.csv(result,file.path(path.DATA,paste0('Ising2D-',N,'x',N,'-c',conv,'.csv')), row.names=FALSE)

ggplot(result, aes(T.J, chi)) +
  geom_smooth(span=0.2)+
  geom_point(col='red', size=2) + 
  ggtitle(paste('N=',N,'x',N,' conv=',conv, ' reInit=',reInit)) + 
  xlab('T/J') +
  ylab(expression(paste(chi))) +
  geom_vline(xintercept=2.27, stroke=5,size=2, alpha=0.5) + 
  theme_bw()
ggsave(file.path(path.FIGS,paste0('Ising2D-',N,'x',N,'-c',conv,'-Chi.png')), width=4, height=3, dpi=220)

