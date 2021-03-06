#########################################
# Compute the temperature dependence of the
# Magnetization for the Ising 2D model 
# with random DEFECTS
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
N = 16           # array size
N.defects = 26    # number of defects, percentage = N.defects / N*N
N.defects / (N*N)*100
J = 1            # interaction strength
conv.eq = 5001   # convergence to equilibrium
conv = 2001      # measurements
reInit = TRUE    # re-initialize for new temperature
TSeq = J*seq(1,3, by=0.05)  # temperature range

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
J.matrix =  matrix(data=J, nrow=N, ncol=N)
x = runif(N.defects,1,N*N)
J.matrix[x]=0
rasterGraph(J.matrix)

# Sample Output
###############
rasterGraph(spin)
beta = 1/1.0
computeIsing2DDefectRand(conv.eq*N*N, J.matrix, beta)
totalEnergy(N,J)
rasterGraph(spin)

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
  computeIsing2DDefectRand(conv.eq*N*N, J.matrix, b)
  Mavg = 0
  M2avg = 0
  for(i in 1:conv) {
    computeIsing2DDefectRand(N*N, J.matrix, b)
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
d.runTime$type = '2D Defect'
d.runTime$step = TSeq[2]-TSeq[1]
d.runTime$TempStart = TSeq[1]
d.runTimeAll = rbind(d.runTimeAll,d.runTime)

write.csv(d.runTimeAll,file=file.runTime,row.names = FALSE)

# Graphing of Data
##################
ggplot(result, aes(T.J, abs(Mavg)/(conv*N*N))) +
  geom_point(col='red', size=2) + 
  ggtitle(paste('N=',N,'x',N,' conv=',conv, ' N.defects=',N.defects)) + 
  xlab('T/J') +
  ylab('|M|') +
  theme_bw()
ggsave(file.path(path.FIGS,paste0('Ising2D-Defect-',N,'x',N,'-c',conv,'.png')), width=4, height=3, dpi=220)
write.csv(result,file.path(path.DATA,paste0('Ising2D-Defect-',N,'x',N,'-c',conv,'.csv')), row.names=FALSE)

ggplot(result, aes(T.J, chi)) +
  geom_smooth(span=0.2)+
  geom_point(col='red', size=2) + 
  ggtitle(paste('N=',N,'x',N,' conv=',conv, ' reInit=',reInit)) + 
  xlab('T/J') +
  ylab(expression(paste(chi))) +
  theme_bw()
ggsave(file.path(path.FIGS,paste0('Ising2D-Defect-',N,'x',N,'-c',conv,'-Chi.png')), width=4, height=3, dpi=220)

