# Energy as a function of MC steps

library(ggplot2)

# Parameters
############ 
N = 32           # array size
J = 1           # interaction strength
conv.eq = 1   # convergence to equilibrium
conv = 550      # measurements
reInit = TRUE  # re-initialize for new temperature
TSeq = J*seq(1.2,3.8, by=0.01)  # temperature range

path.FIGS = 'images'
path.DATA = 'data'
source('ising.func.R')

# Array Initialization
######################
E1 = c()
E2 = c()
M1 = c()
M2 = c()
b = 1/0.1
spin = matrix(data=sign(runif(N*N)-0.5), nrow=N)
E1 = c(E1,totalEnergy(N,J))
E2 = c(E2, totalEnergy2(N,J))

NUM = 1:200
for(i in NUM) {
  computeIsingRandExp(conv.eq*N*N, J, b)
  E1 = c(E1,totalEnergy(N,J))
  E2 = c(E2, totalEnergy2(N,J))
  M1 = c(M1, sum(spin))
  M2 = c(M2, sum(spin*spin))
}


df = data.frame(
  steps = c(0,NUM*conv.eq*N*N),
  E1,
  E2,
  DeltaE = (E2-E1*E1)*b*b
)

ggplot(df, aes(steps, E1)) +
  geom_point(col='red', size=3) + 
  theme_bw()

ggplot(df, aes(steps, sqrt(E2))) +
  geom_point(col='red', size=3) + 
  theme_bw()

rasterGraph(spin)
plot(M2-M1)
plot(E1)
plot(E2/2+E1)
