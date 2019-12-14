# Graphing Data

path.FIGS = 'images'
path.DATA = 'data'
library(ggplot2)

file.list = file.path(path.DATA, dir(path.DATA,'Tri'))
d = read.csv(file.list[3])
names(d)

d[which(d$T.J<3 & d$M2avg<5.5e9),] <- NA

ggplot(d, aes(T.J, abs(Mavg))) +
  geom_path(col='orange', size=1) + 
  geom_point(col='red', size=2) + 
  xlab('T/J') +
  ylab('|M|') +
  theme_bw()

d$chi = c(0,diff(d$Mavg))


ggplot(d, aes(T.J, chi)) +
  geom_point(col='red', size=2) + 
  xlab('T/J') +
  ylab(expression(paste(chi))) +
  theme_bw()
