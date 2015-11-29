
library(ggplot2)
library(reshape2)
library(seqinr)
require(gtable)
library(gridExtra)
library(pracma)
library(scales)

setwd('C:/Users/Jennifer/GoogleDrive/Nina_new/clone_size_paper/code_for_figures/')

df <- read.csv('clone_size_RNA.csv')

df <- df[-which(df$y==0),]

# remove IgA
df <- df[-which(df$Iso=='IGHA'),]
# order datasets
d <- levels(factor(df$Iso))
df$Iso <- factor(df$Iso, levels=d[c(3,2, 1,4)])

ggplot(data=df, aes(x=x, y=y, group=Sample, color = IF)) +
  geom_line()+
  geom_point()+
  xlab('Measured clone size')+
  ylab('Frequency')+
  theme_bw() +  theme(  legend.title=element_blank(),
          legend.position=c(1,1.1), # legend position
          legend.justification=c(1,1.1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank()
  )+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~Iso, ncol=2, scales = 'free')