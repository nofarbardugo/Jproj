
library(ggplot2)
library(reshape2)
library(seqinr)
require(gtable)
library(gridExtra)
library(pracma)
library(scales)

# Pfizer
germline.path <- 'C:/Users/Jennifer/Research/plot_properties/germline/'
organism <- 'human'
chain <- 'VH'

path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/pfizer/analyzed/heavy/Clones/'
setwd(path)
dirs <- list.files(path)
for(i in 1:length(dirs)){
  print(i)
  files <- list.files(paste0(path, dirs[i]))
  # remove bad donors
  files <- files[-grep('136', files)]
  files <- files[-grep('273', files)] 
  D.index <- matrix(0, 100, length(files))
  dimnames(D.index) <- list(as.character(1:100), files)
  for(j in 1:length(files)){
    df <- read.csv(paste0(path, dirs[i], '/', files[j]))
    df <- df[df$FUNCTIONAL==T,]
    
    # get clone size frequencies
    res <- hist(df$CP_NUM, breaks=logspace(0,4,30), plot=F)  
    p.clone.size <- res$counts/sum(res$counts)
    
    # remove freq==0   
    p.clone.size <- p.clone.size[p.clone.size != 0]
    
    # compute hill diversity index for different Qs
    q <- 1:100
    for(k in q){
      if(k==1)
        D.index[k, j] <- exp(-sum(p.clone.size*log(p.clone.size)))
      else
        D.index[k, j] <- sum(p.clone.size^k)^(1/(1-k))
    }
  }
  tmp.df <- data.frame(Q=q, Iso=dirs[i], D=rowMeans(D.index), sem = (apply(D.index,1,sd)/sqrt(length(files))))
  if (i==1)
    plot.df <- tmp.df
  else
    plot.df <- rbind(plot.df, tmp.df)
}

factor(df$Iso)
l <- levels(factor(plot.df$Iso))
plot.df$Iso <- factor(plot.df$Iso, levels=rev(l))
p1 <- ggplot(data=plot.df, aes(x=Q, y=D, color=Iso, group=Iso)) +
  geom_errorbar(aes(ymin=D-sem, ymax=D+sem), width=.03) +
  geom_line() +
  xlab('Order of diversity (q)')+
  ylab('D index')+
  theme_bw() +
  ggtitle('Hill\'s diversity index')+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank()
  )+
  scale_x_log10()
  

setwd('C:/Users/Jennifer/GoogleDrive/Nina_new/clone_size_paper/code_for_figures/')
df <- read.csv('gini_plot_data.csv')

# order datasets
d <- levels(factor(df$Iso))
df$Iso <- factor(df$Iso, levels=d[c(4, 3,2, 1)])

p2 <- ggplot(data=df, aes(x=x, y=y, group=Sample, color = IF)) +
  geom_line()+
  xlab('Fraction of clones')+
  ylab('Fraction of cells')+
  theme_bw() +  theme(  legend.title=element_blank(),
                        legend.position=c(1,1.1), # legend position
                        legend.justification=c(1,1.1),
                        panel.grid.minor=element_blank(), # remove grid
                        panel.grid.major=element_blank()
  )+
  scale_x_log10()+
  facet_wrap(~Iso, ncol=2, scales = 'free')

grid.arrange(p1, p2, ncol = 2)
