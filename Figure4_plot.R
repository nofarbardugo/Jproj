library(ggplot2)
library(reshape2)
library(seqinr)
require(gtable)
library(gridExtra)

path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/BM/analyzed/'
germline.path <- 'C:/Users/Jennifer/Research/plot_properties/germline/'
setwd(path)
organism <- 'human'
chain <- 'VH'

#-------------------------------------------------------------------
# plot number of clones - in decreasing order by 3rd compartment of BM
#-------------------------------------------------------------------
all.rem <- c(10, 15, 19, 23, 27, 45, 46, 48, 49, 50, 51, 52, 13, 32, 14, 47, 68, 31, 12, 63, 44,29,5,62,65)
#-------------------------------
# get order from BM 1st compartment
#-------------------------------
path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/BM/analyzed/'
gene.usage.in <- paste0(path, 'gene_usage/')
V.usage.files <- list.files(gene.usage.in, pattern="V_usage_clone_IF*")
V.usage <- read.csv(paste0(gene.usage.in, V.usage.files[1]), row.names=1) 
V.usage <- as.matrix(V.usage[, -all.rem])
V.usage <- V.usage/rowSums(V.usage)
V.usage.3 <- V.usage[3,]
# get decresing order
V.order <- sort(V.usage.3,  index.return = TRUE, decreasing = T)$ix
#-------------------------------
# get other BM samples
#-------------------------------
all.df <- melt(V.usage[1:5,])
all.df$sem <- 0
colnames(all.df) <- c('Dataset', 'Vgene', 'avg', 'sem')
# average IgM+ in BM
V.usage.IgM.BM <- V.usage[6:7,]
avg <- colMeans(V.usage.IgM.BM)
sem <- apply(V.usage.IgM.BM, 2, sd)/sqrt(nrow(V.usage.IgM.BM))
df <- data.frame(Dataset='BM - IgM+', Vgene=colnames(V.usage.IgM.BM), avg=avg, sem=sem)
all.df <- rbind(all.df, df)
#-------------------------------
# get pfizer
#-------------------------------
path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/pfizer/analyzed/heavy/'
gene.usage.in <- paste0(path, 'gene_usage/')
V.usage.files <- list.files(gene.usage.in, pattern="V_usage_clone_IF*")
for (i in 1:length(V.usage.files)){
  V.usage <- read.csv(paste0(gene.usage.in, V.usage.files[i]), row.names=1) 
  V.usage <- V.usage[-grep('136', row.names(V.usage)),]
  V.usage <- V.usage[-grep('273', row.names(V.usage)),]
  V.usage <- V.usage[, -all.rem]
  V.usage <- V.usage/rowSums(V.usage)
  # compute mean and sem
  avg <- colMeans(V.usage)
  sem <- apply(V.usage, 2, sd)/sqrt(nrow(V.usage))
  df <- data.frame(Dataset = paste0('Pfizer_', substr(V.usage.files[i], 18, nchar(V.usage.files[i])-4)), Vgene=colnames(V.usage), avg=avg, sem=sem)
  all.df <- rbind(all.df, df)
}

# order Vs
l <- levels(all.df$Vgene)
all.df$Vgene <- factor(all.df$Vgene, levels=l[V.order])
# order datasets
d <- levels(all.df$Dataset)
all.df$Dataset <- factor(all.df$Dataset, levels=d[c(1:6, 10,9,8,7)])
all.df.n.clones <- all.df
# plot all together
p1 <- ggplot(data=all.df, aes(x=Vgene, y=avg))+
  geom_bar(stat="identity", position=position_dodge(.9), colour="black")+
  geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.5, alpha = .5, colour="black", position=position_dodge(.9)) +
  theme_bw()+
  facet_wrap(~Dataset,ncol=1, scales = "free_y")+
  ylab('Fraction of clones')+
  theme(  axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank(),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank())


# calculate correlations between datasets
# spearman
df.wide <- dcast(all.df, Dataset ~ Vgene, value.var="avg")
n.mat <- as.matrix(df.wide[, 2:44])
dimnames(n.mat) <- list(df.wide[,1], colnames(df.wide)[2:44])
corr.res <- matrix(0, nrow(n.mat),  nrow(n.mat))
corr.p <- matrix(0, nrow(n.mat),  nrow(n.mat))
for(i in 1:nrow(n.mat)){
  for(j in i:nrow(n.mat)){
    tmp <- cor.test(n.mat[i,], n.mat[j,], method='spearman')
    corr.res[i,j] <- tmp$estimate
    corr.res[j,i] <- tmp$estimate
    corr.p[i,j] <- tmp$p.value
    corr.p[j,i] <- tmp$p.value
  }
}
#----------------------------------------------
# heatmap of correlations 
#----------------------------------------------
df.corr.res <- melt(corr.res)
cols <- rainbow(20)
cols <- rev(cols[1:15])
p2 <- ggplot(df.corr.res, aes(df.corr.res[,1], df.corr.res[,2], fill=value)) +
  geom_tile(shape=1)+
  scale_x_discrete(breaks = 1:10, labels=rownames(n.mat))+
  scale_y_reverse(breaks = 10:1, labels=rev(rownames(n.mat)))+
  scale_fill_gradientn(limits=c(-1,1),colours=cols)+
  coord_cartesian(xlim = c(0.5,10.5), ylim=c(0.5,10.5))+
  ggtitle('Number of clones - Spearman (Rho)')+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=90, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# #----------------------------------------------
# # heatmap of p values of correlations - not needed for plot
# #----------------------------------------------
# df.corr.p <- melt(corr.p)
# p11<-ggplot(df.corr.p, aes(df.corr.p[,1], df.corr.p[,2], fill=value)) +
#   geom_tile(shape=1)+
#   scale_x_discrete(breaks = 1:13, labels=rownames(size.mat))+
#   scale_y_reverse(breaks = 13:1, labels=rev(rownames(size.mat)))+
#   scale_fill_gradientn(limits=c(0,1),colours=cols)+
#   coord_cartesian(xlim = c(0.5,13.5), ylim=c(0.5,13.5))+
#   ggtitle('Number of clones - Spearman (p-value)')+
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle=90, vjust=1),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank())

#-------------------------------------------------------------------
# get IF fraction data - for correlations with number of clones
#-------------------------------------------------------------------
# get BM
#-------------------------------
path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/BM/analyzed/'
gene.usage.in <- paste0(path, 'gene_usage/')
IF.frac.files <- list.files(gene.usage.in, pattern="V_IF_frac_abs_*")
IF.frac <- as.matrix(read.csv(paste0(gene.usage.in, IF.frac.files[1]), row.names=1))
IF.frac <- IF.frac[,-all.rem]
all.df <- melt(IF.frac[1:5,])
colnames(all.df) <- c('Dataset', 'Vgene', 'If.frac')
# average IgM+ in BM
IF.frac.IgM.BM <- IF.frac[6:7,]
avg <- colMeans(IF.frac.IgM.BM)
df <- data.frame(Dataset='BM - IgM+', Vgene=colnames(IF.frac.IgM.BM), If.frac=avg)
all.df <- rbind(all.df, df)
#-------------------------------
# get pfizer
#-------------------------------
path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/pfizer/analyzed/heavy/'
gene.usage.in <- paste0(path, 'gene_usage/')
IF.frac.files <- list.files(gene.usage.in, pattern="V_IF_frac_abs*")
for (i in 1:length(IF.frac.files)){
  IF.frac <- read.csv(paste0(gene.usage.in, IF.frac.files[i]), row.names=1) 
  IF.frac <- IF.frac[-grep('136', row.names(IF.frac)),]
  IF.frac <- IF.frac[-grep('273', row.names(IF.frac)),]
  IF.frac <- IF.frac[, -all.rem]
  # compute mean and sem
  avg <- colMeans(IF.frac)
  df <- data.frame(Dataset = substr(IF.frac.files[i], 15, nchar(IF.frac.files[i])-4), Vgene=colnames(IF.frac), If.frac=avg)
  all.df <- rbind(all.df, df)
}

# order Vs
l <- levels(all.df$Vgene)
all.df$Vgene <- factor(all.df$Vgene, levels=l[V.order])
# order datasets
d <- levels(all.df$Dataset)
all.df$Dataset <- factor(all.df$Dataset, levels=d[c(1:6, 10,9,8,7)])
all.df.IF.frac <- all.df
#-----------------------------------------------------------------------
# Spearman's correlation between number of clones and IF frac between datasets
df.IF.wide <- dcast(all.df.IF.frac, Dataset ~ Vgene, value.var="If.frac")
IF.mat <- as.matrix(df.IF.wide[, 2:44])
dimnames(IF.mat) <- list(df.IF.wide[,1], colnames(df.IF.wide)[2:44])
df.n.clones.wide <- dcast(all.df.n.clones, Dataset ~ Vgene, value.var="avg")
n.clones.mat <- as.matrix(df.n.clones.wide[, 2:44])
dimnames(n.clones.mat) <- list(df.n.clones.wide[,1], colnames(df.n.clones.wide)[2:44])
corr.res <- matrix(0, nrow(n.clones.mat),  nrow(n.clones.mat))
corr.p <- matrix(0, nrow(n.clones.mat),  nrow(n.clones.mat))
for(i in 1:nrow(n.clones.mat)){
  for(j in i:nrow(n.clones.mat)){
    tmp <- cor.test(n.clones.mat[i,], IF.mat[j,], method='spearman')
    corr.res[i,j] <- tmp$estimate
    corr.p[i,j] <- tmp$p.value
    tmp <- cor.test(IF.mat[i,], n.clones.mat[j,],  method='spearman')
    corr.res[j,i] <- tmp$estimate
    corr.p[j,i] <- tmp$p.value
  }
}
df.corr.res <- melt(corr.res)
cols <- rainbow(20)
cols <- rev(cols[1:15])
p3 <- ggplot(df.corr.res, aes(df.corr.res[,1], df.corr.res[,2], fill=value)) +
  geom_tile(shape=1)+
  scale_x_discrete(breaks = 1:10, labels=rownames(n.mat))+
  scale_y_reverse(breaks = 10:1, labels=rev(rownames(n.mat)))+
  scale_fill_gradientn(limits=c(-1,1),colours=cols)+
  coord_cartesian(xlim = c(0.5,10.5), ylim=c(0.5,10.5))+
  ggtitle('Number of clones vs. IF fraction - Spearman (Rho)')+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=90, vjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# #----------------------------------------------
# # heatmap of p values of correlations - not needed for plot
# #----------------------------------------------
# df.corr.p <- melt(corr.p)
# p11<-ggplot(df.corr.p, aes(df.corr.p[,1], df.corr.p[,2], fill=value)) +
#   geom_tile(shape=1)+
#   scale_x_discrete(breaks = 1:13, labels=rownames(size.mat))+
#   scale_y_reverse(breaks = 13:1, labels=rev(rownames(size.mat)))+
#   scale_fill_gradientn(limits=c(0,1),colours=cols)+
#   coord_cartesian(xlim = c(0.5,13.5), ylim=c(0.5,13.5))+
#   ggtitle('Number of clones vs. IF fraction - Spearman (p-value)')+
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle=90, vjust=1),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank())


#----------------------------
# IF fraction in BM plot
#----------------------------

path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/BM/analyzed/Clones/'
germline.path <- 'C:/Users/Jennifer/Research/plot_properties/germline/'
setwd(path)
organism <- 'human'
chain <- 'VH'


files <- list.files(path)
files.names <- sapply(files, function(x) substr(x, 1, nchar(x)-15))
IF.frac <- data.frame(Sample=files.names, IF.frac =numeric(length(files)), sem=0)
for(i in 1:length(files)){
  df <- read.csv(paste0(path, files[i]), row.names=1) 
  IF.frac[i, 'IF.frac'] <- length(which(df$INFRAME==T))/nrow(df)
}

p4 <- ggplot(IF.frac, aes(x=Sample, y=IF.frac))+
  geom_bar(stat="identity")+
  coord_cartesian(ylim = c(0,1))+
  ylab('IF fraction')+
  geom_errorbar(aes(ymin=IF.frac-sem, ymax=IF.frac+sem), width=.1) + 
  theme(  panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank())

grid.arrange(p1, arrangeGrob(p2, p3, p4), ncol = 2)

