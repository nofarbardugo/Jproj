library(ggplot2)
library(reshape2)
library(seqinr)
require(gtable)
library(gridExtra)

germline.path <- 'C:/Users/Jennifer/Research/plot_properties/germline/'
setwd(path)
organism <- 'human'
chain <- 'VH'

#-------------------------------------------------------------------
# plot average (log) - in decreasing order by 3rd compartment of BM
#-------------------------------------------------------------------
all.rem <- c(10, 15, 19, 23, 27, 45, 46, 48, 49, 50, 51, 52, 13, 32, 14, 47, 68, 31, 12, 63, 44,29,5,62,65)
#-------------------------------
# get order from BM 1st compartment
#-------------------------------
path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/BM/analyzed/'
gene.usage.in <- paste0(path, 'gene_usage/')
avg.size.files <- list.files(gene.usage.in, pattern="V_avg_size_IF_BM.csv*")
V.avg.size <- read.csv(paste0(gene.usage.in, avg.size.files[1]), row.names=1) 
V.avg.size <- as.matrix(V.avg.size[, -all.rem])
V.avg.size <- V.usage/rowSums(V.usage)
V.avg.size.3 <- V.avg.size[3,]
# get decresing order
V.order <- sort(V.avg.size.3,  index.return = TRUE, decreasing = T)$ix
#-------------------------------
# get other BM samples
#-------------------------------
all.df <- melt(V.avg.size[1:5,])
all.df$sem <- 0
colnames(all.df) <- c('Dataset', 'Vgene', 'avg', 'sem')
# average IgM+ in BM
V.avg.size.IgM.BM <- V.avg.size[6:7,]
avg <- colMeans(V.avg.size.IgM.BM)
sem <- apply(V.avg.size.IgM.BM, 2, sd)/sqrt(nrow(V.avg.size.IgM.BM))
df <- data.frame(Dataset='BM - IgM+', Vgene=colnames(V.avg.size.IgM.BM), avg=avg, sem=sem)
all.df <- rbind(all.df, df)
#-------------------------------
# get pfizer
#-------------------------------
path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/pfizer/analyzed/heavy/'
gene.usage.in <- paste0(path, 'gene_usage/')
avg.size.files <- list.files(gene.usage.in, pattern="V_avg_size_IF*")
for (i in 1:length(avg.size.files)){
  V.avg.size <- read.csv(paste0(gene.usage.in, avg.size.files[i]), row.names=1) 
  V.avg.size <- V.avg.size[-grep('136', row.names(V.avg.size)),]
  V.avg.size <- V.avg.size[-grep('273', row.names(V.avg.size)),]
  V.avg.size <- V.avg.size[, -all.rem]
  V.avg.size <- V.avg.size/rowSums(V.avg.size)
  # in IGHM remove first 2 samples - too small
  if (i==3)
    V.avg.size <- V.avg.size[-(1:2),]
  # compute mean and sem
  avg <- colMeans(V.avg.size)
  sem <- apply(V.avg.size, 2, sd)/sqrt(nrow(V.avg.size))
  df <- data.frame(Dataset = paste0('Pfizer_', substr(avg.size.files[i], 15, nchar(avg.size.files[i])-4)), Vgene=colnames(V.avg.size), avg=avg, sem=sem)
  all.df <- rbind(all.df, df)
}

# order Vs
l <- levels(all.df$Vgene)
all.df$Vgene <- factor(all.df$Vgene, levels=l[V.order])
# order datasets
d <- levels(all.df$Dataset)
all.df$Dataset <- factor(all.df$Dataset, levels=d[c(1:6, 10,9,8,7)])
all.df.avg.size <- all.df
# plot all together
p1 <- ggplot(data=all.df, aes(x=Vgene, y=avg))+
  geom_bar(stat="identity", position=position_dodge(.9), colour="black")+
  geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.5, alpha = .5, colour="black", position=position_dodge(.9)) +
  theme_bw()+
  facet_wrap(~Dataset,ncol=1, scales = "free_y")+
  ylab('Average clone size')+
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
# get IF fraction data - for correlations with average clone size
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
df.avg.wide <- dcast(all.df.avg.size, Dataset ~ Vgene, value.var="avg")
avg.size.mat <- as.matrix(df.avg.wide[, 2:44])
dimnames(avg.size.mat) <- list(df.avg.wide[,1], colnames(df.avg.wide)[2:44])
corr.res <- matrix(0, nrow(avg.size.mat),  nrow(avg.size.mat))
corr.p <- matrix(0, nrow(avg.size.mat),  nrow(avg.size.mat))
for(i in 1:nrow(avg.size.mat)){
  for(j in i:nrow(avg.size.mat)){
    tmp <- cor.test(avg.size.mat[i,], IF.mat[j,], method='spearman')
    corr.res[i,j] <- tmp$estimate
    corr.p[i,j] <- tmp$p.value
    tmp <- cor.test(IF.mat[i,], avg.size.mat[j,],  method='spearman')
    corr.res[j,i] <- tmp$estimate
    corr.p[j,i] <- tmp$p.value
  }
}
df.corr.res <- melt(corr.res)
cols <- rainbow(20)
cols <- rev(cols[1:15])
p3 <- ggplot(df.corr.res, aes(df.corr.res[,1], df.corr.res[,2], fill=value)) +
  geom_tile(shape=1)+
  scale_x_discrete(breaks = 1:10, labels=rownames(avg.size.mat))+
  scale_y_reverse(breaks = 10:1, labels=rev(rownames(avg.size.mat)))+
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


grid.arrange(p1, arrangeGrob(p2, p3), ncol = 2)
