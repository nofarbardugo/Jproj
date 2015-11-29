''''
library(ggplot2)
library(reshape2)
library(seqinr)
require(gtable)
library(gridExtra)



# VH usage plot - all clones and top 100 clones
# Pfizer and flu data

all.rem <- c(10, 15, 19, 23, 27, 45, 46, 48, 49, 50, 51, 52, 13, 32, 14, 47, 68, 31, 12, 63, 44,29,5,62,65)

# Pfizer
path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/pfizer/analyzed/heavy/'
germline.path <- 'C:/Users/Jennifer/Research/plot_properties/germline/'
setwd(path)
organism <- 'human'
chain <- 'VH'

gene.usage.in <- paste0(path, 'gene_usage/')
V.usage.files <- list.files(gene.usage.in, pattern="V_usage_clone_IF*")
for (i in 1:length(V.usage.files)){
  V.usage <- read.csv(paste0(gene.usage.in, V.usage.files[i]), row.names=1) 
  # remove bad donors
  V.usage <- V.usage[-grep('136', row.names(V.usage)),]
  V.usage <- V.usage[-grep('273', row.names(V.usage)),] 
  # remove unseen V genes
  V.usage <- V.usage[,-all.rem]
  V.genes <- colnames(V.usage)
  V.genes <- gsub(".", "-", V.genes, fixed = TRUE)
  # normalize usage
  V.usage.norm <- V.usage/rowSums(V.usage)  
  avg <- colMeans(V.usage.norm)
  sem <- apply(V.usage.norm, 2, sd)/sqrt(nrow(V.usage.norm))
  tmp.df <- data.frame(Vgene=colnames(V.usage), avg=avg, sem=sem)
  tmp.df$Iso <- substr(V.usage.files[i], 18, nchar(V.usage.files[i])-4)
  if (i==1)
    df <- tmp.df
  else
    df <- rbind(df, tmp.df)
}
# order genes by usage sum
means <- with(df, tapply(avg, Vgene, mean, na.rm = T))
V.order <- order(-means)  
df$Vgene <- factor(df$Vgene)
df$Vgene <- factor(df$Vgene,df$Vgene[V.order])
df$Iso <- factor(df$Iso)
df$Iso <- factor(df$Iso, levels=rev(levels(df$Iso)))
 
# bar plot with error bars
p1<-ggplot(data=df, aes(x=Vgene, y=avg, fill=Iso))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
 # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('All clones - Pfizer')+
#  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=2))

# Pfizer - top 100 clones

path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/pfizer/analyzed/heavy/Clones/'
setwd(path)
dirs <- list.files(path)
for(i in 1:length(dirs)){
  files <- list.files(paste0(path, dirs[i]))
  # remove bad donors
  files <- files[-grep('136', files)]
  files <- files[-grep('273', files)] 
  V.usage <- matrix(0, nrow = length(files), ncol = length(V.genes))
  for(j in 1:length(files)){
    df <- read.csv(paste0(path, dirs[i], '/', files[j]))
    df <- df[df$FUNCTIONAL==T,]
    a <- sort(df$CP_NUM, decreasing = T,  index.return = T)
    df2 <- df[a$ix,]
    df2 <- df2[1:100,]
    # remove allele anotation
    df2$V_CALL <- unlist(lapply(as.character(df2$V_CALL), function(x) substr(x,1,nchar(x)-3)))
    # get V usage for sample
    for(k in 1:length(V.genes))
      V.usage[j,k] <- length(which(df2$V_CALL==V.genes[k]))
  }
  # normalize usage
  V.usage.norm <- V.usage/rowSums(V.usage)  
  avg <- colMeans(V.usage.norm)
  sem <- apply(V.usage.norm, 2, sd)/sqrt(nrow(V.usage.norm))
  tmp.df <- data.frame(Vgene=V.genes, avg=avg, sem=sem)
  tmp.df$Iso <- dirs[i]
  if (i==1)
    all.df <- tmp.df
  else
    all.df <- rbind(all.df, tmp.df)
}

df <- all.df
# order genes by usage sum of all clones
df$Vgene <- factor(df$Vgene,df$Vgene[V.order])
df$Iso <- factor(df$Iso)
df$Iso <- factor(df$Iso, levels=rev(levels(df$Iso)))
# remove gene with low frequencies
if(!is.null(threshold))
  df <- df[-which(df$sum<threshold),]

# bar plot with error bars
p2<-ggplot(data=df, aes(x=Vgene, y=avg, fill=Iso))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - Pfizer')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  theme(  legend.title=element_blank(),
          legend.position='none', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )
# normalize beofre or after removeing bad Vs?????????
#-----------------------------------------------------------------------------------
# Flu
path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/flu/analyzed/'
germline.path <- 'C:/Users/Jennifer/Research/plot_properties/germline/'
setwd(path)
organism <- 'human'
chain <- 'VH'

gene.usage.in <- paste0(path, 'gene_usage/')
V.usage.files <- list.files(gene.usage.in, pattern="V_usage_clone_IF*")
V.usage.files <- V.usage.files[-3]
for (i in 1:length(V.usage.files)){
  V.usage <- read.csv(paste0(gene.usage.in, V.usage.files[i]), row.names=1) 
  # take only pre-vaccine samples
  pre.ind <- c(grep('RL013', row.names(V.usage)), grep('RL014', row.names(V.usage)), grep('RL015', row.names(V.usage)))
  V.usage <- V.usage[pre.ind,]
  # remove unseen V genes
  V.usage <- V.usage[,-all.rem]
  # normalize usage
  V.usage.norm <- V.usage/rowSums(V.usage)  
  avg <- colMeans(V.usage.norm)
  sem <- apply(V.usage.norm, 2, sd)/sqrt(nrow(V.usage.norm))
  tmp.df <- data.frame(Vgene=colnames(V.usage), avg=avg, sem=sem)
  tmp.df$Iso <- substr(V.usage.files[i], 18, nchar(V.usage.files[i])-4)
  if (i==1)
    df <- tmp.df
  else
    df <- rbind(df, tmp.df)
}
# average over IgG1 and 2
tmp.df <- df[which(df$Iso=='IGHG-1'| df$Iso=='IGHG-2'),]
means <- with(tmp.df, tapply(avg, Vgene, mean, na.rm = T))
df[which(all.df$Iso=='IGHG-1'), 'avg'] <- means
df[which(all.df$Iso=='IGHG-1'), 'Iso'] <- 'IGHG'
df <- df[-which(df$Iso=='IGHG-2'),]
# order genes by usage sum
means <- with(df, tapply(avg, Vgene, mean, na.rm = T))
V.order <- order(-means)  
df$Vgene <- factor(df$Vgene)
df$Vgene <- factor(df$Vgene,df$Vgene[V.order])
df$Iso <- factor(df$Iso)
df$Iso <- factor(df$Iso, levels=rev(levels(df$Iso)))
# remove gene with low frequencies
if(!is.null(threshold))
  df <- df[-which(df$sum<threshold),]

# bar plot with error bars
p3<-ggplot(data=df, aes(x=Vgene, y=avg, fill=Iso))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('All clones - Flu')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=2))


# flu - top 100 clones

path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/flu/analyzed/Clones/'
setwd(path)
dirs <- list.files(path)
dirs <- dirs[-3]
for(i in 1:length(dirs)){
  files <- list.files(paste0(path, dirs[i]))
  # take only pre-vaccine samples
  pre.ind <- c(grep('RL013', files), grep('RL014', files), grep('RL015', files))
  files <- files[pre.ind]
  V.usage <- matrix(0, nrow = length(files), ncol = length(V.genes))
  for(j in 1:length(files)){
    df <- read.csv(paste0(path, dirs[i], '/', files[j]))
    df <- df[df$FUNCTIONAL==T,]
    a <- sort(df$CP_NUM, decreasing = T,  index.return = T)
    df2 <- df[a$ix,]
    df2 <- df2[1:100,]
    # remove allele anotation
    df2$V_CALL <- unlist(lapply(as.character(df2$V_CALL), function(x) substr(x,1,nchar(x)-3)))
    # get V usage for sample
    for(k in 1:length(V.genes))
      V.usage[j,k] <- length(which(df2$V_CALL==V.genes[k]))
  }
  # normalize usage
  V.usage.norm <- V.usage/rowSums(V.usage)  
  avg <- colMeans(V.usage.norm)
  sem <- apply(V.usage.norm, 2, sd)/sqrt(nrow(V.usage.norm))
  tmp.df <- data.frame(Vgene=V.genes, avg=avg, sem=sem)
  tmp.df$Iso <- dirs[i]
  if (i==1)
    all.df <- tmp.df
  else
    all.df <- rbind(all.df, tmp.df)
}
# average over IgG1 and 2
tmp.df <- all.df[which(all.df$Iso=='IGHG-1'| all.df$Iso=='IGHG-2'),]
means <- with(tmp.df, tapply(avg, Vgene, mean, na.rm = T))
all.df[which(all.df$Iso=='IGHG-1'), 'avg'] <- means
all.df[which(all.df$Iso=='IGHG-1'), 'Iso'] <- 'IGHG'
all.df <- all.df[-which(all.df$Iso=='IGHG-2'),]
df <- all.df
# order genes by usage sum of all clones
df$Vgene <- factor(df$Vgene,df$Vgene[V.order])
df$Iso <- factor(df$Iso)
df$Iso <- factor(df$Iso, levels=rev(levels(df$Iso)))
# remove gene with low frequencies
if(!is.null(threshold))
  df <- df[-which(df$sum<threshold),]

# bar plot with error bars
p4<-ggplot(data=df, aes(x=Vgene, y=avg, fill=Iso))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - Flu')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
         # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )

grid.arrange(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4), ncol=2)
