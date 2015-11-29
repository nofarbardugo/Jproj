
library(ggplot2)
library(reshape2)
library(seqinr)
require(gtable)
library(gridExtra)
require(plyr)



# VH usage plot - all clones and top 100 clones
# Pfizer and flu data

all.rem <- c(10, 15, 19, 23, 27, 45, 46, 48, 49, 50, 51, 52, 13, 32, 14, 47, 68, 31, 12, 63, 44,29,5,62,65)

# Pfizer
path <-"/home/nofar/Desktop/Lab"
#path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/pfizer/analyzed/heavy/'
germline.path <- 'C:/Users/Jennifer/Research/plot_properties/germline/'
germline.path <- path
setwd(path)
organism <- 'human'
chain <- 'VH'

gene.usage.in <- paste0(path, '/data_for_V_plot/Vusage/naive_b/')
V.usage.files <- list.files(gene.usage.in, pattern="V_usage_clone_IF*")

V.usage <- read.csv(paste0(gene.usage.in, V.usage.files), row.names=1) 
# remove unseen V genes
V.usage <- V.usage[,-all.rem]
V.genes <- colnames(V.usage)
V.genes <- gsub(".", "-", V.genes, fixed = TRUE)
# normalize usage
V.usage.norm <- V.usage/rowSums(V.usage)  
avg <- colMeans(V.usage.norm)
sem <- apply(V.usage.norm, 2, sd)/sqrt(nrow(V.usage.norm))
#sem <- apply(sqrt(nrow(V.usage.norm))) # only one sample - 
tmp.df <- data.frame(Vgene=colnames(V.usage), avg=avg, sem=sem)
tmp.df$Iso <- substr(V.usage.files, 18, nchar(V.usage.files)-4)
df <- tmp.df

# order genes by usage sum
means <- with(df, tapply(avg, Vgene, mean, na.rm = T))
V.order <- order(-means)  
df$Vgene <- factor(df$Vgene)
df$Vgene <- factor(df$Vgene,df$Vgene[V.order])
df$Iso <- factor(df$Iso)
df$Iso <- factor(df$Iso, levels=rev(levels(df$Iso)))

all <- df
# bar plot with error bars
p1<-ggplot(data=df, aes(x=Vgene, y=avg))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  #scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('All clones IF - naive_b')+
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
path <- '/home/nofar/Desktop/Lab/data_for_V_plot/Clones/naive_b/'
setwd(path)

files <- list.files(path)
V.usage <- matrix(0, nrow = length(files), ncol = length(V.genes))
for(j in 1:length(files))
{
  df <- read.csv(paste0(path,'/', files[j]))   
  
  #for cdr3
  dfCDR3.peri <- data.frame("CDR3Length" = df$VJ_DIST,"iso" =  rep("naive_B",nrow(df)))
    
  df <- df[df$FUNCTIONAL==T,]
  a <- sort(df$CP_NUM, decreasing = T,  index.return = T)
  df2 <- df[a$ix,]
  df2 <- df2[1:100,]
    
  #for cdr3 - top 100
  dfCDR3Top100.peri <- data.frame("CDR3Length" = df2$VJ_DIST,"iso" =   rep("naive_B",nrow(df2)))
    
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
tmp.df$Iso <-" naive_b"
all.df <- tmp.df

df <- all.df
# order genes by usage sum of all clones
df$Vgene <- factor(df$Vgene,df$Vgene[V.order])
df$Iso <- factor(df$Iso)
df$Iso <- factor(df$Iso, levels=rev(levels(df$Iso)))
# remove gene with low frequencies
#if(!is.null(threshold))
#  df <- df[-which(df$sum<threshold),]


# bar plot with error bars
p2<-ggplot(data=df, aes(x=Vgene, y=avg))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones IF - naive B')+
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


dfCDR3Top100.peri$DataFrom <- "Top 100"
dfCDR3.peri$DataFrom <- "All clones"

CDR3Peri_Merge <- rbind(dfCDR3Top100.peri,dfCDR3.peri)

CDR3repeat <- ddply(CDR3Peri_Merge, c("CDR3Length ","DataFrom"), summarise,
                    repeatTimes = length(CDR3Length))
                  

ggplot(data=CDR3repeat, aes(x=CDR3Length, y=DataFrom, fill = DataFrom))+
  geom_bar( stat="identity",position=position_dodge())+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('CDR3 length Clone size - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=3))





















CDR3Peri_Merge <- rbind(dfCDR3.peri,dfCDR3Top100.peri)
p3 <- ggplot(data=CDR3Peri_Merge, aes(x=CDR3Length, fill=DataFrom))+
  geom_histogram(aes(y=..density..),binwidth=1,colour="black",position=position_dodge())+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('CDR3 length usage IF - naive IGHM')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)
          #axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=3))

grid.arrange(ggplotGrob(p1), ggplotGrob(p2),ggplotGrob(p3))


