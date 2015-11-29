
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

gene.usage.in <- paste0(path, '/data_for_V_plot/Vusage/pfizer/')
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
path <- '/home/nofar/Desktop/Lab/data_for_V_plot/Clones/pfizer/'
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
    
    #for cdr3
    if(i==1)
    {
      dfCDR3.pfizer <- data.frame("CDR3Length" = df$VJ_DIST,"iso" =  dirs[i] )
    }
    else
    {
      tmp <- data.frame("CDR3Length" = df$VJ_DIST,"iso" =  dirs[i] )
      dfCDR3.pfizer <- rbind(dfCDR3.pfizer,tmp)
    }

    df <- df[df$FUNCTIONAL==T,]
    a <- sort(df$CP_NUM, decreasing = T,  index.return = T)
    df2 <- df[a$ix,]
    df2 <- df2[1:100,]
    
    #for cdr3 - top 100
    if(i==1)
    {
      dfCDR3Top100.pfizer <- data.frame("CDR3Length" = df2$VJ_DIST,"iso" =  dirs[i] )
    }
    else
    {
      tmp <- data.frame("CDR3Length" = df2$VJ_DIST,"iso" =  dirs[i] )
      dfCDR3Top100.pfizer <- rbind(dfCDR3Top100.pfizer,tmp)
    }
    
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
#if(!is.null(threshold))
#  df <- df[-which(df$sum<threshold),]

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

p3.1 <- ggplot(data=dfCDR3, aes(x=CDR3Length, fill=iso))+
  geom_histogram(aes(y=..density..),,binwidth=1,colour="black")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('All clones - Pfizer')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='none', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )


p3.2 <- ggplot(data=dfCDR3Top100, aes(x=CDR3Length, fill=iso))+
  geom_histogram(aes(y=..density..),,binwidth=1,colour="black")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('All clones - Pfizer')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='none', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )

dfCDR3Top100.pfizer$DataFrom <- "Top 100"
dfCDR3.pfizer$DataFrom <- "All clones"

CDR3Pfizer_Merge <- rbind(dfCDR3.pfizer,dfCDR3Top100.pfizer)
#CDR3Pfizer_Merge <- CDR3Pfizer_Merge[CDR3Pfizer_Merge$iso =="naive_IGHM",]
p3 <- ggplot(data=CDR3Pfizer_Merge, aes(x=CDR3Length, fill=DataFrom))+
  facet_wrap(~iso) +
  geom_histogram(aes(y=..density..),binwidth=1,colour="black",position=position_dodge())+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('CDR3 length usage - naive IGHM - Pfizer')+
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



                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
# normalize beofre or after removeing bad Vs?????????
#-----------------------------------------------------------------------------------
# Flu
path <-"/home/nofar/Desktop/Lab"
#path <- 'C:/Users/Jennifer/Research/Datasets/Bcells/human/pfizer/analyzed/heavy/'
germline.path <- 'C:/Users/Jennifer/Research/plot_properties/germline/'
germline.path <- path
setwd(path)
organism <- 'human'
chain <- 'VH'

gene.usage.in <- paste0(path, '/data_for_V_plot/Vusage/flu/')
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
#if(!is.null(threshold))
#  df <- df[-which(df$sum<threshold),]

# bar plot with error bars
p4<-ggplot(data=df, aes(x=Vgene, y=avg, fill=Iso))+
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

path <- '/home/nofar/Desktop/Lab/data_for_V_plot/Clones/flu/'
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
    
    #for cdr3
    if(i==1)
    {
      dfCDR3.flu <- data.frame("CDR3Length" = df$VJ_DIST,"iso" =  dirs[i] )
    }
    else
    {
      tmp <- data.frame("CDR3Length" = df$VJ_DIST,"iso" =  dirs[i] )
      dfCDR3.flu <- rbind(dfCDR3.flu,tmp)
    }
    
    
    df <- df[df$FUNCTIONAL==T,]
    a <- sort(df$CP_NUM, decreasing = T,  index.return = T)
    df2 <- df[a$ix,]
    df2 <- df2[1:100,]
    
    #for cdr3 - top 100
    if(i==1)
    {
      dfCDR3Top100.flu <- data.frame("CDR3Length" = df2$VJ_DIST,"iso" =  dirs[i] )
    }
    else
    {
      tmp <- data.frame("CDR3Length" = df2$VJ_DIST,"iso" =  dirs[i] )
      dfCDR3Top100.flu <- rbind(dfCDR3Top100.flu,tmp)
    }
    
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
#if(!is.null(threshold))
#  df <- df[-which(df$sum<threshold),]

# bar plot with error bars
p5<-ggplot(data=df, aes(x=Vgene, y=avg, fill=Iso))+
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

# cdr3 length
p6.1 <- ggplot(data=dfCDR3.flu, aes(x=CDR3Length, fill=iso))+
  geom_histogram(aes(y=..density..),,binwidth=1,colour="black")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('All clones - Flu')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='none', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )


p6.2 <- ggplot(data=dfCDR3Top100.flu, aes(x=CDR3Length, fill=iso))+
  geom_histogram(aes(y=..density..),,binwidth=1,colour="black")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('All clones - Pfizer')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='none', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )


dfCDR3Top100.flu$DataFrom <- "Top 100"
dfCDR3.flu$DataFrom <- "All clones"

CDR3Flu_Merge <- rbind(dfCDR3.flu,dfCDR3Top100.flu)
write.csv(CDR3Flu_Merge ,file = "flu_cdrDetails.csv")
write.csv(CDR3Pfizer_Merge ,file = "pfizer_cdrDetails.csv")
CDR3Flu_Merge <- CDR3Flu_Merge[CDR3Flu_Merge$iso =="IGHM",]
p6 <- ggplot(data=CDR3Flu_Merge, aes(x=CDR3Length, fill=DataFrom))+
  facet_wrap(~iso) +
  geom_histogram(aes(y=..density..),binwidth=1,colour="black",position=position_dodge())+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('CDR3 length usage- IGHM - Flu')+
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
  guides(fill=guide_legend(nrow=2))





# keck 

setwd("/home/nofar/Desktop/Lab")
Keck98.top500.df <- read.csv(file = 'Keck98.top500EditVnames.csv' ,header= T, stringsAsFactors= F)
Keck98.top500.df$frequencyNormalized <- Keck98.top500.df$frequencyNormalized/100


# get all onframe data
IF_df_all500 <- Keck98.top500.df[Keck98.top500.df$sequenceStatus=="In",]

#sample_num_in_out <-ddply(IF_df_all500, c("sample"), summarise,
#                         IF = length(sample),
#                         OF = 500 - IF)

#Keck98.top500.df$IF_Num_to_sample
#Keck98.top500.df$OF_Num_to_sample
#for(i in 1:nrow(Keck98.top500.df))
#{
#  sampleName <-  Keck98.top500.df[i,"sample"]
#  Keck98.top500.df[i,"IF_Num_to_sample"] <-sample_num_in_out[sample_num_in_out$sample==sampleName,"IF"]
#  Keck98.top500.df[i,"OF_Num_to_sample"] <-sample_num_in_out[sample_num_in_out$sample==sampleName,"OF"]
#}
#write.csv(Keck98.top500.df,file = 'Keck98.top500EditVnames.csv')

IF_df_all500 <- Keck98.top500.df[Keck98.top500.df$sequenceStatus=="In",]


IF_Vrepeat_all500_by_sample <- ddply(IF_df_all500, c("EditV","sequenceStatus","sample"), summarise,
                                     repeatTimes = length(EditV),
                                     #  frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"IF"],
                                    frecuency = repeatTimes/unique(IF_Num_to_sample),
                                    Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
                                    )

IF_Vrepeat_all500 <- ddply(IF_Vrepeat_all500_by_sample, c("EditV","sequenceStatus"), summarise,
                           repeatTimes = sum(repeatTimes),
                           meanFrecuency = mean(frecuency),
                           SDfrecuency = sd(frecuency),
                           SEfrecuency = SDfrecuency/sqrt(98),
                           meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                           SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                           SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
                          )

IF_Vrepeat_all500$frecuency <- IF_Vrepeat_all500$repeatTimes/sum(IF_Vrepeat_all500$repeatTimes)

IF_Vrepeat_all500 <- arrange(IF_Vrepeat_all500, -frecuency)
IF_Vrepeat_all500$EditV <- factor(as.character(IF_Vrepeat_all500$EditV),levels=as.character(IF_Vrepeat_all500$EditV))

# plot - only first 40
  p7<- ggplot(data=IF_Vrepeat_all500[1:40,], aes(x=EditV, y=frecuency ))+
  geom_bar( position="stack", stat="identity", fill="gray")+
   geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 500 clones - Keck')+
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

# keck - top 100 clones
sample.vector <- unique(IF_df_all500$sample)
df <-IF_df_all500[IF_df_all500$sample==sample.vector[1],]
df <-arrange(df,-frequencyNormalized)
dfFirst100 <- head(df,100)

for(i in 2:length(sample.vector))
{
  df <-IF_df_all500[IF_df_all500$sample==sample.vector[i],]
  df <-arrange(df,-frequencyNormalized)
  dfFirst100 <-  rbind(dfFirst100,head(df,100))
}

IF_Vrepeat_first500_by_sample <- ddply(dfFirst100, c("EditV","sequenceStatus","sample"), summarise,
                                     repeatTimes = length(EditV),
                                     #  frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"IF"],
                                     frecuency = repeatTimes/unique(IF_Num_to_sample),
                                     Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
)

IF_Vrepeat_first500 <- ddply(IF_Vrepeat_first500_by_sample, c("EditV","sequenceStatus"), summarise,
                           repeatTimes = sum(repeatTimes),
                           meanFrecuency = mean(frecuency),
                           SDfrecuency = sd(frecuency),
                           SEfrecuency = SDfrecuency/sqrt(98),
                           meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                           SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                           SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
)

IF_Vrepeat_first500$frecuency <- IF_Vrepeat_first500$repeatTimes/sum(IF_Vrepeat_first500$repeatTimes)

IF_Vrepeat_first500_40 <- arrange(IF_Vrepeat_first500, -frecuency)
IF_Vrepeat_first500_40 <- IF_Vrepeat_first500_40[1:40,]
IF_Vrepeat_first500_40$EditV <- factor(as.character(IF_Vrepeat_first500_40$EditV),levels=as.character(IF_Vrepeat_all500$EditV))

# plot - only first 40 
p8<-ggplot(data=IF_Vrepeat_first500_40, aes(x=EditV, y=frecuency))+
    geom_bar( position="stack", stat="identity",fill="gray")+
    geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
    # scale_fill_manual(values=cols) +
    theme_bw() +
    ggtitle('Top 100 clones - Keck')+
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

grid.arrange(ggplotGrob(p7), ggplotGrob(p8), ncol=2)

# cdr3 length

IN_CDR3repeat_by_sample <- ddply(IF_df_all500, c("EditCDR3Length","sequenceStatus","sample"), summarise,
                                       repeatTimes = length(EditCDR3Length),
                                 #  frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"IF"],
                                 frecuency = repeatTimes/unique(IF_Num_to_sample),
                                       Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
)

IN_CDR3repeat <- ddply(IN_CDR3repeat_by_sample, c("EditCDR3Length","sequenceStatus"), summarise,
                             repeatTimes = sum(repeatTimes),
                             meanFrecuency = mean(frecuency),
                             SDfrecuency = sd(frecuency),
                             SEfrecuency = SDfrecuency/sqrt(98),
                             meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                             SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                             SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
)

IN_CDR3repeat$frecuency <- IN_CDR3repeat$repeatTimes/sum(IN_CDR3repeat$repeatTimes)

p9<-ggplot(data=IN_CDR3repeat, aes(x=EditCDR3Length, y=frecuency))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 500 clones - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )


p18 <-ggplot(data=CDR3keck_Merge, aes(x=EditCDR3Length, y=frecuency, fill=DataFrom))+
  geom_bar( stat="identity",position=position_dodge())+
  
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('CDR3 length usage - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4,position=position_dodge(.9)) +
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)
          #axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=3))

#########################################################################################################################

grid.arrange(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4),ggplotGrob(p5),ggplotGrob(p6),ggplotGrob(p7),ggplotGrob(p8),ggplotGrob(p18), ncol=3)

#########################################################################################################################

# V_usage naive_B_cells data;(no error bar becuase it only one sample)

setwd("/home/nofar/Desktop/Lab")
Edit_naive_B_cells.df <- read.csv(file = 'Edit_naive_B_cells.csv' ,header= T, stringsAsFactors= F)

# get all inframe data
IF_df_naive_B <- Edit_naive_B_cells.df[Edit_naive_B_cells.df$sequenceStatus=="Productive",]

IF_Vrepeat_naive_B <- ddply(IF_df_naive_B, c("EditV","sequenceStatus"), summarise,
                           repeatTimes = length(EditV)
                  
                            )

s <-sum(IF_Vrepeat_naive_B$repeatTimes)
IF_Vrepeat_naive_B$frecuency <- IF_Vrepeat_naive_B$repeatTimes / s

IF_Vrepeat_naive_B <- arrange(IF_Vrepeat_naive_B, -frecuency)
IF_Vrepeat_naive_B$EditV <- factor(as.character(IF_Vrepeat_naive_B$EditV),levels=as.character(IF_Vrepeat_naive_B$EditV))


x1 <-ggplot(data=IF_Vrepeat_naive_B, aes(x=EditV, y=frecuency))+
  geom_bar( position="stack", stat="identity")+
  #geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('IF  - naive_B_cells')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )


# get all outframe data

OF_df_naive_B <- Edit_naive_B_cells.df[Edit_naive_B_cells.df$sequenceStatus!="Productive",]
OF_Vrepeat_naive_B <- ddply(OF_df_naive_B, c("EditV","sequenceStatus"), summarise,
                            repeatTimes = length(EditV)                  
)

s <-sum(OF_Vrepeat_naive_B$repeatTimes)
OF_Vrepeat_naive_B$frecuency <- OF_Vrepeat_naive_B$repeatTimes / s

OF_Vrepeat_naive_B <- arrange(OF_Vrepeat_naive_B, -frecuency)
OF_Vrepeat_naive_B$EditV <- factor(as.character(OF_Vrepeat_naive_B$EditV),levels=as.character(IF_Vrepeat_naive_B$EditV))

x2 <-ggplot(data=OF_Vrepeat_naive_B, aes(x=EditV, y=frecuency))+
  geom_bar( position="stack", stat="identity")+
  #geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('OF  - naive_B_cells')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )

grid.arrange(ggplotGrob(x1), ggplotGrob(x2))






#~~~~~~~ print V_usage keck 500 in/100 in / out

# v usage - outfrmae/stop

OF_df_all <- Keck98.top500.df[Keck98.top500.df$sequenceStatus!="In",]


OF_Vrepeat_all_by_sample <- ddply(OF_df_all, c("EditV","sequenceStatus","sample"), summarise,
                                     repeatTimes = length(EditV),
                                    # frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"OF"],
                                     frecuency = repeatTimes/unique(OF_Num_to_sample),
                                     Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
)

OF_Vrepeat_all <- ddply(OF_Vrepeat_all_by_sample, c("EditV","sequenceStatus"), summarise,
                           repeatTimes = sum(repeatTimes),
                           meanFrecuency = mean(frecuency),
                           SDfrecuency = sd(frecuency),
                           SEfrecuency = SDfrecuency/sqrt(98),
                           meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                           SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                           SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
                      )

OF_Vrepeat_all <- ddply(OF_Vrepeat_all_by_sample, c("EditV"), summarise,
                        repeatTimes = sum(repeatTimes),
                        meanFrecuency = mean(frecuency),
                        SDfrecuency = sd(frecuency),
                        SEfrecuency = SDfrecuency/sqrt(98),
                        meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                        SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                        SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
)
OF_Vrepeat_all$frecuency <- OF_Vrepeat_all$repeatTimes/sum(OF_Vrepeat_all$repeatTimes)

to_remove <- c()
index <- 1
for(i in 1:nrow(OF_Vrepeat_all))
{
  if(OF_Vrepeat_all[i,"EditV"]%in%IF_Vrepeat_all500[1:40,"EditV"] ==F)
  {
    to_remove[index]=i
    index <- 1 + index
  }

}
OF_Vrepeat_all_40 <- OF_Vrepeat_all[-c(to_remove),]

OF_Vrepeat_all_40$EditV <- factor(as.character(OF_Vrepeat_all_40$EditV),levels=as.character(IF_Vrepeat_all500$EditV))

# plot - only first 40
p10<- ggplot(data=OF_Vrepeat_all_40, aes(x=EditV, y=frecuency))+#, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('OF - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  theme(  legend.title=element_blank(),
          #legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=2))

grid.arrange(ggplotGrob(p7), ggplotGrob(p8),ggplotGrob(p10))

#~~~~~~~ print J_usage keck 500 in/100 in / out

# j usage - 500 in  
IF_Jrepeat_all500_by_sample <- ddply(IF_df_all500, c("EditJ","sequenceStatus","sample"), summarise,
                                     repeatTimes = length(EditJ),
                                     #  frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"IF"],
                                     frecuency = repeatTimes/unique(IF_Num_to_sample),
                                     Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
                                    )

IF_Jrepeat_all500 <- ddply(IF_Jrepeat_all500_by_sample, c("EditJ","sequenceStatus"), summarise,
                        repeatTimes = sum(repeatTimes),
                        meanFrecuency = mean(frecuency),
                        SDfrecuency = sd(frecuency),
                        SEfrecuency = SDfrecuency/sqrt(98),
                        meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                        SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                        SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
)

IF_Jrepeat_all500$frecuency <- IF_Jrepeat_all500$repeatTimes/sum(IF_Jrepeat_all500$repeatTimes)

IF_Jrepeat_all500 <- arrange(IF_Jrepeat_all500, -frecuency)
IF_Jrepeat_all500$EditJ <- factor(as.character(IF_Jrepeat_all500$EditJ),levels=as.character(IF_Jrepeat_all500$EditJ))


# plot -
p11<- ggplot(data=IF_Jrepeat_all500, aes(x=EditJ, y=frecuency))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 500 clones - Keck')+
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

# keck - j usage - top 100 in
sample.vector <- unique(IF_df_all500$sample)
df <-IF_df_all500[IF_df_all500$sample==sample.vector[1],]
df <-arrange(df,-frequencyNormalized)
dfFirst100 <- head(df,100)

for(i in 2:length(sample.vector))
{
  df <-IF_df_all500[IF_df_all500$sample==sample.vector[i],]
  df <-arrange(df,-frequencyNormalized)
  dfFirst100 <-  rbind(dfFirst100,head(df,100))
}

IF_Jrepeat_first500_by_sample <- ddply(dfFirst100, c("EditJ","sequenceStatus","sample"), summarise,
                                     repeatTimes = length(EditJ),
                                     #  frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"IF"],
                                     frecuency = repeatTimes/unique(IF_Num_to_sample),
                                     Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
)

IF_Jrepeat_first500 <- ddply(IF_Jrepeat_first500_by_sample, c("EditJ","sequenceStatus"), summarise,
                           repeatTimes = sum(repeatTimes),
                           meanFrecuency = mean(frecuency),
                           SDfrecuency = sd(frecuency),
                           SEfrecuency = SDfrecuency/sqrt(98),
                           meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                           SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                           SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
)

IF_Jrepeat_first500$frecuency <- IF_Jrepeat_first500$repeatTimes/sum(IF_Jrepeat_first500$repeatTimes)

IF_Jrepeat_first500$EditJ <- factor(as.character(IF_Jrepeat_first500$EditJ),levels=as.character(IF_Jrepeat_all500$EditJ))

# plot - 
p12<-ggplot(data=IF_Jrepeat_first500, aes(x=EditJ, y=frecuency))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - J - Keck')+
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

# keck - j usage - all out

OF_Jrepeat_all_by_sample <- ddply(OF_df_all, c("EditJ","sequenceStatus","sample"), summarise,
                                       repeatTimes = length(EditJ),
                                       #  frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"IF"],
                                       frecuency = repeatTimes/unique(OF_Num_to_sample),
                                       Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
)

OF_Jrepeat_all <- ddply(OF_Jrepeat_all_by_sample, c("EditJ"), summarise,
                             repeatTimes = sum(repeatTimes),
                             meanFrecuency = mean(frecuency),
                             SDfrecuency = sd(frecuency),
                             SEfrecuency = SDfrecuency/sqrt(98),
                             meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                             SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                             SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
)

OF_Jrepeat_all$frecuency <- OF_Jrepeat_all$repeatTimes/sum(OF_Jrepeat_all$repeatTimes)

OF_Jrepeat_all$EditJ <- factor(as.character(OF_Jrepeat_all$EditJ),levels=as.character(IF_Jrepeat_all500$EditJ))
OF_Jrepeat_all <- OF_Jrepeat_all[!is.na(OF_Jrepeat_all$EditJ),]  

# plot - 
p13<- ggplot(data=OF_Jrepeat_all, aes(x=EditJ, y=frecuency))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('OF - Keck')+
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

grid.arrange(ggplotGrob(p11), ggplotGrob(p12), ggplotGrob(p13))


#~~~~~~~ print CDR3_usage keck 500 in/100 in / out

# keck - cdr3 length usage - top 500 in <- p9
# keck - cdr3 length usage - top 100 in 

sample.vector <- unique(IF_df_all500$sample)
df <-IF_df_all500[IF_df_all500$sample==sample.vector[1],]
df <-arrange(df,-frequencyNormalized)
dfFirst100 <- head(df,100)

for(i in 2:length(sample.vector))
{
  df <-IF_df_all500[IF_df_all500$sample==sample.vector[i],]
  df <-arrange(df,-frequencyNormalized)
  dfFirst100 <-  rbind(dfFirst100,head(df,100))
}

IN_CDR3repeat_first100_by_sample <- ddply(dfFirst100, c("EditCDR3Length","sequenceStatus","sample"), summarise,
                                 repeatTimes = length(EditCDR3Length),
                                 #  frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"IF"],
                                 frecuency = repeatTimes/unique(IF_Num_to_sample),
                                 Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
)

IN_CDR3repeat_first100 <- ddply(IN_CDR3repeat_first100_by_sample, c("EditCDR3Length","sequenceStatus"), summarise,
                       repeatTimes = sum(repeatTimes),
                       meanFrecuency = mean(frecuency),
                       SDfrecuency = sd(frecuency),
                       SEfrecuency = SDfrecuency/sqrt(98),
                       meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                       SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                       SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
)

IN_CDR3repeat_first100$frecuency <- IN_CDR3repeat_first100$repeatTimes/sum(IN_CDR3repeat_first100$repeatTimes)

p14<-ggplot(data=IN_CDR3repeat_first100, aes(x=EditCDR3Length, y=frecuency))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )

# keck - cdr3 length usage - of 

OF_df_all <- Keck98.top500.df[Keck98.top500.df$sequenceStatus!="In",]

OF_CDR3repeat_all_by_sample <- ddply(OF_df_all, c("EditCDR3Length","sequenceStatus","sample"), summarise,
                                          repeatTimes = length(EditCDR3Length),
                                          #  frecuency = repeatTimes/sample_num_in_out[sample_num_in_out$sample==sample,"IF"],
                                          frecuency = repeatTimes/unique(OF_Num_to_sample),
                                          Avg.frequencyNormalized = exp(mean(log(frequencyNormalized)))
)

OF_CDR3repeat_all <- ddply(OF_CDR3repeat_all_by_sample, c("EditCDR3Length"), summarise,
                                repeatTimes = sum(repeatTimes),
                                meanFrecuency = mean(frecuency),
                                SDfrecuency = sd(frecuency),
                                SEfrecuency = SDfrecuency/sqrt(98),
                                meanOfAvg.frequencyNormalized = exp(mean(log(Avg.frequencyNormalized))),
                                SDofAvg.frequencyNormalized = sd(meanOfAvg.frequencyNormalized),
                                SEofAvg.frequencyNormalized = SDofAvg.frequencyNormalized/sqrt(98)
)

OF_CDR3repeat_all$frecuency <- OF_CDR3repeat_all$repeatTimes/sum(OF_CDR3repeat_all$repeatTimes)


# plot - 
p15<- ggplot(data=OF_CDR3repeat_all, aes(x=EditCDR3Length, y=frecuency))+
  geom_bar( position="stack", stat="identity")+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('OF - Keck')+
  scale_x_continuous(limits = c(0, 50))+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  theme(  legend.title=element_blank(),
          #legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=2))

grid.arrange(ggplotGrob(p9), ggplotGrob(p14), ggplotGrob(p15))

######################################################################### merge into one plot - Vusage keck

#OF_Vrepeat_all_40
#IF_Vrepeat_all500[1:40,]
#IF_Vrepeat_first500_40

OF_V.df <- OF_Vrepeat_all_40
OF_V.df$DataFrom <- "OF"
OF_V.df$sequenceStatus <- "OF"
IN_500_V.df <- IF_Vrepeat_all500[1:40,]
IN_500_V.df$DataFrom <- "IF - All"
IN_100_V.df <- IF_Vrepeat_first500_40
IN_100_V.df$DataFrom <- "IF - Top 100"


Vkeck_Merge <- rbind(IN_500_V.df,IN_100_V.df,OF_V.df)
p16 <-ggplot(data=Vkeck_Merge, aes(x=EditV, y=frecuency, fill=DataFrom))+
  geom_bar( stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4,position=position_dodge(.9)) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Vusage - Keck')+
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
  guides(fill=guide_legend(nrow=3))

######################################################################### merge into one plot - Jusage keck

#IF_Jrepeat_first500
#IF_Jrepeat_all500
#OF_Jrepeat_all

OF_J.df <- OF_Jrepeat_all

OF_J.df$DataFrom <- "OF"
OF_J.df$sequenceStatus <- "OF"
IN_500_J.df <- IF_Jrepeat_all500[1:40,]
IN_500_J.df$DataFrom <- "IF - All"
IN_100_J.df <- IF_Jrepeat_first500
IN_100_J.df$DataFrom <- "IF - Top 100"


Jkeck_Merge <- rbind(IN_500_J.df,IN_100_J.df,OF_J.df)
Jkeck_Merge <- Jkeck_Merge[!is.na(Jkeck_Merge$EditJ),]
p17 <-ggplot(data=Jkeck_Merge, aes(x=EditJ, y=frecuency, fill=DataFrom))+
  geom_bar( stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4,position=position_dodge(.9)) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Jusage - Keck')+
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
  guides(fill=guide_legend(nrow=3))

######################################################################### merge into one plot - CDR3 Length - keck
#IN_CDR3repeat
#OF_CDR3repeat_all
#IN_CDR3repeat_first100

OF_CDR3.df <- OF_CDR3repeat_all
OF_CDR3.df$DataFrom <- "OF"
OF_CDR3.df$sequenceStatus <- "OF"
IN_500_CDR3.df <- IN_CDR3repeat
IN_500_CDR3.df$DataFrom <- "IF - All"
IN_100_CDR3.df <- IN_CDR3repeat_first100
IN_100_CDR3.df$DataFrom <- "IF - Top 100"


CDR3keck_Merge <- rbind(IN_500_CDR3.df,IN_100_CDR3.df,OF_CDR3.df)
p18 <-ggplot(data=CDR3keck_Merge, aes(x=EditCDR3Length, y=frecuency, fill=DataFrom))+
  geom_bar( stat="identity",position=position_dodge())+
  
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('CDR3 length usage - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  scale_x_continuous(limits = c(0, 50))+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4,position=position_dodge(.9)) +
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=3))


ggplot(data=CDR3keck_Merge, aes(x=EditCDR3Length, fill=DataFrom))+
         geom_histogram(aes(y=..density..),binwidth=1,colour="black",position=position_dodge())+
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('CDR3 length usage - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Frequency")+
  scale_x_continuous(limits = c(0, 50))+
  geom_errorbar(aes(ymin=frecuency-SEfrecuency, ymax=frecuency+SEfrecuency), width=.4,position=position_dodge(.9)) +
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=3))

ggplot(data=CDR3Pfizer_Merge, aes(x=CDR3Length, fill=DataFrom))+
  geom_histogram(aes(y=..density..),binwidth=1,colour="black",position=position_dodge())+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('CDR3 length usage - naive IGHM - Pfizer')+
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

grid.arrange(ggplotGrob(p16), ggplotGrob(p17), ggplotGrob(p18))

################################ merge into one plot - v - CLONE SIZE - keck

# IF_Vrepeat_all500
# ~~ all IN
IF_Vrepeat_all500_CloneSize <- arrange(IF_Vrepeat_all500, -meanOfAvg.frequencyNormalized)
IF_Vrepeat_all500_CloneSize$EditV <- factor(as.character(IF_Vrepeat_all500_CloneSize$EditV),
                                            levels=as.character(IF_Vrepeat_all500_CloneSize$EditV))

x1 <-ggplot(data=IF_Vrepeat_all500_CloneSize[1:40,], aes(x=EditV, y=meanOfAvg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 500 clones - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )

# ~~ Top 100 IN
#IF_Vrepeat_first500

#IF_Vrepeat_first100_CloneSize <- arrange(IF_Vrepeat_first500, -Avg.frequencyNormalized)
IF_Vrepeat_first100_CloneSize <- arrange(IF_Vrepeat_first500, -meanOfAvg.frequencyNormalized)
IF_Vrepeat_first100_CloneSize <- IF_Vrepeat_first100_CloneSize[1:40,]
IF_Vrepeat_first100_CloneSize$EditV <- factor(as.character(IF_Vrepeat_first100_CloneSize$EditV),
                                              levels=as.character(IF_Vrepeat_all500_CloneSize$EditV))

x2 <- ggplot(data=IF_Vrepeat_first100_CloneSize, aes(x=EditV, y=meanOfAvg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )

grid.arrange(ggplotGrob(x1), ggplotGrob(x2))
#~~~~~~~~ OF 
#OF_Vrepeat_all

to_remove <- c()
index <- 1
for(i in 1:nrow(OF_Vrepeat_all))
{
  if(OF_Vrepeat_all[i,"EditV"]%in%IF_Vrepeat_all500_CloneSize[1:40,"EditV"] ==F)
  {
    to_remove[index]=i
    index <- 1 + index
  }
  
}

OF_Vrepeat_all_CloneSize_40 <- OF_Vrepeat_all[-c(to_remove),]
OF_Vrepeat_all_CloneSize_40$EditV <- factor(as.character(OF_Vrepeat_all_CloneSize_40$EditV),
                                            levels=as.character(IF_Vrepeat_all500_CloneSize$EditV))

ggplot(data=OF_Vrepeat_all_CloneSize_40, aes(x=EditV, y=meanOfAvg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )

#IF_Vrepeat_all500_CloneSize
#IF_Vrepeat_first100_CloneSize
#OF_Vrepeat_all_CloneSize_40


IF_500_cloneSize_V.df <- IF_Vrepeat_all500_CloneSize[1:40,]
IF_500_cloneSize_V.df$DataFrom <- "IF - All"
IF_100_cloneSize_V.df <- IF_Vrepeat_first100_CloneSize
IF_100_cloneSize_V.df$DataFrom <- "IF - Top 100"
OF_cloneSize_v.df <- OF_Vrepeat_all_CloneSize_40
OF_cloneSize_v.df$DataFrom <- "OF"


IF_500_cloneSize_V.df <- data.frame(Vgene = IF_500_cloneSize_V.df$EditV,Avg.frequencyNormalized =IF_500_cloneSize_V.df$meanOfAvg.frequencyNormalized,
                                    DataFrom = IF_500_cloneSize_V.df$DataFrom)
IF_100_cloneSize_V.df <- data.frame(Vgene = IF_100_cloneSize_V.df$EditV,Avg.frequencyNormalized =IF_100_cloneSize_V.df$meanOfAvg.frequencyNormalized,
                                    DataFrom = IF_100_cloneSize_V.df$DataFrom)
OF_cloneSize_v.df <- data.frame(Vgene = OF_cloneSize_v.df$EditV,Avg.frequencyNormalized =OF_cloneSize_v.df$meanOfAvg.frequencyNormalized,
                                    DataFrom = OF_cloneSize_v.df$DataFrom)

# add missing vGene into the 100 in and to of:
#top 100
vGeneToAdd100 <-IF_500_cloneSize_V.df[IF_500_cloneSize_V.df$Vgene%in%IF_100_cloneSize_V.df$Vgene==F,"Vgene"]
IF_100_cloneSize_V.df <- rbind(IF_100_cloneSize_V.df,data.frame(Vgene = vGeneToAdd100,Avg.frequencyNormalized = rep(0, length(vGeneToAdd100)),
           DataFrom = rep("IF - Top 100", length(vGeneToAdd100))))
IF_100_cloneSize_V.df$Vgene <-factor(as.character(IF_100_cloneSize_V.df$Vgene),
                                     levels=as.character(IF_500_cloneSize_V.df$Vgene))
IF_100_cloneSize_V.df <- IF_100_cloneSize_V.df[!is.na(IF_100_cloneSize_V.df$Vgene),]

#of 
vGeneToOF <-IF_500_cloneSize_V.df[IF_500_cloneSize_V.df$Vgene%in%OF_cloneSize_v.df$Vgene==F,"Vgene"]
OF_cloneSize_v.df <-rbind(OF_cloneSize_v.df,data.frame(Vgene = vGeneToOF,Avg.frequencyNormalized = rep(0, length(vGeneToOF)),
                                                           DataFrom = rep("OF", length(vGeneToOF))))
OF_cloneSize_v.df$Vgene <-factor(as.character(OF_cloneSize_v.df$Vgene),
                                     levels=as.character(IF_500_cloneSize_V.df$Vgene))
OF_cloneSize_v.df <- OF_cloneSize_v.df[!is.na(OF_cloneSize_v.df$Vgene),]

Vkeck_cloneSize_Merge <- rbind(IF_500_cloneSize_V.df,IF_100_cloneSize_V.df,OF_cloneSize_v.df)
p19 <-ggplot(data=Vkeck_cloneSize_Merge, aes(x=Vgene, y=Avg.frequencyNormalized, fill=DataFrom))+
  geom_bar( stat="identity",position=position_dodge())+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Vgene Clone size - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=3))


################################ merge into one plot - J - CLONE SIZE - keck

# IF_Jrepeat_all500
# ~~ all IN
IF_Jrepeat_all500_CloneSize <- arrange(IF_Jrepeat_all500, -meanOfAvg.frequencyNormalized)
IF_Jrepeat_all500_CloneSize$EditJ <- factor(as.character(IF_Jrepeat_all500_CloneSize$EditJ),
                                            levels=as.character(IF_Jrepeat_all500_CloneSize$EditJ))

ggplot(data=IF_Jrepeat_all500_CloneSize, aes(x=EditJ, y=meanOfAvg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 500 clones - J - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )


# ~~ Top 100 IN
#IF_Jrepeat_first500

IF_Jrepeat_first100_CloneSize<-IF_Jrepeat_first500
IF_Jrepeat_first100_CloneSize$EditJ <- factor(as.character(IF_Jrepeat_first100_CloneSize$EditJ),
                                              levels=as.character(IF_Jrepeat_all500_CloneSize$EditJ))

ggplot(data=IF_Jrepeat_first100_CloneSize, aes(x=EditJ, y=meanOfAvg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - J - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )


#~~~~~~~~ OF 
#OF_Jrepeat_all

OF_Jrepeat_CloneSize <-OF_Jrepeat_all[!is.na(OF_Jrepeat_all$EditJ),]
OF_Jrepeat_CloneSize$EditJ <- factor(as.character(OF_Jrepeat_CloneSize$EditJ),
                                     levels=as.character(IF_Jrepeat_all500$EditJ))


ggplot(data=OF_Jrepeat_CloneSize, aes(x=EditJ, y=meanOfAvg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - J - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )

#IF_Jrepeat_all500_CloneSize
#IF_Jrepeat_first100_CloneSize
#OF_Jrepeat_CloneSize

IF_500_cloneSize_J.df <- IF_Jrepeat_all500_CloneSize
IF_500_cloneSize_J.df$DataFrom <- "IF - All"
IF_100_cloneSize_J.df <- IF_Jrepeat_first100_CloneSize
IF_100_cloneSize_J.df$DataFrom <- "IF - Top 100"
OF_cloneSize_J.df <- OF_Jrepeat_CloneSize
OF_cloneSize_J.df$sequenceStatus <-"OF"
OF_cloneSize_J.df$DataFrom <- "OF"


# add missing jGene into the 100 in and to of:
#top 100
jGeneToAdd100 <-IF_500_cloneSize_J.df[IF_500_cloneSize_J.df$EditJ%in%IF_100_cloneSize_J.df$EditJ==F,]
jGeneToAdd100$meanOfAvg.frequencyNormalized <- 0
jGeneToAdd100$DataFrom <- "IF - Top 100"

OF_cloneSize_J.df <- rbind(OF_cloneSize_J.df,jGeneToAdd100)

#of
jGeneToAddOF <-IF_500_cloneSize_J.df[IF_500_cloneSize_J.df$EditJ%in%OF_cloneSize_J.df$EditJ==F,]
jGeneToAddOF$meanOfAvg.frequencyNormalized <- 0
jGeneToAddOF$DataFrom <- "OF"

OF_cloneSize_J.df <- rbind(OF_cloneSize_J.df,jGeneToAddOF)


Jkeck_cloneSize_Merge <- rbind(IF_500_cloneSize_J.df,IF_100_cloneSize_J.df,OF_cloneSize_J.df)
p20 <-ggplot(data=Jkeck_cloneSize_Merge, aes(x=EditJ, y=meanOfAvg.frequencyNormalized, fill=DataFrom))+
  geom_bar( stat="identity",position=position_dodge())+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Jgene Clone size - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=3))

################################ merge into one plot - CDR3 - CLONE SIZE - keck

#IN_CDR3repeat
# ~~ all IN

IF_JCDR3repeat_all_CloneSize <-IN_CDR3repeat

ggplot(data=IF_JCDR3repeat_all_CloneSize, aes(x=EditCDR3Length, y=Avg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 500 clones - CDR3 length - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )

#~~~~~~~~ Top 100
#IN_CDR3repeat_first100
IF_JCDR3repeat_first100_CloneSize <- IN_CDR3repeat_first100

ggplot(data=IF_JCDR3repeat_all_CloneSize, aes(x=EditCDR3Length, y=meanOFAvg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('Top 100 clones - CDR3 length - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )

#~~~~~~~~ OF 
#OF_CDR3repeat_all

OF_JCDR3repeat_all_CloneSize <- OF_CDR3repeat_all

ggplot(data=OF_JCDR3repeat_all_CloneSize, aes(x=EditCDR3Length, y=meanOFAvg.frequencyNormalized, fill=sequenceStatus))+
  geom_bar( position="stack", stat="identity")+
  # geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.4) +
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('OF clones - CDR3 length - Keck')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  ylab("Avg frequency normalized")+
  xlab("CDR3 length")+
  scale_x_continuous(limits = c(0, 50))+
  theme(  legend.title=element_blank(),
          legend.position='non', #c(1,1), # legend position
          # legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(angle=90, vjust=1)#,
          #axis.title.x=element_blank()
  )


#IF_JCDR3repeat_all_CloneSize
#IF_JCDR3repeat_first100_CloneSize
#OF_JCDR3repeat_all_CloneSize

IF_500_cloneSize_CDR3.df <- IF_JCDR3repeat_all_CloneSize
IF_500_cloneSize_CDR3.df$DataFrom <- "IF - All"
IF_100_cloneSize_CDR3.df <- IF_JCDR3repeat_first100_CloneSize
IF_100_cloneSize_CDR3.df$DataFrom <- "IF - Top 100"
OF_cloneSize_CDR3.df <- OF_JCDR3repeat_all_CloneSize
OF_cloneSize_CDR3.df$DataFrom <- "OF"
OF_cloneSize_CDR3.df$sequenceStatus <- "OF"

CDR3keck_cloneSize_Merge <- rbind(IF_500_cloneSize_CDR3.df,IF_100_cloneSize_CDR3.df,OF_cloneSize_CDR3.df)
p21 <-ggplot(data=CDR3keck_cloneSize_Merge, aes(x=EditCDR3Length, y=meanOfAvg.frequencyNormalized, fill=DataFrom))+
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


grid.arrange(ggplotGrob(p19), ggplotGrob(p20), ggplotGrob(p21))

grid.arrange(ggplotGrob(p16), ggplotGrob(p17), ggplotGrob(p18),ggplotGrob(p19), ggplotGrob(p20), ggplotGrob(p21),ncol=3)


