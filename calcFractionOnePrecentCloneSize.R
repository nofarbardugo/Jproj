
library(ggplot2)
library(reshape2)
library(seqinr)
require(gtable)
library(gridExtra)
require(plyr)
require(beeswarm)


# Pfizer

path <- '/home/nofar/Desktop/Lab/data_for_V_plot/Clones/pfizer/'
setwd(path)

cloneData.df = data.frame(Iso = c(0),sample = c(0),cloneSize = c(0))

isotype.dirs <- list.files(path)

for(i in 1:length(isotype.dirs))
{
  sample.files <- list.files(isotype.dirs[i])
  # remove bad donors
  sample.files <- sample.files[-grep('136', sample.files)]
  sample.files <- sample.files[-grep('273', sample.files)] 
  
   for(j in 1:length(sample.files))
   {
     df <- read.csv(paste0( isotype.dirs[i], '/', sample.files[j]),header= T, stringsAsFactors= F)
     lPlacesCloneId <- gregexpr(pattern ='_',df[1,"CLONE_ID"]) # get places of the char '_' in "CLONE_ID" column
     patientNumber <- substr(df[1,"CLONE_ID"],start = lPlacesCloneId[[1]][5] +1,stop = lPlacesCloneId[[1]][6]-1)
  
     #for pfizer - remove last char
     
     if(patientNumber== "91B" || patientNumber== "91A")
     {
       patientNumber <- substr(patientNumber, 1,2)
     }else
     {
       patientNumber <- substr(patientNumber, 1,3)
     }
     
     # check IF/OF
   #  is.IF <- c()
    rowNum <- nrow(df)
  #   for(n in 1:rowNum)
   #  {
  #     if(df[n,"INFRAME"]==T && df[n,"STOP"]==F)
  #     {
  #       is.IF[n]<-"IF"
  #     }
  #     else
  #     {
  #       is.IF[n]<-"OF"
  #     }
  #   }
     
     if(i==1)
     {
       df$INFRAME
       
       cloneData.df <- data.frame(sample = rep(patientNumber,rowNum),cloneSize = df$CP_NUM, sequenceStatus=  df$FUNCTIONAL,Iso = rep(isotype.dirs[i],rowNum))
     }else
     {
       cloneData.df <- rbind(cloneData.df,
                             data.frame(sample = rep(patientNumber,rowNum),cloneSize = df$CP_NUM,sequenceStatus=  df$FUNCTIONA,Iso = rep(isotype.dirs[i],rowNum)))
     }
     
   }
}

# save the pfizer data - 
setwd("/home/nofar/Desktop/Lab")
write.csv(cloneData.df,file = 'pfizerCloneData.csv')


# Flu

path <- '/home/nofar/Desktop/Lab/data_for_V_plot/Clones/flu/'
setwd(path)

cloneData_flu.df = data.frame(Iso = c(0),sample = c(0),cloneSize = c(0))
isotype.dirs <- list.files(path)
isotype.dirs <- isotype.dirs[-3] # remove IGHE
for(i in 1:length(isotype.dirs))
{
  sample.files <- list.files(isotype.dirs[i])
  # remove bad donors
  sample.files <- sample.files[-grep('IB', sample.files)]
  
  # take only pre-vaccine samples
  pre.ind <- c(grep('RL013', sample.files), grep('RL014', sample.files), grep('RL015', sample.files))
  sample.files <- sample.files[pre.ind]
  
  for(j in 1:length(sample.files))
  {
    df <- read.csv(paste0( isotype.dirs[i], '/', sample.files[j]),header= T, stringsAsFactors= F)
    lPlacesCloneId <- gregexpr(pattern ='_',df[1,"CLONE_ID"]) # get places of the char '_' in "CLONE_ID" column
    patientNumber <- substr(df[1,"CLONE_ID"],start = 1,stop = lPlacesCloneId[[1]][1]-1)
    
    
    # check IF/OF
     #is.IF <- c()
      rowNum <- nrow(df)
    #  for(n in 1:rowNum)
    #  {
    #   if(df[n,"INFRAME"]==T && df[n,"STOP"]==F)
    #    {
    #    is.IF[n]<-"IF"
    #  }
    #  else
    #  {
    #    is.IF[n]<-"OF"
    #  }
    #   }
    
    if(i==1)
    {
      df$INFRAME
      
      cloneData_flu.df <- data.frame(sample = rep(patientNumber,rowNum),cloneSize = df$CP_NUM, sequenceStatus=  df$FUNCTIONAL,Iso = rep(isotype.dirs[i],rowNum))
    }else
    {
      cloneData_flu.df <- rbind(cloneData_flu.df,
                            data.frame(sample = rep(patientNumber,rowNum),cloneSize = df$CP_NUM,sequenceStatus= df$FUNCTIONAL,Iso = rep(isotype.dirs[i],rowNum)))
    }
    
  }
}


# save the flu data - 
setwd("/home/nofar/Desktop/Lab")
write.csv(cloneData_flu.df,file = 'fluCloneData.csv')


# Keck data
setwd("/home/nofar/Desktop/Lab")
Keck98.top500.df <- read.csv(file = 'Keck98.top500EditVnames.csv' ,header= T, stringsAsFactors= F)

cloneData_keck.df <- data.frame(sample = Keck98.top500.df$sample,cloneSize = Keck98.top500.df$count,sequenceStatus = Keck98.top500.df$Productive)
cloneData_keck.df$Iso <- "naive_IGHM"

# save the keck data - 
setwd("/home/nofar/Desktop/Lab")
write.csv(cloneData_keck.df,file = 'keckCloneData.csv')

#############################################################################################

setwd("/home/nofar/Desktop/Lab")
pfizer_cloneData.df <- read.csv('pfizerCloneData.csv',header= T, stringsAsFactors= F)
pfizer_cloneData.df <-pfizer_cloneData.df[pfizer_cloneData.df$sequenceStatus==T,]
pfizer_cloneData.df <- arrange(pfizer_cloneData.df, -cloneSize)
pfizer_cloneData.df$dataSet <- "Pfizer"

flu_cloneData.df <- read.csv('fluCloneData.csv',header= T, stringsAsFactors= F)
flu_cloneData.df$dataSet <- "Flu"
flu_cloneData.df <-flu_cloneData.df[flu_cloneData.df$sequenceStatus==T,]
flu_cloneData.df <- arrange(flu_cloneData.df, -cloneSize)

keck_cloneData.df <- read.csv('keckCloneData.csv',header= T, stringsAsFactors= F)
keck_cloneData.df$dataSet <- "Keck"
keck_cloneData.df <-keck_cloneData.df[keck_cloneData.df$sequenceStatus=="productive",]
keck_cloneData.df <- arrange(keck_cloneData.df, -cloneSize)

cloneData.df <- rbind(pfizer_cloneData.df,flu_cloneData.df,keck_cloneData.df)


#arrenge clone fron bigesst to lower
#cloneData.df <- arrange(cloneData.df, -cloneSize)
#pfizerAllCloneSum <- sum(cloneData.df$cloneSize)
#pfizerComulativeSum <-cumsum(cloneData.df$cloneSize)
#pfizerComulativeSumDIVsum <- pfizerComulativeSum/pfizerAllCloneSum
# get the Y of the 1 precent cells
#onePrecPlace <- length(cloneData.df$cloneSize)/100
#Y<- pfizerComulativeSumDIVsum[onePrecPlace]
#cloneData.df$dataSet <- "Pfizer"

b <- ddply(cloneData.df, c("Iso","sample","dataSet"), summarise,
           clonesNumber = length(sample),
           allCloneSize = sum(cloneSize),
           onePrecPlace = clonesNumber/100,
           onePrecPlaceComulativeSum = sum(cloneSize[1:onePrecPlace])/allCloneSize             
)

bA <- ddply(cloneData.df, c("Iso","dataSet"), summarise,
           clonesNumber = length(sample),
           allCloneSize = sum(cloneSize),
           onePrecPlace = clonesNumber/100,
           onePrecPlaceComulativeSum = sum(cloneSize[1:onePrecPlace])/allCloneSize             
)

ggplot(data=bA, aes(x=dataSet, y=onePrecPlaceComulativeSum, fill=Iso))+
  geom_bar(stat="identity",position=position_dodge())+
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle("")+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  xlab("")+
  ylab("///")+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=7))



################################################################# not relevant

b$dataSet <- factor(b$dataSet)
b$Iso <- factor(b$Iso)

beeswarm <- beeswarm(onePrecPlaceComulativeSum ~dataSet, 
                     data = b, method = 'swarm', 
                     corral = 'wrap',
                     pwcol = Iso)[, c(1, 2, 4, 6)]
colnames(beeswarm) <- c("x", "y", "Iso", "DataSet") 


ggplot(data=beeswarm, aes(x=x, y=y, fill=Iso))+
  geom_point(colour="black", size=4, shape=21)+
  # scale_fill_manual(values=cols) +
  theme_bw() +
  ggtitle('')+
  #  coord_cartesian(ylim=c(0,ceiling(max(df.mean$sum)))) +
  xlab("Data")+
  ylab("///")+
  scale_x_continuous(breaks = c(1:3), 
                     labels = c("Flu","Keck","Pfizer" ), expand = c(0, 0.5,1,1.5))+
  theme(  legend.title=element_blank(),
          legend.position=c(1,1), # legend position
          legend.justification=c(1,1),
          panel.grid.minor=element_blank(), # remove grid
          panel.grid.major=element_blank(),
          axis.text.x = element_text(vjust=1),
          axis.title.x=element_blank()
  )+
  guides(fill=guide_legend(nrow=4))




beeswarm(onePrecPlaceComulativeSum ~dataSet , data = b,pch=20,
         pwcol = Iso, xlab = '', 
         ylab = '', 
         labels = c('Pfizer','Flu','Keck'))
legend('right', legend = levels(b$Iso), title = 'Isotype', 
       pch = 16, col = 1:4)


beeswarm.plot <- ggplot(beeswarm, aes(x, y)) +
  xlab("Data") +
  ylab("ka")+
  scale_y_continuous(expression("bbbb"))
beeswarm.plot2 <- beeswarm.plot + geom_boxplot(aes(x, y,
                                                   group = round_any(x, 1, round)), outlier.shape = NA)
beeswarm.plot3 <- beeswarm.plot + geom_point(aes(colour = ER)) +
  scale_colour_manual(values = c("black", "red")) + 
  scale_x_continuous(breaks = c(1:3), 
                     labels = c("Pfizer", "Flu","Keck"), expand = c(0, 0.5,1,1.5))



  
