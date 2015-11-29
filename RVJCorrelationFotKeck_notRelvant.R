
setwd("/home/nofar/Desktop/Lab")

# imports
require(seqinr)
require(ggplot2)
require(plyr)


# read the file into data frame
Keck98.top500.df <- read.table(file = 'Keck98.top500.tsv' ,header= T, stringsAsFactors= F,sep = "\t")
#Keck98.top500IN.df <- Keck98.top500.df[Keck98.top500.df$sequenceStatus=="In",] 
#Keck98.top500Out.df <- Keck98.top500.df[Keck98.top500.df$sequenceStatus=="Out",] 
Keck98.top500.df$frequencyNormalized <- Keck98.top500.df$frequencyNormalized/100 

# orginaize the V and J Gene for use - remove lines that the vGene is unsure (have 'or'  )
Vgene.vector <- c()
Jgene.vector <- c()
length <- nrow(Keck98.top500.df)
for(i in 1:length)
{
  
  if(grepl( "or",Keck98.top500.df[i,"V"])==F)
  {
    Vgene.vector[i]<- strsplit(Keck98.top500.df[i,"V"],",")[[1]][1] # get Vgene form clone file - first one
    Jgene.vector[i]<- strsplit(Keck98.top500.df[i,"J"],",")[[1]][1] # get Vgene form clone file - first one
    
    l <-nchar(Vgene.vector[i]) # get the string length
    charToCheck <-substr(Vgene.vector[i],start = l-1,stop = l-1)
     
    if(charToCheck == "_")# remove Allele char
    {
      Vgene.vector[i]<-substr(Vgene.vector[i],start = 1,stop = l-2)
    }
    
  }#if(grepl(Keck98.top500.df[i,"V1"], "or")==F)
  else
  {
    Vgene.vector[i]<-NA
  }

}

Vgene.vector <- Vgene.vector[!is.na(Vgene.vector)]
Vgene <- unique(Vgene.vector)
d <- data.frame("Vgene" = Vgene)
write.csv(d,"Keck_unique_Vgene.csv")

Keck98.top500.df$EditV <- gsub('V0','V',Vgene.vector)
Keck98.top500.df$EditJ <- gsub('J0','J',Jgene.vector)

# remove those that Vgene are not certain - NA
Keck98.top500.df <- Keck98.top500.df[is.na(Keck98.top500.df$EditV)==F,] 


##### for the cdr3 length #######
IgHV_CDR3.df <-  read.csv(file = 'IgHV_CDR3.csv' ,header= T, stringsAsFactors= F)

# add column of the v_gene without the allele in order to match it with the clones file
IgHV_CDR3.df$V_GENE_no_Allele <- IgHV_CDR3.df$V_GENE
for(i in 1:273)
{
  IgHV_CDR3.df[i,"V_GENE_no_Allele"] <- paste0(substr(x = IgHV_CDR3.df[i,"V_GENE"],start = 1,stop =  nchar(IgHV_CDR3.df[i,"V_GENE"])-3)
                                              #,substr(x = IgHV_CDR3.df[i,"V_GENE"] , start = 5, stop = nchar(IgHV_CDR3.df[i,"V_GENE"])-3)
                                              #,sep = "0")
}


# find the cdr3 length 
# CDR3_length.df <- data.frame('CDR3_length'= vector(length=nrow(Keck98.top500.df)));


Vgene.vector <- c()
Keck98.top500.df$EditCDR3Length 
for(i in 1:nrow(Keck98.top500.df))
{
  # find VEnd - get the first one
  if(is.na(Keck98.top500.df[i,"EditV"]))
  {
    Vgene<- strsplit(Keck98.top500.df[i,"V"],",")[[1]][1] # get Vgene form clone file - first one
    Vgene.vector[i]<- Vgene
  }
  else
  {
    Vgene<- Keck98.top500.df[i,"EditV"] # get Vgene form clone file - first one
    Vgene.vector[i]<- Vgene
  }

  IgHV_CDR3Length <- nchar(IgHV_CDR3.df[IgHV_CDR3.df$V_GENE_no_Allele==Vgene,][1,2]) # get the length of the IgHV_CDR3
  
  VEND <- Keck98.top500.df[i,"vIndex"] + IgHV_CDR3Length #VEND: vIndex + length(IgHV_CDR3(V))
  
  #find JSTART
  JSTART <- Keck98.top500.df[i,"jIndex"] - Keck98.top500.df[i,"jDeletion"] +1 # JSTART:jIndex - jDeletion +1
  
  # calculte the cdr3Length
  Keck98.top500.df[i,"EditCDR3Length"]  <- JSTART-VEND-1
}


#plot distrbution of the cdr3 length
#ggplot(  CDR3_length.df, aes(x=CDR3_length)) + geom_histogram(binwidth=.5, colour="black", fill="white")

#CDR3repeat <- ddply(CDR3_length.df, c("CDR3_length"), summarise,repeatTimes = length(CDR3_length))
#########################################################################################################

###### calc RVJ ######

# calc the p(v)
Vrepeat <- ddply(Keck98.top500.df, c("EditV","sequenceStatus"), summarise,repeatTimes = length(EditV),Avg.frequencyNormalized = exp(mean(log(frequencyNormalized))))
vNum <- sum(Vrepeat[,"repeatTimes"]) 
Vrepeat$frecuency <- Vrepeat[,"repeatTimes"]/vNum
#calc the p(J)
Jrepeat <- ddply(Keck98.top500.df, c("EditJ"), summarise,repeatTimes = length(EditJ),Avg.frequencyNormalized =exp(mean(log(frequencyNormalized))))
#Jrepeat<- Jrepeat[1:nrow(Jrepeat)-1,]
jNum <- sum(Jrepeat[,"repeatTimes"]) 
Jrepeat$frecuency <- Jrepeat[,"repeatTimes"]/jNum
#calc the p(VJ)
VJrepeat <- ddply(Keck98.top500.df, c("EditV","EditJ"), summarise,repeatTimes = length(EditJ),Avg.frequencyNormalized =exp(mean(log(frequencyNormalized))))
#VJrepeat<- VJrepeat[1:nrow(VJrepeat)-1,]
VJrepeat$dup <-paste0(VJrepeat[,"EditV"],VJrepeat[,"EditJ"])
vjNum <- sum(VJrepeat[,"repeatTimes"]) 
VJrepeat$frecuency <- VJrepeat[,"repeatTimes"]/vjNum



#calc RVJ for each line 
Keck98.top500.df$RVJ <- 0
for(i in 1:nrow(Keck98.top500.df))
{
  print(i)
  if(!(is.na(Keck98.top500.df[i,"EditV"])))
  {
    pv <- Vrepeat[Vrepeat[,"EditV"]==Keck98.top500.df[i,"EditV"],"frecuency"][1]
    pj <- Jrepeat[Jrepeat[,"EditJ"]==Keck98.top500.df[i,"EditJ"],"frecuency"][1]
    pVJ <- VJrepeat[VJrepeat[,"dup"]==paste0(Keck98.top500.df[i,"EditV"],Keck98.top500.df[i,"EditJ"]),"frecuency"][1]
    
    if(pv !=0 && pj!=0)
    {
      Keck98.top500.df[i,"RVJ"] <- pVJ/(pv*pj)
    }
  }

}

# calc coralation between RVJ to "frequencyNormalized"
require(Hmisc)
cor=rcorr(x=Keck98.top500.df$RVJ,y=Keck98.top500.df$frequencyNormalized, type="spearman")
R=cor$r
pVal=cor$P

corCDR3=rcorr(x=Keck98.top500.df$cdr3Length,y=Keck98.top500.df$frequencyNormalized, type="spearman")

x <-ggplot( Keck98.top500.df, aes(x=RVJ,y= frequencyNormalized, group=1,color = sample)) + geom_point()+
facet_wrap(~ sample)

#####################################  V ##############  

#df <- arrange(Keck98.top500.df, sample,-frequencyNormalized)


IF_df_all500 <- Keck98.top500.df[Keck98.top500.df$sequenceStatus=="In",]

IF_Vrepeat_all500 <- ddply(IF_df_all500, c("EditV","sequenceStatus"), summarise,
                 repeatTimes = length(EditV),
                 Avg.frequencyNormalized = exp(mean(log(frequencyNormalized))))

IF_Vrepeat_all500$frecuency <- IF_Vrepeat_all500$repeatTimes/sum(IF_Vrepeat_all500$repeatTimes)
# 
IF_Vrepeat_all500 <- arrange(IF_Vrepeat_all500, -frecuency)
IF_Vrepeat_all500$EditV <- factor(as.character(IF_Vrepeat_all500$EditV),levels=as.character(IF_Vrepeat_all500$EditV))
# plot
vUse_IF_Vrepeat_all500.plot <- ggplot(IF_Vrepeat_all500, aes(x=EditV, y=frecuency, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)+ theme(axis.text.x  = element_text(angle=90))


# save only the first top 100 in clone for each sample 
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

IF_Vrepeat_first500 <- ddply(dfFirst100, c("EditV","sequenceStatus"), summarise,
                 repeatTimes = length(EditV),
                 Avg.frequencyNormalized = exp(mean(log(frequencyNormalized))))

IF_Vrepeat_first500$frecuency <- IF_Vrepeat_first500$repeatTimes/sum(IF_Vrepeat_first500$repeatTimes)


IF_Vrepeat_first500$EditV <- factor(as.character(IF_Vrepeat_first500$EditV),levels=as.character(IF_Vrepeat_all500$EditV))

vUse_IF_Vrepeat_first500.plot<- ggplot(IF_Vrepeat_first500, aes(x=EditV, y=frecuency, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)+ theme(axis.text.x  = element_text(angle=90))


-----------------------------------------------------------------------------
########  get plot of The percentage of using each #######
vNumIn <- sum(Vrepeat[Vrepeat$sequenceStatus=="In","repeatTimes"]) 
vNumOut <- sum(Vrepeat[Vrepeat$sequenceStatus=="Out","repeatTimes"]) 
vNumStop <- sum(Vrepeat[Vrepeat$sequenceStatus=="Stop","repeatTimes"])
Vrepeat$frecuency <- 0
for(i in 1:nrow(Vrepeat))
{
  if(Vrepeat[i,"sequenceStatus"]=="In")
  {
    Vrepeat[i,"frecuency"] <- Vrepeat[i,"repeatTimes"]/vNumIn 
  }
  else if(Vrepeat[i,"sequenceStatus"]=="Out")
  {
    Vrepeat[i,"frecuency"] <- Vrepeat[i,"repeatTimes"]/vNumOut 
  }
  else
  {
    Vrepeat[i,"frecuency"] <- Vrepeat[i,"repeatTimes"]/vNumStop
  }
 
}

--------------------------------------------------------------
v <- gsub("IGHV0?","",Vrepeat$EditV)
v <- head(v,length(levels))
levels <- gsub("IGHV0?","",as.character(VrepeatAll$VGeneName))
arrange(Vrepeat, -frecuency)


Vrepeat$EditV <-
  Vrepeat$xx <-factor(as.character(Vrepeat$EditV),levels=as.character(VrepeatAll$VGeneName))
------------------------------------------------------------------------------------------------
  
IF_Vrepeat_all500$EditV <- factor(as.character(IF_Vrepeat_all500$EditV),levels=as.character(IF_Vrepeat_all500$EditV))

vUse.plot <- ggplot(Vrepeat, aes(x=EditV, y=frecuency, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
             facet_wrap(~ sequenceStatus, ncol = 1)+ theme(axis.text.x  = element_text(angle=90))+
  


########  get plot of The percentage of geometric frequency each #######
vFrequency.plot <- ggplot(Vrepeat, aes(x=EditV, y=Avg.frequencyNormalized, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1) + theme(axis.text.x  = element_text(angle=90))

---------------------------------------------------------------------------------------------------------------
##################################### get plot of The percentage of using each J ##############  

Jrepeat <- ddply(Keck98.top500.df, c("EditJ","sequenceStatus"), summarise,
                 repeatTimes = length(EditJ),
                 Avg.frequencyNormalized = exp(mean(log(frequencyNormalized))))

jNumIn <- sum(Jrepeat[Jrepeat$sequenceStatus=="In","repeatTimes"]) 
jNumOut <- sum(Jrepeat[Jrepeat$sequenceStatus=="Out","repeatTimes"]) 
jNumStop <- sum(Jrepeat[Jrepeat$sequenceStatus=="Stop","repeatTimes"])
Jrepeat$frecuency <- 0
for(i in 1:nrow(Jrepeat))
{
  if(Jrepeat[i,"sequenceStatus"]=="In")
  {
    Jrepeat[i,"frecuency"] <- Jrepeat[i,"repeatTimes"]/jNumIn 
  }
  else if(Jrepeat[i,"sequenceStatus"]=="Out")
  {
    Jrepeat[i,"frecuency"] <- Jrepeat[i,"repeatTimes"]/jNumOut 
  }
  else
  {
    Jrepeat[i,"frecuency"] <- Jrepeat[i,"repeatTimes"]/jNumStop
  }
  
}

#----------------------------------
v <- gsub("IGHV0?","",Jrepeat$EditJ)
v <- head(v,length(levels))
levels <- gsub("IGHV0?","",as.character(JrepeatAll$JGeneName))
arrange(Vrepeat, -frecuency)
Vrepeat$EditV <-
  Vrepeat$xx <-factor(as.character(v),levels=as.character(levels))
#------------------------------------------------------------------


as.character(JrepeatAll$JGeneName)

jUse.plot <- ggplot(data=Jrepeat, aes(x=EditJ, y=frecuency, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)

########  get plot of The percentage of geometric frequency each #######
jFrequency.plot <- ggplot(Jrepeat, aes(x=EditJ, y=Avg.frequencyNormalized, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)


##################################### get plot of The percentage of using each cdr3 length ##############  

CDR3repeat <- ddply(Keck98.top500.df, c("EditCDR3Length ","sequenceStatus"), summarise,
                 repeatTimes = length(EditCDR3Length),
                 Avg.frequencyNormalized = exp(mean(log(frequencyNormalized))))


CDR3NumIn <- sum(CDR3repeat[CDR3repeat$sequenceStatus=="In","repeatTimes"]) 
CDR3NumOut <- sum(CDR3repeat[CDR3repeat$sequenceStatus=="Out","repeatTimes"]) 
CDR3NumStop <- sum(CDR3repeat[CDR3repeat$sequenceStatus=="Stop","repeatTimes"])
CDR3repeat$frecuency <- 0
for(i in 1:nrow(CDR3repeat))
{
  if(CDR3repeat[i,"sequenceStatus"]=="In")
  {
    CDR3repeat[i,"frecuency"] <- CDR3repeat[i,"repeatTimes"]/CDR3NumIn 
  }
  else if(CDR3repeat[i,"sequenceStatus"]=="Out")
  {
    CDR3repeat[i,"frecuency"] <- CDR3repeat[i,"repeatTimes"]/CDR3NumOut 
  }
  else
  {
    CDR3repeat[i,"frecuency"] <- CDR3repeat[i,"repeatTimes"]/CDR3NumStop
  }
  
}


cdr3Use.plot <- ggplot(data=CDR3repeat, aes(x=EditCDR3Length, y=frecuency, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)

########  get plot of The percentage of geometric frequency each #######
cdr3Frequency.plot <- ggplot(CDR3repeat, aes(x=EditCDR3Length, y=Avg.frequencyNormalized, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)







########  get plot of The percentage of geometric frequency each #######




##############################################################################
# גרף לציור קורלציה 
p <- ggplot(Keck98.top500.df, aes(RVJ, frequencyNormalized)) + 
  geom_smooth(method = "lm",color="black")  + 
  geom_point(aes(color=Keck98.top500.df$group),size=2) + 
  xlab("ADAR1") + ylab("CD8 betta") +
  theme(legend.title=element_blank(), legend.text = element_text(size = 6), panel.grid.major = element_line(size = 0.1)) + 
  scale_color_manual("",values=c("red", "red","red", "red", "red", "red", "blue"))

cor=rcorr(x=Keck98.top500.df$RVJ,y=Keck98.top500.df$frequencyNormalized, type="spearman")
cor=rcorr(x=expression_ADARs_matrix$ADAR,y=expression_ADARs_matrix$FOXP3, type="pearson")

R=cor$r[1]
pVal=cor$P[1]

p + annotate("text", x = 2000, y = 0.45, label = paste("Pearson R =" , round(R,digits = 3)),size=5) + 
  annotate("text", x = 2000, y = 0.4, label = paste("pVal = " , format(pVal,scientific = TRUE)), size=5)

x <- ggplot( Keck98.top500.df, aes(x=RVJ,y= frequencyNormalized, group=1)) + geom_point()
x+ coord_trans(y="log2")

####################################### inframe details ###

Infreame<- ddply(Keck98.top500.df, c("sample","sequenceStatus"), summarise,
                 inframeSum = length(sequenceStatus))

Infreame<-Infreame[Infreame$sequenceStatus=="In",]
Infreame$prec <- Infreame$inframeSum/500
meanIF <-mean(Infreame$prec)
GmeanIf <-mean(exp(mean(log(Infreame$prec))))
SD <-sd(Infreame$prec)
SE <-SD/(sqrt(98))

#########################################################