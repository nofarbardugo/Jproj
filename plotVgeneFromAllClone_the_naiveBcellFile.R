
setwd("/home/nofar/Desktop/Lab")

# imports
require(seqinr)
require(ggplot2)
require(plyr)

# read the file into data frame
#naiveBcell.df <- read.table(file = 'naive_B_cells.txt' ,header= T, stringsAsFactors= F,sep = "\t")


# orginaize the V and J Gene for use - remove lines that the vGene is unsure (have 'or')
Vgene.vector <- c()
Jgene.vector <- c()
length <- nrow(naiveBcell.df)
for(i in 1:length)
{
  if(naiveBcell.df[i,"VGeneName"]!="(undefined)")
  {
    Vgene.vector[i] <- naiveBcell.df[i,"VGeneName"]
    Jgene.vector[i] <- naiveBcell.df[i,"JGeneName"]
    l <-nchar(Vgene.vector[i]) # get the string length

    findMakaf <- !is.na(gregexpr(pattern = '-',naiveBcell.df[i,"VGeneName"])[[1]][2])
    
    if(!is.na(gregexpr(pattern = '-',naiveBcell.df[i,"VGeneName"])[[1]][2]))# remove Allele char
    {
      Vgene.vector[i]<-substr(Vgene.vector[i],start = 1,stop = l-2)
    }
    
  }#if(grepl(Keck98.top500.df[i,"V1"], "or")==F)
  else
  {
    Vgene.vector[i]<-NA
  }
  
}

naiveBcell.df$EditV <- Vgene.vector
naiveBcell.df$EditJ <- Jgene.vector

# remove those that Vgene are not certain - NA
naiveBcell.df <- naiveBcell.df[is.na(naiveBcell.df$EditV)==F,]
write.csv(naiveBcell.df ,"Edit_naive_B_cells.csv")

################################################################

naiveBcell.df <- read.table(file = 'naive_B_cells.txt' ,header= T, stringsAsFactors= F,sep = "\t")


####### V
# calc the p(v)
VrepeatAll <- ddply(naiveBcell.df, c("VGeneName","sequenceStatus"), summarise,repeatTimes = length(VGeneName))


vNumIn <- sum(VrepeatAll[VrepeatAll$sequenceStatus=="Productive","repeatTimes"]) 
vNumOut <- sum(VrepeatAll[VrepeatAll$sequenceStatus=="Out of frame","repeatTimes"]) 
vNumStop <- sum(VrepeatAll[VrepeatAll$sequenceStatus=="Has stop","repeatTimes"])
VrepeatAll$frecuency <- 0
for(i in 1:nrow(VrepeatAll))
{
  if(VrepeatAll[i,"sequenceStatus"]=="Productive")
  {
    VrepeatAll[i,"frecuency"] <- VrepeatAll[i,"repeatTimes"]/vNumIn 
  }
  else if(VrepeatAll[i,"sequenceStatus"]=="Out of frame")
  {
    VrepeatAll[i,"frecuency"] <- VrepeatAll[i,"repeatTimes"]/vNumOut 
  }
  else
  {
    VrepeatAll[i,"frecuency"] <- VrepeatAll[i,"repeatTimes"]/vNumStop
  }
  
}

VrepeatAll <- arrange(VrepeatAll, -frecuency)
sort(VrepeatAll, decreasing=TRUE)

VrepeatAll$VGeneName <- factor(as.character(VrepeatAll$VGeneName),levels=as.character(VrepeatAll$VGeneName))
VrepeatAll$VGeneName

vUse.plot <- ggplot(data=VrepeatAll, aes(x=VGeneName, y=frecuency, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)+
  scale_x_discrete() +theme(axis.text.x  = element_text(angle=90))

########  get plot of The percentage of geometric frequency each #######
vFrequency.plot <- ggplot(VrepeatAll, aes(x=VGeneName, y=Avg.frequencyNormalize, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)+ theme(axis.text.x  = element_text(angle=90))



####### J ############
# calc the p(j)
JrepeatAll <- ddply(naiveBcell.df, c("JGeneName","sequenceStatus"), summarise,repeatTimes = length(JGeneName))


jNumIn <- sum(JrepeatAll[JrepeatAll$sequenceStatus=="Productive","repeatTimes"]) 
jNumOut <- sum(JrepeatAll[JrepeatAll$sequenceStatus=="Out of frame","repeatTimes"]) 
jNumStop <- sum(JrepeatAll[JrepeatAll$sequenceStatus=="Has stop","repeatTimes"])
JrepeatAll$frecuency <- 0
for(i in 1:nrow(JrepeatAll))
{
  if(JrepeatAll[i,"sequenceStatus"]=="Productive")
  {
    JrepeatAll[i,"frecuency"] <- JrepeatAll[i,"repeatTimes"]/jNumIn 
  }
  else if(JrepeatAll[i,"sequenceStatus"]=="Out of frame")
  {
    JrepeatAll[i,"frecuency"] <- JrepeatAll[i,"repeatTimes"]/jNumOut 
  }
  else
  {
    JrepeatAll[i,"frecuency"] <- JrepeatAll[i,"repeatTimes"]/jNumStop
  }
  
}

JrepeatAll <- arrange(JrepeatAll, -frecuency)
JrepeatAll$JGeneName <- factor(as.character(JrepeatAll$JGeneName),levels=as.character(JrepeatAll$JGeneName))


jUse.plot <- ggplot(data=JrepeatAll, aes(x=JGeneName, y=frecuency, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)+
  scale_x_discrete() +theme(axis.text.x  = element_text(angle=90))

########  get plot of The percentage of geometric frequency each #######
jFrequency.plot <- ggplot(JrepeatAll, aes(x=JGeneName, y=Avg.frequencyNormalize, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)+ theme(axis.text.x  = element_text(angle=90))



##################################### get plot of The percentage of using each cdr3 length ##############  

##### for the cdr3 length #######
IgHV_CDR3.df <-  read.csv(file = 'IgHV_CDR3.csv' ,header= T, stringsAsFactors= F)

# add column of the v_gene without the allele in order to match it with the clones file
IgHV_CDR3.df$V_GENE_no_Allele <- IgHV_CDR3.df$V_GENE
for(i in 1:273)
{
  IgHV_CDR3.df[i,"V_GENE_no_Allele"] <- paste0(substr(x = IgHV_CDR3.df[i,"V_GENE"],start = 1,stop =  nchar(IgHV_CDR3.df[i,"V_GENE"])-3))
                                               #,substr(x = IgHV_CDR3.df[i,"V_GENE"] , start = 5, stop = nchar(IgHV_CDR3.df[i,"V_GENE"])-3)
                                               #,sep = "0")
}

# find the cdr3 length 

# CDR3_length.df <- data.frame('CDR3_length'= vector(length=nrow(Keck98.top500.df)));
Vgene.vector <- c()
naiveBcell.df$EditCDR3Length 
for(i in 1:nrow(naiveBcell.df))
{  
  IgHV_CDR3Length <- nchar(IgHV_CDR3.df[IgHV_CDR3.df$V_GENE_no_Allele==naiveBcell.df[i,"EditV"],][1,2]) # get the length of the IgHV_CDR3
  
  VEND <- naiveBcell.df[i,"VIndex"] + IgHV_CDR3Length #VEND: vIndex + length(IgHV_CDR3(V))
  
  #find JSTART
  JSTART <- naiveBcell.df[i,"JIndex"] - naiveBcell.df[i,"JDeletion"] +1 # JSTART:jIndex - jDeletion +1
  
  # calculte the cdr3Length
  naiveBcell.df[i,"EditCDR3Length"]  <- JSTART-VEND-1
}





















CDR3repeatAll <- ddply(naiveBcell.df, c("EditCDR3Length ","sequenceStatus"), summarise,
                    repeatTimes = length(EditCDR3Length),
                    Avg.frequencyNormalized = exp(mean(log(frequencyNormalized))))


CDR3NumIn <- sum(CDR3repeatAll[CDR3repeatAll$sequenceStatus=="In","repeatTimes"]) 
CDR3NumOut <- sum(CDR3repeatAll[CDR3repeatAll$sequenceStatus=="Out","repeatTimes"]) 
CDR3NumStop <- sum(CDR3repeatAll[CDR3repeatAll$sequenceStatus=="Stop","repeatTimes"])
CDR3repeatAll$frecuency <- 0
for(i in 1:nrow(CDR3repeatAll))
{
  if(CDR3repeatAll[i,"sequenceStatus"]=="In")
  {
    CDR3repeatAll[i,"frecuency"] <- CDR3repeatAll[i,"repeatTimes"]/CDR3NumIn 
  }
  else if(CDR3repeatAll[i,"sequenceStatus"]=="Out")
  {
    CDR3repeatAll[i,"frecuency"] <- CDR3repeatAll[i,"repeatTimes"]/CDR3NumOut 
  }
  else
  {
    CDR3repeatAll[i,"frecuency"] <- CDR3repeatAll[i,"repeatTimes"]/CDR3NumStop
  }
  
}


cdr3Use.plot <- ggplot(data=CDR3repeatAll, aes(x=EditCDR3Length, y=frecuency, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)

########  get plot of The percentage of geometric frequency each #######
jFrequency.plot <- ggplot(CDR3repeatAll, aes(x=EditCDR3Length, y=Avg.frequencyNormalized, fill=sequenceStatus)) + geom_bar(colour="black", stat="identity")+
  facet_wrap(~ sequenceStatus, ncol = 1)






