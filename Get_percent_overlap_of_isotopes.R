library(ggplot2)
library(reshape2)
library(seqinr)
require(gtable)
library(gridExtra)
require(plyr)
require(beeswarm)


# Pfizer

path <- '/home/nofar/Desktop/Lab/CDR3_check/'
setwd(path)

pfizerData.df <- read.csv(file = 'PFIZER_ALL_totalLeavesType_combineAB.tab' ,header= T, stringsAsFactors= F,sep = '\t')
pfizerData.df$Patient <- substr(x = pfizerData.df$treeNUM, start = 13,stop = 15)

Patient.vector = unique(pfizerData.df$Patient)

#for each isotype create table of overlap: 
for(i in 1:length(Patient.vector))
{
  # for each isotype save the top 100 clone
  pfizerData_tmp.df <- pfizerData.df[pfizerData.df$Patient==Patient.vector[i],]
  # naive_IGHM
  naive_IGHM.df <- data.frame(treeNUM = pfizerData_tmp.df$treeNUM, naive_IGHM = pfizerData_tmp.df$naive_IGHM)
  naive_IGHM.df <- arrange(naive_IGHM.df, -naive_IGHM)
  naive_IGHM.df <- naive_IGHM.df[1:100,]
  
  # IGHM
  IGHM.df <- data.frame(treeNUM = pfizerData_tmp.df$treeNUM, IGHM = pfizerData_tmp.df$IGHM)
  IGHM.df <- arrange(IGHM.df, -IGHM)
  IGHM.df <- IGHM.df[1:100,]
  
  # IGHG
  IGHG.df <- data.frame(treeNUM = pfizerData_tmp.df$treeNUM, IGHG = pfizerData_tmp.df$IGHG)
  IGHG.df <- arrange(IGHG.df, -IGHG)
  IGHG.df <- IGHG.df[1:100,]
  
  # IGHA
  IGHA.df <- data.frame(treeNUM = pfizerData_tmp.df$treeNUM, IGHA = pfizerData_tmp.df$IGHA)
  IGHA.df <- arrange(IGHA.df, -IGHA)
  IGHA.df <- IGHA.df[1:100,]
  
  # enter all trees into one table
  top100cloneTreeForEachIsotype.df <- data.frame(naive_IGHM = naive_IGHM.df$treeNUM,
                                                 IGHM = IGHM.df$treeNUM,
                                                 IGHG = IGHG.df$treeNUM,
                                                 IGHA = IGHA.df$treeNUM, stringsAsFactors= F)
  
  # create tabel with the percentage of overlap
  overlapTable <- matrix(data = 0 ,nrow = 4, ncol = 4 ) 
  colnames(overlapTable)  <- colnames(top100cloneTreeForEachIsotype.df)
  rownames(overlapTable)  <- colnames(top100cloneTreeForEachIsotype.df)
  
  # fill matrix
  for(j in 1:nrow(overlapTable))
  {
    for(k in 1:nrow(overlapTable))
    {
      overlapTable[j,k] <- length(which(top100cloneTreeForEachIsotype.df[,j]%in%top100cloneTreeForEachIsotype.df[,k]))/100
    }
  }
  
  write.csv(overlapTable,paste0(Patient.vector[i],"_overlapTop.csv"))
}


########################################################################################################

#~~~~~~ Flu


# Pfizer

path <- '/home/nofar/Desktop/Lab/CDR3_check/'
setwd(path)

fluData.df <- read.csv(file = 'FLU_ALL_totalLeavesType.csv' ,header= T, stringsAsFactors= F)

Patient.vector = unique(fluData.df$patient_name)

#for each isotype create table of overlap: 
for(i in 1:length(Patient.vector))
{
  # for each isotype save the top 100 clone
  fluData_tmp.df <- fluData.df[fluData.df$patient_name==Patient.vector[i],]
  # naive_IGHM
  naive_IGHD.df <- data.frame(treeNUM = fluData_tmp.df$treeNUM, naive_IGHD = fluData_tmp.df$naive_IGHD)
  naive_IGHD.df <- arrange(naive_IGHD.df, -naive_IGHD)
  naive_IGHD.df <- naive_IGHD.df[1:100,]
  
  # IGHM
  IGHM.df <- data.frame(treeNUM = fluData_tmp.df$treeNUM, IGHM = fluData_tmp.df$IGHM)
  IGHM.df <- arrange(IGHM.df, -IGHM)
  IGHM.df <- IGHM.df[1:100,]
  
  IGHG1.df <- data.frame(treeNUM = fluData_tmp.df$treeNUM, IGHG1 = fluData_tmp.df$IGHG.1)
  IGHG1.df <- arrange(IGHG1.df, -IGHG1)
  IGHG1.df <- IGHG1.df[1:100,]
  
  IGHG2.df <- data.frame(treeNUM = fluData_tmp.df$treeNUM, IGHG2 = fluData_tmp.df$IGHG.2)
  IGHG2.df <- arrange(IGHG2.df, -IGHG2)
  IGHG2.df <- IGHG2.df[1:100,]
  
  # IGHG
  
  #IGHG.df <- data.frame(treeNUM = c(IGHG1.df$treeNUM,IGHG2.df$treeNUM), IGHG = c(IGHG1.df$IGHG1,IGHG2.df$IGHG2))
  #IGHG.df <- arrange(IGHG.df, -IGHG)
  #IGHG.df <- IGHG.df[1:100,]
  
  # IGHA
  IGHA.df <- data.frame(treeNUM = fluData_tmp.df$treeNUM, IGHA = fluData_tmp.df$IGHA)
  IGHA.df <- arrange(IGHA.df, -IGHA)
  IGHA.df <- IGHA.df[1:100,]
  
  # IGHE
  IGHE.df <- data.frame(treeNUM = fluData_tmp.df$treeNUM, IGHE = fluData_tmp.df$IGHE)
  IGHE.df <- arrange(IGHE.df, -IGHE)
  IGHE.df <- IGHE.df[1:100,]
  
  # enter all trees into one table
  top100cloneTreeForEachIsotype.df <- data.frame(naive_IGHD = naive_IGHD.df$treeNUM,
                                                 IGHM = IGHM.df$treeNUM,
                                                 IGHG1 = IGHG1.df$treeNUM,
                                                 IGHG2 = IGHG2.df$treeNUM,
                                              #   IGHG = IGHG.df$treeNUM,
                                                 IGHA = IGHA.df$treeNUM,
                                                 IGHE = IGHE.df$treeNUM,stringsAsFactors= F)
  
  # create tabel with the percentage of overlap
  overlapTable <- matrix(data = 0 ,nrow = 6, ncol = 6 ) 
  colnames(overlapTable)  <- colnames(top100cloneTreeForEachIsotype.df)
  rownames(overlapTable)  <- colnames(top100cloneTreeForEachIsotype.df)
  
  # fill matrix
  for(j in 1:nrow(overlapTable))
  {
    for(k in 1:nrow(overlapTable))
    {
      overlapTable[j,k] <- length(which(top100cloneTreeForEachIsotype.df[,j]%in%top100cloneTreeForEachIsotype.df[,k]))/100
    }
  }
  
  write.csv(overlapTable,paste0(Patient.vector[i],"_overlapTop.csv"))
}












