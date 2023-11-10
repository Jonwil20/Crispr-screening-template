## By comparing the abundance of CrispR/shRNA in test conditions versus control conditions we can 
## assess which genes have a functional role related to the test condition. 
## i.e. If a set sgRNA guides targeting are enriched after drug treatment, that gene is important to 
## susceptiabilty to the drug.

install.packages("BiocManager", repos = "https://cran.r-project.org")
BiocManager::install()
BiocManager::install("ShortRead")
BiocManager::install("GenomeInfoDbData")
BiocManager::install("GenomeInfoDb", force = TRUE)
library(GenomeInfoDb)
library(GenomeInfoDbData)
library(ShortRead)

setwd("C:/Users/jrichards/Downloads/PinAPL-py_demo_data")
#Example data download, Download of all fasta files further down in the code (lines 56-73)
fqSample <- FastqSampler("ToxA_R1_98_S2_L008_R1_001_x.fastq.gz",n=10^6)
fastq <- yield(fqSample)


params <-
list(isSlides = "no")

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
AsSlides <- TRUE
library(DESeq2)
library(GenomicAlignments)
library(Rsubread)


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Crispr Screenining analysis in R

---
"    
  )
  
}



## ----downloadFile,message=FALSE-----------------------------------------------
download.file("https://github.com/LewisLabUCSD/PinAPL-Py/archive/master.zip","pieappleData.zip")
unzip("pieappleData.zip")
download.file("http://pinapl-py.ucsd.edu/example-data","TestData.zip")
unzip("TestData.zip")
download.file("http://pinapl-py.ucsd.edu/run/download/example-run","Results.zip")
unzip("Results.zip")


## ----shortreada,include=FALSE-------------------------------------------------
library(ShortRead)

## ----crsiprRep1Reads,echo=T,eval=T--------------------------------------------
fqSample1 <- FastqSampler("ToxA_R1_98_S2_L008_R1_001_x.fastq.gz",n=10^6)
fastq1 <- yield(fqSample1)

fqSample2 <- FastqSampler("ToxA_R2_S2_L005_R1_001_x.fastq.gz",n=10^6)
fastq2 <- yield(fqSample2)

fqSample3 <- FastqSampler("ToxB_R1_98_S4_L008_R1_001_x.fastq.gz",n=10^6)
fastq3 <- yield(fqSample3)

fqSample4 <- FastqSampler("ToxB_R2_S4_L005_R1_001_x.fastq.gz",n=10^6)
fastq4 <- yield(fqSample4)


fqSample5 <- FastqSampler("Control_R2_S15_L008_R1_001_x.fastq.gz",n=10^6)
fastq5 <- yield(fqSample5)

fqSample6 <- FastqSampler("Control_R1_S14_L008_R1_001_x.fastq.gz",n=10^6)
fastq6 <- yield(fqSample6)


## ----mycRep1ReadsShortReadQ---------------------------------------------------
fastq1
fastq2
fastq3
fastq4
fastq5
fastq6


## ----mycRep1ReadsAccessor-----------------------------------------------------
readSequences1 <- sread(fastq1)
readQuality1 <- quality(fastq1)
#readIDs1 <- vctrs::vec_group_id(fastq1)
readSequences1

readSequences2 <- sread(fastq2)
readQuality2 <- quality(fastq2)
#readIDs1 <- vctrs::vec_group_id(fastq1)
readSequences2

readSequences3 <- sread(fastq3)
readQuality3 <- quality(fastq3)
#readIDs1 <- vctrs::vec_group_id(fastq1)
readSequences3

readSequences4 <- sread(fastq4)
readQuality4 <- quality(fastq4)
#readIDs1 <- vctrs::vec_group_id(fastq1)
readSequences4

readSequences5 <- sread(fastq5)
readQuality5 <- quality(fastq5)
#readIDs1 <- vctrs::vec_group_id(fastq1)
readSequences1

readSequences6 <- sread(fastq6)
readQuality6 <- quality(fastq6)
#readIDs1 <- vctrs::vec_group_id(fastq1)
readSequences6

## ----mycRep1ReadsQScores------------------------------------------------------

readQualities1 <- alphabetScore(readQuality1)
readQualities1[1:10]

readQualities2 <- alphabetScore(readQuality2)
readQualities2[1:10]

readQualities3 <- alphabetScore(readQuality3)
readQualities3[1:10]

readQualities4 <- alphabetScore(readQuality4)
readQualities4[1:10]

readQualities5 <- alphabetScore(readQuality5)
readQualities5[1:10]

readQualities6 <- alphabetScore(readQuality6)
readQualities6[1:10]

## ----mycRep1ReadsQScoresPlot,fig.height=3,fig.width=8-------------------------
library(ggplot2)
toPlot1 <- data.frame(ReadQ=readQualities1)
ggplot(toPlot1,aes(x=ReadQ))+geom_histogram(colour = "black", fill = "cyan")+theme_minimal()+ ggtitle("ToxA_R1")

toPlot2 <- data.frame(ReadQ=readQualities2)
ggplot(toPlot2,aes(x=ReadQ))+geom_histogram(colour = "black", fill = "forestgreen")+theme_minimal()+ ggtitle("ToxA_R2")

toPlot3 <- data.frame(ReadQ=readQualities3)
ggplot(toPlot3,aes(x=ReadQ))+geom_histogram(colour = "black", fill = "magenta")+theme_minimal()+ ggtitle("ToxB_R1")

toPlot4 <- data.frame(ReadQ=readQualities4)
ggplot(toPlot4,aes(x=ReadQ))+geom_histogram(colour = "black", fill = "lightgray")+theme_minimal()+ ggtitle("ToxB_R2")

toPlot5<- data.frame(ReadQ=readQualities5)
ggplot(toPlot5,aes(x=ReadQ))+geom_histogram(colour = "black", fill = "midnightblue")+theme_minimal()+ ggtitle("Control_R2")

toPlot6 <- data.frame(ReadQ=readQualities6)
ggplot(toPlot6,aes(x=ReadQ))+geom_histogram(colour = "black", fill = "red")+theme_minimal() + ggtitle("Control_R1")



## ----mycRep1ReadsAlpFreq------------------------------------------------------
readSequences1 <- sread(fastq1)
readSequences_AlpFreq1 <- alphabetFrequency(readSequences1)
readSequences_AlpFreq1[1:3,]

readSequences2 <- sread(fastq2)
readSequences_AlpFreq2 <- alphabetFrequency(readSequences2)
readSequences_AlpFreq2[1:3,]

readSequences3 <- sread(fastq3)
readSequences_AlpFreq3 <- alphabetFrequency(readSequences3)
readSequences_AlpFreq3[1:3,]

readSequences4 <- sread(fastq4)
readSequences_AlpFreq4 <- alphabetFrequency(readSequences4)
readSequences_AlpFreq4[1:3,]

readSequences5 <- sread(fastq5)
readSequences_AlpFreq5 <- alphabetFrequency(readSequences5)
readSequences_AlpFreq5[1:3,]

readSequences6 <- sread(fastq6)
readSequences_AlpFreq6 <- alphabetFrequency(readSequences6)
readSequences_AlpFreq6[1:3,]


## ----mycRep1ReadsAlpFreqSum---------------------------------------------------
summed__AlpFreq1  <- colSums(readSequences_AlpFreq1)
summed__AlpFreq1[c("A","C","G","T","N")]

summed__AlpFreq2  <- colSums(readSequences_AlpFreq2)
summed__AlpFreq2[c("A","C","G","T","N")]

summed__AlpFreq3  <- colSums(readSequences_AlpFreq3)
summed__AlpFreq3[c("A","C","G","T","N")]

summed__AlpFreq4  <- colSums(readSequences_AlpFreq4)
summed__AlpFreq4[c("A","C","G","T","N")]

summed__AlpFreq5  <- colSums(readSequences_AlpFreq5)
summed__AlpFreq5[c("A","C","G","T","N")]

summed__AlpFreq6  <- colSums(readSequences_AlpFreq6)
summed__AlpFreq6[c("A","C","G","T","N")]

## ----mycRep1ReadsAlpByCycle---------------------------------------------------
readSequences_AlpbyCycle1 <- alphabetByCycle(readSequences1)
readSequences_AlpbyCycle1[1:4,1:10]

readSequences_AlpbyCycle2 <- alphabetByCycle(readSequences2)

readSequences_AlpbyCycle3 <- alphabetByCycle(readSequences3)

readSequences_AlpbyCycle4 <- alphabetByCycle(readSequences4)

readSequences_AlpbyCycle5 <- alphabetByCycle(readSequences5)

readSequences_AlpbyCycle6 <- alphabetByCycle(readSequences6)



## ----mycRep1ReadsAlpByCyclePlot-----------------------------------------------
AFreq1 <- readSequences_AlpbyCycle1["A",]
CFreq1 <- readSequences_AlpbyCycle1["C",]
GFreq1 <- readSequences_AlpbyCycle1["G",]
TFreq1 <- readSequences_AlpbyCycle1["T",]
toPlot1b <- data.frame(Count=c(AFreq1,CFreq1,GFreq1,TFreq1),
                     Cycle=rep(1:max(width(readSequences1)),4),
                     Base=rep(c("A","C","G","T"),each=max(width(readSequences1))))

AFreq2 <- readSequences_AlpbyCycle2["A",]
CFreq2 <- readSequences_AlpbyCycle2["C",]
GFreq2 <- readSequences_AlpbyCycle2["G",]
TFreq2 <- readSequences_AlpbyCycle2["T",]
toPlot2b <- data.frame(Count=c(AFreq2,CFreq2,GFreq2,TFreq2),
                       Cycle=rep(1:max(width(readSequences2)),4),
                       Base=rep(c("A","C","G","T"),each=max(width(readSequences2))))

AFreq3 <- readSequences_AlpbyCycle3["A",]
CFreq3 <- readSequences_AlpbyCycle3["C",]
GFreq3 <- readSequences_AlpbyCycle1["G",]
TFreq3 <- readSequences_AlpbyCycle3["T",]
toPlot3b <- data.frame(Count=c(AFreq3,CFreq3,GFreq3,TFreq3),
                       Cycle=rep(1:max(width(readSequences3)),4),
                       Base=rep(c("A","C","G","T"),each=max(width(readSequences3))))

AFreq4 <- readSequences_AlpbyCycle4["A",]
CFreq4 <- readSequences_AlpbyCycle4["C",]
GFreq4 <- readSequences_AlpbyCycle4["G",]
TFreq4 <- readSequences_AlpbyCycle4["T",]
toPlot4b <- data.frame(Count=c(AFreq4,CFreq4,GFreq4,TFreq4),
                       Cycle=rep(1:max(width(readSequences4)),4),
                       Base=rep(c("A","C","G","T"),each=max(width(readSequences4))))

AFreq5 <- readSequences_AlpbyCycle5["A",]
CFreq5 <- readSequences_AlpbyCycle5["C",]
GFreq5 <- readSequences_AlpbyCycle5["G",]
TFreq5 <- readSequences_AlpbyCycle5["T",]
toPlot5b <- data.frame(Count=c(AFreq5,CFreq5,GFreq5,TFreq5),
                       Cycle=rep(1:max(width(readSequences5)),4),
                       Base=rep(c("A","C","G","T"),each=max(width(readSequences5))))

AFreq6 <- readSequences_AlpbyCycle6["A",]
CFreq6 <- readSequences_AlpbyCycle6["C",]
GFreq6 <- readSequences_AlpbyCycle6["G",]
TFreq6 <- readSequences_AlpbyCycle6["T",]
toPlot6b <- data.frame(Count=c(AFreq6,CFreq6,GFreq6,TFreq6),
                       Cycle=rep(1:max(width(readSequences6)),4),
                       Base=rep(c("A","C","G","T"),each=max(width(readSequences6))))

## ----mycRep1ReadsAlpByCyclePlot2,fig.height=4,fig.width=8---------------------

ggplot(toPlot1b,aes(y=Count,x=Cycle,colour=Base))+geom_line()+
  theme_bw()+ggtitle("ToxA_R1")

ggplot(toPlot2b,aes(y=Count,x=Cycle,colour=Base))+geom_line()+
  theme_bw()+ggtitle("ToxA_R2")

ggplot(toPlot3b,aes(y=Count,x=Cycle,colour=Base))+geom_line()+
  theme_bw()+ggtitle("ToxB_R1")

ggplot(toPlot4b,aes(y=Count,x=Cycle,colour=Base))+geom_line()+
  theme_bw()+ggtitle("ToxB_R2")

ggplot(toPlot5b,aes(y=Count,x=Cycle,colour=Base))+geom_line()+
  theme_bw()+ggtitle("Control_R2")

ggplot(toPlot6b,aes(y=Count,x=Cycle,colour=Base))+geom_line()+
  theme_bw()+ggtitle("Control_R1")

## ----mycRep1ReadsQByCycle-----------------------------------------------------
qualAsMatrix1 <- as(readQuality1,"matrix")
qualAsMatrix1[1:2,]

qualAsMatrix2 <- as(readQuality2,"matrix")

qualAsMatrix3 <- as(readQuality3,"matrix")

qualAsMatrix4 <- as(readQuality4,"matrix")

qualAsMatrix5 <- as(readQuality5,"matrix")

qualAsMatrix6 <- as(readQuality6,"matrix")


## ----mycRep1ReadsQByCyclePlotfig.width=8,fig.height=4-------------------------
toPlot1c <- colMeans(qualAsMatrix1)
plot(toPlot1c)

toPlot2c <- colMeans(qualAsMatrix2)
plot(toPlot2c)

toPlot3c <- colMeans(qualAsMatrix3)
plot(toPlot3c)

toPlot4c <- colMeans(qualAsMatrix4)
plot(toPlot4c)

toPlot5c <- colMeans(qualAsMatrix5)
plot(toPlot5c)

toPlot6c <- colMeans(qualAsMatrix6)
plot(toPlot6c)



## ----fa2----------------------------------------------------------------------
GeCKO <- read.delim("GeCKOv21_Human.tsv")
GeCKO[1:2,]


## ----fa3----------------------------------------------------------------------
require(Biostrings)
sgRNAs <- DNAStringSet(GeCKO$seq)
names(sgRNAs) <- GeCKO$UID


## ----fa4----------------------------------------------------------------------
writeXStringSet(sgRNAs,
                file="GeCKO.fa")


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
library(Rsubread)
buildindex("GeCKO","GeCKO.fa", 
           indexSplit=FALSE)



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs1 <- "ToxA_R1_98_S2_L008_R1_001_x.fastq.gz"
myMapped1 <- align("GeCKO",myFQs1,output_file = gsub(".fastq.gz",".bam",myFQs1),
                  nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA")



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped1 



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs2 <- "ToxA_R2_S2_L005_R1_001_x.fastq.gz"
myMapped2 <- align("GeCKO",myFQs2,output_file = gsub(".fastq.gz",".bam",myFQs2),
                  nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1)



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped2 



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs3 <- "ToxB_R1_98_S4_L008_R1_001_x.fastq.gz"
myMapped3 <- align("GeCKO",myFQs3,output_file = gsub(".fastq.gz",".bam",myFQs3),
                  nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1,
                  maxMismatches = 0,indels = 0)
## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped3 



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs4 <- "ToxB_R2_S4_L005_R1_001_x.fastq.gz"
myMapped4 <- align("GeCKO",myFQs4,output_file = gsub(".fastq.gz",".bam",myFQs4),
                   nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA")



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped4 



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs5 <- "Control_R1_S14_L008_R1_001_x.fastq.gz"
myMapped5 <- align("GeCKO",myFQs5,output_file = gsub(".fastq.gz",".bam",myFQs5),
                   nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1)



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped5 



## ---- echo=TRUE, eval=TRUE----------------------------------------------------
myFQs6 <- "Control_R2_S15_L008_R1_001_x.fastq.gz"
myMapped6 <- align("GeCKO",myFQs6,output_file = gsub(".fastq.gz",".bam",myFQs6),
                   nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1,
                   maxMismatches = 0,indels = 0)



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------

myMapped6 



## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
require(GenomicAlignments)
temp1 <- readGAlignments(gsub(".fastq.gz",".bam",myFQs1))
temp1

temp2 <- readGAlignments(gsub(".fastq.gz",".bam",myFQs2))
temp2

temp3 <- readGAlignments(gsub(".fastq.gz",".bam",myFQs3))
temp3

temp4 <- readGAlignments(gsub(".fastq.gz",".bam",myFQs4))
temp4

temp5 <- readGAlignments(gsub(".fastq.gz",".bam",myFQs5))
temp5

temp6 <- readGAlignments(gsub(".fastq.gz",".bam",myFQs6))
temp6

## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
temp1


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
cigar1 <- cigar(temp1)
cigar1[1:5]

cigar2 <- cigar(temp2)
cigar2[1:5]

cigar3 <- cigar(temp3)

cigar4 <- cigar(temp4)

cigar5 <- cigar(temp5)

cigar6 <- cigar(temp6)

## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
cigarRLE2 <- cigarToRleList(cigar2)
cigarRLE2[1]

cigarRLE3 <- cigarToRleList(cigar3)
cigarRLE3[1]

cigarRLE4 <- cigarToRleList(cigar4)
##cigarRLE4[1]

cigarRLE5 <- cigarToRleList(cigar5)
cigarRLE5[1]

cigarRLE6 <- cigarToRleList(cigar6)
cigarRLE6[1]


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
cigarMat2 <- matrix(as.vector(unlist(cigarRLE2)),ncol=50,byrow = TRUE)
cigarMat2[1:2,]

cigarMat3 <- matrix(as.vector(unlist(cigarRLE3)),ncol=50,byrow = TRUE)
cigarMat4 <- matrix(as.vector(unlist(cigarRLE4)),ncol=50,byrow = TRUE)
cigarMat5 <- matrix(as.vector(unlist(cigarRLE5)),ncol=50,byrow = TRUE)
cigarMat6 <- matrix(as.vector(unlist(cigarRLE6)),ncol=50,byrow = TRUE)

## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
cigarFreq2 <- apply(cigarMat2,2,table)
cigarFreq2[1:2,]

cigarFreq3 <- apply(cigarMat3,2,table)
cigarFreq4 <- apply(cigarMat4,2,table)
cigarFreq5 <- apply(cigarMat5,2,table)
cigarFreq6 <- apply(cigarMat6,2,table)




## ---- echo=TRUE,eval=TRUE,fig.height=3,fig.width=8----------------------------
require(ggplot2)
toPlot2d <- data.frame(Freq=c(cigarFreq2["S",],cigarFreq2["M",]),
                      Cigar=rep(c("S","M"),each=ncol(cigarFreq2)),
                      Cycle=rep(seq(1,ncol(cigarFreq2)),2))
ggplot(toPlot2d,aes(x=Cycle,y=Freq,colour=Cigar))+geom_line()+theme_bw()+ggtitle("ToxA_R2")

toPlot3d <- data.frame(Freq=c(cigarFreq3["S",],cigarFreq3["M",]),
                       Cigar=rep(c("S","M"),each=ncol(cigarFreq3)),
                       Cycle=rep(seq(1,ncol(cigarFreq3)),2))
ggplot(toPlot3d,aes(x=Cycle,y=Freq,colour=Cigar))+geom_line()+theme_bw()+ggtitle("ToxB_R1")

##toPlot4d <- data.frame(Freq=c(cigarFreq4["S",],cigarFreq4["M",]),
                       ##Cigar=rep(c("S","M"),each=ncol(cigarFreq4)),
                       ##Cycle=rep(seq(1,ncol(cigarFreq4)),2))
##ggplot(toPlot4d,aes(x=Cycle,y=Freq,colour=Cigar))+geom_line()+theme_bw()+ggtitle("ToxB_R2")

toPlot5d <- data.frame(Freq=c(cigarFreq5["S",],cigarFreq5["M",]),
                       Cigar=rep(c("S","M"),each=ncol(cigarFreq5)),
                       Cycle=rep(seq(1,ncol(cigarFreq5)),2))
ggplot(toPlot5d,aes(x=Cycle,y=Freq,colour=Cigar))+geom_line()+theme_bw()+ggtitle("Control_R2")

toPlot6d <- data.frame(Freq=c(cigarFreq6["S",],cigarFreq6["M",]),
                       Cigar=rep(c("S","M"),each=ncol(cigarFreq6)),
                       Cycle=rep(seq(1,ncol(cigarFreq6)),2))
ggplot(toPlot6d,aes(x=Cycle,y=Freq,colour=Cigar))+geom_line()+theme_bw()+ggtitle("Control_R1")


## ---- echo=TRUE,eval=TRUE-----------------------------------------------------
counts1 <- data.frame(table(seqnames(temp1)),row.names = "Var1")
counts1[1:2,,drop=FALSE]

counts2 <- data.frame(table(seqnames(temp2)),row.names = "Var1")
counts2[1:2,,drop=FALSE]

counts3 <- data.frame(table(seqnames(temp3)),row.names = "Var1")
counts3[1:2,,drop=FALSE]

counts4 <- data.frame(table(seqnames(temp5)),row.names = "Var1")
counts4[1:2,,drop=FALSE]

counts5 <- data.frame(table(seqnames(temp5)),row.names = "Var1")
counts5[1:2,,drop=FALSE]

counts6 <- data.frame(table(seqnames(temp6)),row.names = "Var1")
counts6[1:2,,drop=FALSE]




## ---- echo=TRUE,eval=TRUE,fig.height=3,fig.width=8----------------------------
setwd("C:/Users/jrichards/Downloads/New_Example_1580342908_4/Analysis/01_Alignment_Results/ReadCounts_per_sgRNA")

ss1 <- read.delim("ToxA_1_GuideCounts.txt")
new1 <- merge(ss1,counts1,by.x=1,by.y=0)
##plot(new1$X0,new1$Freq)

ss2 <- read.delim("ToxA_2_GuideCounts.txt")
new2 <- merge(ss2,counts2,by.x=1,by.y=0)
##plot(new2$X0,new2$Freq)


ss3 <- read.delim("ToxB_1_GuideCounts.txt")
new3 <- merge(ss3,counts3,by.x=1,by.y=0)
##plot(new3$X0,new3$Freq)


ss4 <- read.delim("ToxB_2_GuideCounts.txt")
new1 <- merge(ss1,counts1,by.x=1,by.y=0)
##plot(new1$X0,new1$Freq)

ss5 <- read.delim("Control_1_GuideCounts.txt")
new5 <- merge(ss5,counts5,by.x=1,by.y=0)
plot(new5$X0,new5$Freq)


ss6 <- read.delim("Control_2_GuideCounts.txt")
new6 <- merge(ss6,counts6,by.x=1,by.y=0)
##plot(new6$X0,new6$Freq)

## ---- echo=TRUE,eval=TRUE,fig.height=3,fig.width=8----------------------------

counts5b <- data.frame(table(seqnames(temp5[width(temp5) == 20])),row.names = "Var1")
new5b <- merge(ss5,counts5b,by.x=1,by.y=0)
plot(new5b$X0,new5b$Freq)


## ---- echo=TRUE---------------------------------------------------------------
setwd("C:/Users/jrichards/Downloads")
myData <- read.delim("PinAPL-py_demo_data/GeCKOv21_Human.tsv",sep="\t")
require(Biostrings)
sgRNALib <- DNAStringSet(myData$seq)
names(sgRNALib) <- myData$UID
require(Rsubread)
writeXStringSet(sgRNALib,"sgRNALib.fa")
buildindex("sgRNALib",reference = "sgRNALib.fa")
myFQs <- c("PinAPL-py_demo_data/Control_R1_S14_L008_R1_001_x.fastq.gz", "PinAPL-py_demo_data/Control_R2_S15_L008_R1_001_x.fastq.gz", "PinAPL-py_demo_data/ToxA_R1_98_S2_L008_R1_001_x.fastq.gz", "PinAPL-py_demo_data/ToxA_R2_S2_L005_R1_001_x.fastq.gz", "PinAPL-py_demo_data/ToxB_R1_98_S4_L008_R1_001_x.fastq.gz", "PinAPL-py_demo_data/ToxB_R2_S4_L005_R1_001_x.fastq.gz")
require(GenomicAlignments)
counts <- list()
stats <- list()

for(f in 1:length(myFQs)){
  stats[[f]] <- align("sgRNALib",myFQs[f],output_file = gsub(".fastq.gz",".bam",myFQs[f]),
                      nthreads=2,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1,maxMismatches = 0,indels = 0)
  temp <- readGAlignments(gsub(".fastq.gz",".bam",myFQs[f]))
  counts[[f]] <- data.frame(table(seqnames(temp[width(temp) == "20"])),row.names = "Var1")
}
myRes <- do.call(cbind,counts)
colnames(myRes) <- c("Control_1","Control_2","ToxA_1","ToxA_2","ToxB_1","ToxB_2")


## ---- echo=TRUE---------------------------------------------------------------
library(DESeq2)
metadata <- DataFrame(Group=factor(c("Control","Control","ToxA","ToxA","ToxB","ToxB"),
                                   levels = c("Control","ToxA","ToxB")),
                      row.names = colnames(myRes))
dds <- DESeqDataSetFromMatrix(myRes,colData = metadata,design = ~Group)


## ---- echo=TRUE---------------------------------------------------------------
dds <- DESeq(dds)


## ---- echo=TRUE,fig.height=3,fig.width=8--------------------------------------
normCounts <- counts(dds,normalized=TRUE)
boxplot(log2(normCounts+0.1))


## ---- echo=TRUE---------------------------------------------------------------
ToxBvsControl <- results(dds,contrast=c("Group","ToxB","Control"))
ToxBvsControl <- ToxBvsControl[order(ToxBvsControl$pvalue),]
ToxBvsControl

ToxAvsControl <- results(dds,contrast=c("Group","ToxA","Control"))
ToxAvsControl <- ToxAvsControl[order(ToxAvsControl$pvalue),]
ToxAvsControl

## ---- echo=TRUE,fig.height=3,fig.width=8--------------------------------------
ssb <- read.delim("New_Example_1580342908_4/Analysis/02_sgRNA-Ranking_Results/sgRNA_Rankings/ToxB_avg_0.01_Sidak_sgRNAList.txt")
toPlt1 <- merge(ssb,as.data.frame(ToxBvsControl),by.x=1,by.y=0)
corr1 <- cor(log2(toPlt1[!is.na(toPlt1$padj),]$fold.change),toPlt1[!is.na(toPlt1$padj),]$log2FoldChange)
plot(log2(toPlt1[!is.na(toPlt1$padj),]$fold.change),toPlt1[!is.na(toPlt1$padj),]$log2FoldChange,main=corr1)

ssc <- read.delim("New_Example_1580342908_4/Analysis/02_sgRNA-Ranking_Results/sgRNA_Rankings/ToxA_avg_0.01_Sidak_sgRNAList.txt")
toPlt2 <- merge(ssc,as.data.frame(ToxAvsControl),by.x=1,by.y=0)
corr2 <- cor(log2(toPlt2[!is.na(toPlt2$padj),]$fold.change),toPlt2[!is.na(toPlt2$padj),]$log2FoldChange)
plot(log2(toPlt2[!is.na(toPlt2$padj),]$fold.change),toPlt2[!is.na(toPlt2$padj),]$log2FoldChange,main=corr2)


## -----------------------------------------------------------------------------
ToxBvsControl2 <- as.data.frame(ToxBvsControl)[order(ToxBvsControl$pvalue),]
ToxBvsControl2$Enriched <- !is.na(ToxBvsControl$padj) & ToxBvsControl$pvalue < 0.05 &ToxBvsControl$log2FoldChange > 0
ToxBvsControl2[1:2,]

ToxAvsControl2 <- as.data.frame(ToxAvsControl)[order(ToxAvsControl$pvalue),]
ToxAvsControl2$Enriched <- !is.na(ToxAvsControl$padj) & ToxAvsControl$pvalue < 0.05 &ToxAvsControl$log2FoldChange > 0
ToxAvsControl2[1:2,]


## -----------------------------------------------------------------------------
ToxBvsControl2 <- merge(GeCKO,ToxBvsControl2,by.x=2,by.y=0)
ToxBvsControl2[1:2,]

ToxAvsControl2 <- merge(GeCKO,ToxAvsControl2,by.x=2,by.y=0)
ToxAvsControl2[1:2,]


## ----eval=FALSE---------------------------------------------------------------
 genesB <- unique(ToxBvsControl2$gene_id)
 listofGene <- list()
 for(i in 1:length(genesB)){
   tempRes <- ToxBvsControl2[ToxBvsControl2$gene_id %in% genesB[i],]
   meanLogFC <- mean(tempRes$log2FoldChange,na.rm=TRUE)
   logFCs <- paste0(tempRes$log2FoldChange,collapse=";")
   minPvalue <- min(tempRes$pvalue,na.rm=TRUE)
   pvalues <- paste0(tempRes$pvalue,collapse=";")
   nEnriched <- sum(tempRes$Enriched,na.rm=TRUE)
   listofGene[[i]] <- data.frame(Gene=genesB[i],meanLogFC,logFCs,minPvalue,pvalues,nEnriched)
 }
 geneTable <- do.call(rbind,listofGene)
 geneTable <- geneTable[order(geneTable$nEnriched,decreasing=TRUE),]

 genesA <- unique(ToxAvsControl2$gene_id)
 listofGene_A <- list()
 for(i in 1:length(genesA)){
   tempRes_A <- ToxBvsControl2[ToxAvsControl2$gene_id %in% genesA[i],]
   meanLogFC_A <- mean(tempRes_A$log2FoldChange,na.rm=TRUE)
   logFCs_A <- paste0(tempRes_A$log2FoldChange,collapse=";")
   minPvalue_A <- min(tempRes_A$pvalue,na.rm=TRUE)
   pvalues_A <- paste0(tempRes_A$pvalue,collapse=";")
   nEnriched_A <- sum(tempRes_A$Enriched,na.rm=TRUE)
   listofGene_A[[i]] <- data.frame(Gene=genesA[i],meanLogFC_A,logFCs_A,minPvalue_A,pvalues_A,nEnriched_A)
 }
 geneTable_A <- do.call(rbind,listofGene)
 geneTable_A <- geneTable[order(geneTable$nEnriched,decreasing=TRUE),]
 

## -----------------------------------------------------------------------------
geneTable[1:3,]
 
geneTable_A[1:3,]


