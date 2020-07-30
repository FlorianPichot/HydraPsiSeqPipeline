#SCRIPT TO calculate Normalized U profile
rm(list=ls(all=TRUE)) #REMOVE ALL Variables

####### Input area #######
#extDataDir <- "Path/to/input"	//!\\ # output will be put there also

if (!requireNamespace("ggplot2", quietly = TRUE)) {install.packages("ggplot2")}
library(ggplot2)   

###### Treatment area ######
setwd(extDataDir)
#Select the Folder Name and use as Project Name
project <-sub('.*\\/', '', extDataDir) #FIND THE LAST "/" in the path and return the right string
project

####### MAIN BLOCK #######
#LIST FOLDERS TO TREAT to check if OK
DataDirs <- list.dirs(extDataDir, full.names = TRUE, recursive = FALSE)
DataDirs					#FULL path format
#INDICATE pattern for COUNTand COVERAGE  files, the list of analyzed RNAs is the same for all folders
count_files <- "^UCount5prime_" #indicate just a pattern to use
coverage_files <- "^coverage_" #indicate just a pattern to use

#z=DataDirs[1] #for tests

ReadsStat<- data.frame(matrix(ncol = 0, nrow = 0))

for (z in DataDirs) {
print(z)
  dir.create (paste0(z,"/StatGRAPHS"), showWarnings = FALSE)
  
# READ RNA LIST, create the list of RNA sequences to treat
Sample_files <- list.files(z, pattern=count_files, recursive=F, full.names = TRUE)
RNA_list<-basename(Sample_files)
type <- gsub(".csv$","",gsub("^UCount5prime_","", as.list(RNA_list)))
type #PRINT RNA list to check

ReadsStatSample<- data.frame(matrix(ncol = 1, nrow = 0))
colnames(ReadsStatSample)<-c("RNA")

#y=type[26] #for tests

for (y in type) {
Count_file    <- list.files(z, pattern=paste0(count_files,y,".csv"), recursive=F, full.names = TRUE) #".*" one ore more any characters
Coverage_file <- list.files(z, pattern=paste0(coverage_files,y,".csv"), recursive=F, full.names = TRUE)  
statRNA <-NA

if (!file.size(Count_file) == 0) {  #check if file is NOT empty, otherwise skip to next
#reading file from csv with names of colums and rows
sample<-read.table(Count_file, sep=" ", header = FALSE, colClasses = (rep("numeric",2)))
sample <- sample[,c(2:3)] #keep only col 2 and 3
sample <- na.omit(sample)
colnames(sample) <-c("position","counts")

Ncounts <- sum(sample$counts)

#READING COVERAGE FILE$
coverage<-read.csv(Coverage_file, sep="\t", header = FALSE)
coverage<-coverage[,c(1:3)] #skip last column with coverage value
colnames(coverage) <-c("RNA_name","position","base")
coverage$base<-gsub("t","T",coverage$base)  #replace eventual t in the coverage file sequence

RNA_length <- nrow(coverage)


#VERIFY reads COVERAGE if <5 skip all calculations and go to next RNA 
if (Ncounts/RNA_length >5) {

FullRNA<-merge(coverage,sample, by="position", all=TRUE) #merge coverage and counts
FullRNA<-FullRNA[-1,] #SKIP First row
FullRNA[is.na(FullRNA)] <- 0 #replace missing values by zero

statRNA[1]<-sum(FullRNA$counts) #NUMBER OF COUNTS for RNA
statRNA[2] <- nrow(FullRNA[which(FullRNA$counts == 0),])  #ZERO COVERAGE positions


#NORMALIZATION by window of 11 with non-T signals (special treatment for 5 first and 5 last positions)

FullRNA$NormMedian <- NA
FullRNA$Nwindow <- NA

for (x in 6:(nrow(FullRNA)-5)) {
 region <-FullRNA[c((x-5):(x+5)),]

 NonT_region<-region[which(region$base!="T"),]
 
 if(nrow(NonT_region)!=0) {  #AVOID BLOCKAGE if >10 T in a row...
 FullRNA$NormMedian[x]<-0.4  #place 0.4 as FullRNA$NormMedian value by default
 if (median(NonT_region[,4], na.rm = T) > 0) {FullRNA$NormMedian[x]<-median(NonT_region[,4], na.rm = T)} #if >0, replace by real value
 }
 
 FullRNA$Nwindow[x] <- nrow(NonT_region)
                              }

#TREAT FIRST AND LAST ROWs
regionS <- FullRNA[c(1:6),]

FullRNA$NormMedian[1:5]<-0.4
if(nrow(regionS[which(regionS$base!="T"),])!=0) {  #AVOID BLOCKAGE if >10 T in a row...
if (median(regionS[which(regionS$base!="T"),][,4], na.rm = T)>0) {FullRNA$NormMedian[1:5]<-median(regionS[which(regionS$base!="T"),][,4], na.rm = T)}
}

regionL <- FullRNA[c(nrow(FullRNA)-4):nrow(FullRNA),]
FullRNA$NormMedian[(nrow(FullRNA)-4):(nrow(FullRNA))] <-0.4

if(nrow(regionS[which(regionL$base!="T"),])!=0) {  #AVOID BLOCKAGE if >10 T in a row...
if (median(regionL[which(regionL$base!="T"),][,4], na.rm = T)>0) {FullRNA$NormMedian[(nrow(FullRNA)-4):(nrow(FullRNA))] <- median(regionL[which(regionL$base!="T"),][,4], na.rm = T)}
}
FullRNA$NormUcount <- FullRNA$counts/FullRNA$NormMedian

#drop non T data
Uprofile <-FullRNA[which(FullRNA$base=="T"),]
Uprofile$NormUcount[is.infinite(Uprofile$NormUcount)] <- 0 #REMOVE inf
Uprofile$NormUcount[is.nan(Uprofile$NormUcount)] <- 0 #REMOVE NAN

statRNA[3]<-sum(Uprofile$counts) #NUMBER OF COUNTS for RNA
statRNA[4] <- nrow(Uprofile[which(Uprofile$counts == 0),])  #ZERO COVERAGE positions

dataStatRNA<-data.frame(matrix(ncol = 1, nrow = 4))
colnames(dataStatRNA)<-c("stat")
row.names(dataStatRNA) <-c(paste0(y,"_total_reads"),paste0(y,"_pos_Zero_cov"),paste0(y,"_totalU_reads"), paste0(y,"_posU_Zero_cov") )
dataStatRNA$stat <- statRNA

ReadsStatSample<-rbind(ReadsStatSample,dataStatRNA )

#CREATE GRAPHICS
#DISTRIBUTION OF different values

df <- FullRNA[,c(4:7)]
library(ggplot2)

msize.x=1 #parameters used in inch
msize.y=1
pdf.width=2.5
pdf.height=2.5

p<-ggplot(df, aes(x=Nwindow)) + 
  geom_histogram(binwidth=1,color="black", fill="white")+
  labs(title=y)

pdf(paste0(z,"/StatGRAPHS/NwindowDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()

p<-ggplot(df, aes(x=NormMedian)) + 
  geom_histogram(binwidth=1,color="black", fill="lightblue")+
  labs(title=y)

pdf(paste0(z,"/StatGRAPHS/NormMedianDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()

msize.x=1 #parameters used in inch
msize.y=1
pdf.width=4
pdf.height=4
p<-ggplot(df, aes(x=NormUcount)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(title=y)
pdf(paste0(z,"/StatGRAPHS/NormUcountDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()

df <- Uprofile[,c(4:7)]
p<-ggplot(df, aes(x=Nwindow)) + 
  geom_histogram(binwidth=1,color="black", fill="white")+
  labs(title=y)
pdf(paste0(z,"/StatGRAPHS/U_NwindowDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()

p<-ggplot(df, aes(x=NormMedian)) + 
  geom_histogram(binwidth=1,color="black", fill="lightblue")+
  labs(title=y)
pdf(paste0(z,"/StatGRAPHS/U_NormMedianDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()

msize.x=1 #parameters used in inch
msize.y=1
pdf.width=4
pdf.height=4
p<-ggplot(df, aes(x=NormUcount)) + 
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  labs(title=y)
pdf(paste0(z,"/StatGRAPHS/U_NormUcountDistribution_",y,".pdf"),width=pdf.width,height=pdf.height,paper='special')
print(p)
dev.off()


#SAVE DATA TO DISK
output_file <- paste0(z,"/Norm_Uprofile_",y,".csv")
write.csv(Uprofile,file=output_file, quote =F)
output_file <- paste0(z,"/FullRNA_profile_",y,".csv")
write.csv(FullRNA,file=output_file, quote =F)
} #end for 1st if (file size)
} #end for 2nd if  (coverage good?)
} #end for RNA

colnames(ReadsStatSample) <- c(paste0(basename(z)))
output_file <- paste0(z,"/Stat_",basename(z),".csv")
write.csv(ReadsStatSample,file=output_file, quote =F)

ReadsStat <-merge(ReadsStat,ReadsStatSample, by="row.names", all=T)
row.names(ReadsStat)<-ReadsStat$Row.names
ReadsStat$Row.names <-NULL
}  #end for Dirnames

output_file <- paste0(extDataDir,"/GlobalStat.csv")
write.csv(ReadsStat,file=output_file, quote =F)

#END


