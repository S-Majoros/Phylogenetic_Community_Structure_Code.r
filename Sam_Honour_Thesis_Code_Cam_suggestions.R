##pipeline for my honour's thesis 
#Centroid and alignment code adapted from Matthew Orton https://github.com/m-orton/Evolutionary-Rates-Analysis-Pipeline/blob/master/EvolutionaryComparisonPipelineSmallTaxa.R
#Sequence trimming and outlier removal functions adapted from Jacqueline May https://github.com/jmay29/phylo/blob/master/refSeqTrim.R

#Part 1: Inputting and Filtering Data ----

#Packages
#install.packages("readr")
library(readr)
#install.packages("plyr")
library(plyr)
#install.packages("dplyr")
library(dplyr)
#install.packages("foreach")
library(foreach)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("stringr")
library(stringr)
#install.packages("stringi")
library(stringi)
#install.packages("ape")
library(ape)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(Biostrings)
#TODO:
#source("https://bioconductor.org/biocLite.R")
#biocLite("muscle")
library(muscle)
#install.packages("phangorn")
library(phangorn)
#install.packages("picante")
library(picante)
#install.packages("data.table")
library(data.table)
#install.packages("ggtree")
library(ggtree)
#install.packages("phytools")
library(phytools)

#Upload order data into R
#dfOrder <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Coleoptera&geo=Alaska|Canada&format=tsv")
#Write file to hard disk
#write_tsv(dfOrder, "Coleoptera_download_Oct26")
#Read in order
#dfOrder <- read_tsv("Coleoptera_download_Oct26")

# me replicating file using BOLD-CLI
dfOrder = read_tsv('majoros_bold_data.txt')
#loads, with some column erros?

#Filtering the data
dfOrder <- dfOrder %>%
  #Filter out those without bin_uri
  filter(str_detect(bin_uri, ":")) %>%
  #Filter out those without a sequence 
  filter(str_detect(nucleotides, "[ACTG]")) %>%
  #Filter for COI-5P
  filter(markercode == "COI-5P") %>%
  #Filter out sequences with fewer than 500 base pairs 
  filter(nchar(gsub("-", "", nucleotides)) > 499) %>%
  #Filter out records without a family name
  filter(!is.na(family_name))

#Filter out high gap/N content. A threshold of 1% was chosen because species often differ by more than 2% divergence. By filtering out records with > 1% N and gap content, we are likely to get a high-quality data set, given typical patterns of variability in COI in animals. 
startNGap <- sapply(regmatches(dfOrder$nucleotides, gregexpr("^[-N]", dfOrder$nucleotides)), length)
startNGap <- foreach(i=1:nrow(dfOrder)) %do% 
  if (startNGap[[i]]>0) { 
    split <- strsplit(dfOrder$nucleotides[i], "^[-N]+") 
    dfOrder$nucleotides[i] <- split[[1]][2]
  }
endNGap <- sapply(regmatches(dfOrder$nucleotides, gregexpr("[-N]$", dfOrder$nucleotides)), length)
endNGap <- foreach(i=1:nrow(dfOrder)) %do%
  if (endNGap[[i]]>0) {
    split <- strsplit(dfOrder$nucleotides[i], "[-N]+$")
    dfOrder$nucleotides[i] <- split[[1]][1]
  }
internalNGap <- sapply(regmatches(dfOrder$nucleotides, gregexpr("[-N]", dfOrder$nucleotides)), length)
internalNGap <- foreach(i=1:nrow(dfOrder)) %do%
  which((internalNGap[[i]]/nchar(dfOrder$nucleotides[i]) > 0.01))
nGapCheck <- sapply(internalNGap, function(x)length(x))
nGapCheck <- which(nGapCheck>0)
dfOrder <- dfOrder[-nGapCheck, ]
#Remove redundant "BOLD" section from BIN column 
dfOrder$bin_uri <- substr(dfOrder$bin_uri, 6, 13)
#Filter out sequences without coordinates 
containLatLon <- grep ("[0-9]", dfOrder$lat)
dfOrder <- dfOrder[containLatLon, ]

#Create subset filter using coordinates 
#Filter for Churchill, Manitoba 
SubsetFilter_Churchill <- which(dfOrder$lat > 58.6 &
                           dfOrder$lon > -94.2 & dfOrder$lat < 58.7 &
                           dfOrder$lon < -93.8)
#Apply filter
dfOrder_Churchill <- dfOrder[SubsetFilter_Churchill, ]

#Find total number of BINs per family in the regional subset of the order
dfOrder_Churchill <- as.data.table(dfOrder_Churchill)
total_species_number <- dfOrder_Churchill[ , .(.N),by=.(family_name)]
number_of_unique_species <- dfOrder_Churchill[ , .(number_of_species=length(unique(bin_uri))), by=family_name]
number_of_unique_species <- as.data.frame(number_of_unique_species)
#Filter down to families with more than 3 or more species
number_of_unique_species <- filter(number_of_unique_species, number_of_unique_species$number_of_species > 2)

#Create filter to filter down the order to families with three or more species in the subset  
dfOrder_filter <- which(dfOrder$family_name %in% number_of_unique_species$family_name)
#Apply filter 
dfOrder <- dfOrder[dfOrder_filter, ]

#Remove uneeded variables 
rm(number_of_unique_species, total_species_number, SubsetFilter_Churchill, containLatLon, endNGap, internalNGap, nGapCheck, split, startNGap, dfOrder_filter, i)

#Part 2: Choosing a Centroid----

#In this section we find a centroid sequence for each BIN present in the order (Not the subset)
#Create smaller dataframe with needed info
dfBinList <- (dfOrder[, c("processid", "bin_uri", "nucleotides")])
#Create groupings by BIN, each with different bin_uri
binList <- lapply(unique(dfOrder$bin_uri), function(x) dfOrder[dfOrder$bin_uri==x, ])
#Number of processids in each bin 
binSize <- sapply(binList, function(x)length(x$processid))
#Create new data frame with bin_uri and bin size
dfOrder_bins <- data.frame(binSize)
dfOrder_bins$bin_uri <- c(unique(dfOrder$bin_uri))
#Merge dfBinList and dfOrder_bins
dfBinList <- merge(dfBinList, dfOrder_bins, by.x="bin_uri", by.y="bin_uri")
#Reorder dfFamily_bins by bin_uri
dfOrder_bins <- dfOrder_bins[order(dfOrder_bins$bin_uri), ]

#Find bins with more than one member
largeBin <- which(dfBinList$binSize > 1)   
#Create dataframe with only bins with more than one member 
if (length(largeBin) > 0) {
  dfCentroid <- dfBinList[largeBin, ]
}

#Subset dfOrder_bins down to number of bins in dfOrder
dfOrder_bins <- subset(dfOrder_bins, dfOrder_bins$bin_uri %in% dfCentroid$bin_uri)

#Find number of unique bins in dfCentroid 
binNumberCentroid <- unique(dfCentroid$bin_uri)   
binNumberCentroid <- length(binNumberCentroid) 

#Create dataframe with bins with only one sequence 
dfNonCentroid <- dfBinList[-largeBin, ]

#Create list from dfCentroid 
largeBinList <- lapply(unique(dfCentroid$bin_uri), function(x) dfCentroid[dfCentroid$bin_uri == x, ])
#Extract process Id from each bin
largeBinProcessid <- sapply(largeBinList, function(x) (x$processid))

#Convert sequences to dnaStringSet
dnaStringSet1 <- sapply(largeBinList, function(x) DNAStringSet(x$nucleotides))  
#Name dnaStringSet with processids   
for(i in seq(from=1, to=binNumberCentroid, by=1)) {
  names(dnaStringSet1[[i]]) <- largeBinProcessid[[i]]
}           

#Run multiple sequence alignment for sequences in each BIN in dnaStringSet1   
alignment1 <- foreach(i=1:binNumberCentroid) %do%
  muscle::muscle(dnaStringSet1[[i]], maxiters=3, diags=TRUE, gapopen=-3000)

#Convert to DNAbin format 
dnaBINCentroid <- foreach(i=1:binNumberCentroid) %do% as.DNAbin(alignment1[[i]])

#Calculate a pairwise distance matrix for each BIN
geneticDistanceCentroid <- foreach(i=1:binNumberCentroid) %do%
  dist.dna(dnaBINCentroid[[i]], model="TN93", as.matrix = TRUE,
           pairwise.deletion = TRUE)

#Determine centroid sequence; The sequence with the minimum average distance to all other sequences in the BIN.
centroidSeq <- foreach(i=1:binNumberCentroid) %do% which.min(rowSums(geneticDistanceCentroid[[i]]))
centroidSeq <- centroidSeq %>%
  unlist() %>%
  names()

#Subset dfCentroid by the processid on the list 
dfCentroid <- subset(dfCentroid, processid %in% centroidSeq)

#Merge with dfNonCentroid
dfAllSeq <- rbind(dfCentroid, dfNonCentroid)
#Merge with the original data set 
dfAllSeq <- merge(dfAllSeq, dfOrder, by.x="processid", by.y="processid")
#Reorganize and clean up 
dfAllSeq <- (dfAllSeq[, c("bin_uri.x", "binSize", "processid", "family_taxID", "family_name", "species_taxID", "species_name", "nucleotides.x", "lat", "lon", "subfamily_name", "order_name")])
colnames(dfAllSeq)[1] <- "bin_uri"
colnames(dfAllSeq)[8] <- "nucleotides"
#Delete any possible duplicate entries 
dfAllSeq <- (by(dfAllSeq, dfAllSeq["bin_uri"], head, n=1))
dfAllSeq <- Reduce(rbind, dfAllSeq)
#Add an index column
dfAllSeq$ind <- row.names(dfAllSeq)

#Remove unneeded dataframes and variables 
rm(alignment1, binList, dfBinList, dfCentroid, dfOrder_bins, dfNonCentroid, dnaBINCentroid, dnaStringSet1, geneticDistanceCentroid, largeBinList, largeBinProcessid, binNumberCentroid, binSize,  centroidSeq, i, largeBin)

#Part 3: Alignment----

#Create a function to trim the sequences 
RefSeqTrim <- function(x) {
  #Create data frame for reference sequence 
  #This reference sequence was taken from BOLD for Coleoptera. Process id: AEDNA549-12. Species:Colymbetes dolabratus.
  dfRefSeq <- data.frame(taxa=c("Coleoptera"), nucleotides=c("TAACTTTATATTTTATTTTTGGTGCATGGGCTGGAATGGTAGGAACATCTTTAAGTATGTTGATTCGAGCCGAATTAGGAAATCCTGGTTCTCTGATTGGAGATGATCAAATTTATAATGTTATTGTAACAGCACATGCTTTTGTAATAATTTTTTTCATAGTAATACCTATTATAATTGGGGGATTTGGAAATTGATTAGTTCCATTAATATTGGGGGCCCCAGATATAGCTTTTCCCCGAATAAATAATATAAGTTTTTGACTTCTTCCGCCTTCTTTAACTCTTCTATTAATAAGAAGAATAGTTGAAAGTGGGGCCGGGACAGGATGAACAGTTTACCCCCCTCTATCTTCAGGAATTGCACACGGAGGAGCTTCAGTTGATCTAGCAATTTTTAGTCTTCATTTAGCTGGAATTTCATCTATTTTAGGGGCTGTAAATTTCATTACAACTATTATTAATATACGATCAGTGGGAATAACATTCGACCGAATGCCTCTATTTGTATGATCCGTAGGAATTACAGCTTTATTACTATTATTATCTTTACCTGTATTAGCGGGAGCTATTACTATATTATTAACTGATCGTAATCTAAACACCTCATTCTTCGACCCGGCAGGAGGGGGAGATCCAATTTTATATCAACATTTATT"))
  colnames(dfRefSeq)[2] <- "nucleotides"
  #Convert to datatable
  dfRefSeq <- setDT(dfRefSeq)
  dfRefSeq[, "nucleotides":=as.character(nucleotides)]
  #Trim sequences to 620bp
  dfRefSeq[, nucleotides:=substr(nucleotides, 20, nchar(nucleotides)-19)]
  #Check sequence length
  dfRefSeq[, seqLength:=nchar(nucleotides)]
  #Ensure sequences are of character type
  alignmentSeqs <- as.character(x$nucleotides)
  #Name according to bin_uri
  names(alignmentSeqs) <- x$bin_uri
  alignmentref <- as.character(dfRefSeq$nucleotides[1])
  #Name reference sequence 
  names(alignmentref) <- "Reference"
  #Put sequences together
  alignmentSeqsPlusRef <- append(alignmentref, alignmentSeqs)
  #Convert to DNAStringSet 
  DNAStringSet2 <- DNAStringSet(alignmentSeqsPlusRef)
  #Run alignment 
  alignment2 <- muscle::muscle(DNAStringSet2, diags=TRUE, gapopen=-3000)
  #Check alignment 
  classFileNames <- foreach(i=1:nrow(dfRefSeq)) %do% 
    paste("alignmentUntrimmed", dfRefSeq$taxa[i], ".fas", sep="")
  alignmentUntrimmed <- DNAStringSet(alignment2)
  writeXStringSet(alignmentUntrimmed, file=classFileNames[[1]],
                  format = "fasta", width=1500)
  #Find stop and start postions in reference 
  refSeqPos <- which(alignment2@unmasked@ranges@NAMES=="Reference")
  refSeqPos <- alignment2@unmasked[refSeqPos]
  refSeqPosStart <- regexpr("[ACTG]", refSeqPos)
  refSeqPosStart <- as.numeric(refSeqPosStart)
  refSeqPosEnd <- nchar(dfRefSeq$nucleotides[1]) + refSeqPosStart
  refSeqPosEnd <- as.numeric(refSeqPosEnd)
  #Trim sequence
  alignment2Trimmed <- substr(alignment2, refSeqPosStart, refSeqPosEnd)
  #Convert to DNAStringSet
  DNAStringSet3 <- DNAStringSet(alignment2Trimmed)
  #Check alignment 
  classFileNames <- foreach(i=1:nrow(dfRefSeq)) %do% 
    paste("alignmentTrimmed", dfRefSeq$taxa[i], ".fas", sep="")
  writeXStringSet(DNAStringSet3, file=classFileNames[[1]],
                  format = "fasta", width=1500)
  #Remove reference sequence 
  refSeqRm <- which(DNAStringSet3@ranges@NAMES=="Reference")
  dnaStringSet3 <- subset(DNAStringSet3[-refSeqRm])
  alignmentOrder <- DNAStringSet3@ranges@NAMES
  #Reorder based on alignment
  x <- x[match(alignmentOrder, x$bin_uri), ]
  #Replace old sequences with new ones 
  trimmedSeqs <- as.character(DNAStringSet3)
  x$nucleotides <- trimmedSeqs
  #Return datafrmae with new sequences 
  return(x)
}

#Trim centroid sequences to reference sequence 
dfAllSeq2 <- RefSeqTrim(dfAllSeq)

#Convert sequences to DNAbin format 
DNABin <- DNAStringSet(dfAllSeq2$nucleotides)
names(DNABin) <- dfAllSeq2$bin_uri
DNABin <- as.DNAbin(DNABin)
#Construct a distance matrix 
distanceMatrix <- dist.dna(DNABin, model="TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
#Visualizing the values in the distance matrix using a histogram
hist(distanceMatrix)

#Using upper threshold of IQR to detect outliers 
lowerQuantile <- quantile(distanceMatrix)[2]
upperQuantile <- quantile(distanceMatrix)[4]
iqr <-upperQuantile - lowerQuantile
#For now leaving as 1.5, could change to 3 for extreme outliers only 
upperThreshold <- (iqr*1.5) + upperQuantile
#Remove 0 values 
distanceMatrix[distanceMatrix==0] <- NA
#Convert to data table 
dfOutliers <- as.data.table(distanceMatrix, keep.rownames = T)
#Change the "rn" column to bin_uri
setnames(dfOutliers, "rn", "bin_uri")
#Identify divergent bins 
dfOutliers <- dfOutliers[, outlier := apply(.SD, 1, function(x)all(x>upperThreshold, na.rm=T))][outlier==TRUE]

#Create remove sequences function
RemoveSequences<-function(x, y){
  if(length(y)==0){
    print("There are no sequences to remove!")
  }
  else if(length(y)>0){
    x <- x[!x$bin_uri%in%y]
  }
  return(x)
}

#Remove outliers 
#Outliers should be blasted prior to removal 
dfAllSeq <- RemoveSequences(dfAllSeq, dfOutliers$bin_uri)

#Create final alignment of sequences 
#Create RefSeq data frame
#Sequence was taken from BOLD and manually put in
dfRefSeq <- data.frame(taxa=c("Coleoptera"), nucleotides=c("TAACTTTATATTTTATTTTTGGTGCATGGGCTGGAATGGTAGGAACATCTTTAAGTATGTTGATTCGAGCCGAATTAGGAAATCCTGGTTCTCTGATTGGAGATGATCAAATTTATAATGTTATTGTAACAGCACATGCTTTTGTAATAATTTTTTTCATAGTAATACCTATTATAATTGGGGGATTTGGAAATTGATTAGTTCCATTAATATTGGGGGCCCCAGATATAGCTTTTCCCCGAATAAATAATATAAGTTTTTGACTTCTTCCGCCTTCTTTAACTCTTCTATTAATAAGAAGAATAGTTGAAAGTGGGGCCGGGACAGGATGAACAGTTTACCCCCCTCTATCTTCAGGAATTGCACACGGAGGAGCTTCAGTTGATCTAGCAATTTTTAGTCTTCATTTAGCTGGAATTTCATCTATTTTAGGGGCTGTAAATTTCATTACAACTATTATTAATATACGATCAGTGGGAATAACATTCGACCGAATGCCTCTATTTGTATGATCCGTAGGAATTACAGCTTTATTACTATTATTATCTTTACCTGTATTAGCGGGAGCTATTACTATATTATTAACTGATCGTAATCTAAACACCTCATTCTTCGACCCGGCAGGAGGGGGAGATCCAATTTTATATCAACATTTATT"))

colnames(dfRefSeq)[2] <- "nucleotides"
dfRefSeq$nucleotides <- as.character(dfRefSeq$nucleotides) 
#Trim references to standard 620
dfRefSeq$nucleotides <- substr(dfRefSeq$nucleotides, 20, nchar(dfRefSeq$nucleotides)-19) 
#Check sequence length
dfRefSeq$seqLength <- nchar(dfRefSeq$nucleotides)
#Subset centroid sequences by those found in reference sequence dataframe
dfAllSeq <- subset(dfAllSeq, dfAllSeq$order_name %in% dfRefSeq$taxa) 
#Break down dataframe into families
taxalistcomplete <- lapply(unique(dfAllSeq$family_taxID), function(x) dfAllSeq[dfAllSeq$family_taxID==x, ])

#Extract sequences and bin_uri
familyBin <- foreach(i=1:length(taxalistcomplete)) %do% taxalistcomplete[[i]]$bin_uri
familySequences <- foreach(i=1:length(taxalistcomplete)) %do% taxalistcomplete[[i]]$nucleotides
familySequenceNames <- familyBin

#Take reference sequences
alignmentref <- as.character(dfRefSeq$nucleotides)
dfRefSeq$reference <- "reference"
#Name reference as a reference
alignmentRefNames <- dfRefSeq$reference
#Merge reference with other sequences
alignmentSequencesPlusRef <- foreach(i=1:length(taxalistcomplete)) %do%
  append(familySequences[[i]], alignmentref[[1]])

#Merge names together
alignmentNames <- foreach(i=1:length(taxalistcomplete)) %do%
  append(familySequenceNames[[i]], alignmentRefNames[[1]])

#Convert sequences to DNAStringSet format 
dnaStringSet3 <- foreach(i=1:length(alignmentSequencesPlusRef)) %do%
  DNAStringSet(alignmentSequencesPlusRef[[i]])

#Then name each one individually for each family
names(dnaStringSet3[[1]]) <- alignmentNames[[1]]
names(dnaStringSet3[[2]]) <- alignmentNames[[2]]
names(dnaStringSet3[[3]]) <- alignmentNames[[3]]
names(dnaStringSet3[[4]]) <- alignmentNames[[4]]
names(dnaStringSet3[[5]]) <- alignmentNames[[5]]
names(dnaStringSet3[[6]]) <- alignmentNames[[6]]
names(dnaStringSet3[[7]]) <- alignmentNames[[7]]
names(dnaStringSet3[[8]]) <- alignmentNames[[8]]
names(dnaStringSet3[[9]]) <- alignmentNames[[9]]
names(dnaStringSet3[[10]]) <- alignmentNames[[10]]
names(dnaStringSet3[[11]]) <- alignmentNames[[11]]
names(dnaStringSet3[[12]]) <- alignmentNames[[12]]
names(dnaStringSet3[[13]]) <- alignmentNames[[13]]
names(dnaStringSet3[[14]]) <- alignmentNames[[14]]
names(dnaStringSet3[[15]]) <- alignmentNames[[15]]
names(dnaStringSet3[[16]]) <- alignmentNames[[16]]


############################
#Cam's refactoring of above change to dataframes

#note this can be done with an lapply like my solutions before (could be a take home challenge)
#but since they both rely on the index, we can keep it simple.

for(i in 1:16){
	names(dnaStringSet3[[i]]) = alignmentNames[[i]] 
}


###########################



#Multiple sequence alignment 
alignmentFinal <- foreach(i=1:length(dnaStringSet3)) %do%
  muscle(dnaStringSet3[[i]], diags=TRUE, gapopen=-3000)
#Check Alignment
familyFileNames2 <- foreach(i=1:length(alignmentFinal)) %do%
  paste("alignmentFinal", dfRefSeq$taxa[i], ".fas", sep="")
alignmentFinalFasta <- foreach(i=1:length(alignmentFinal)) %do%
  DNAStringSet(alignmentFinal[[i]])
foreach(i=1:length(alignmentFinal)) %do%
  writeXStringSet(alignmentFinalFasta[[i]], file=familyFileNames2[[i]], format="fasta", width=1500)

#Convert to dnaStringSet format
dnaStringSet4 <- foreach(i=1:length(alignmentFinal)) %do%
  DNAStringSet(alignmentFinal[[i]])

#Remove unneeded info
rm(alignmentFinal, alignmentNames, alignmentSequencesPlusRef, dnaStringSet3, familyBin, alignmentref, alignmentRefNames, i, dfAllSeq2, dfOutliers, distanceMatrix, DNABin, familySequences, iqr, lowerQuantile, upperQuantile, upperThreshold, dfRefSeq, familySequenceNames)

#Part 4: Create Maximum Liklihood tree----

#Would first need to make a separate dataframe for each family and then rbind into one dataframe
FamilyDNA1 <- as.data.frame(dnaStringSet4[[1]])
FamilyDNA2 <- as.data.frame(dnaStringSet4[[2]])
FamilyDNA3 <- as.data.frame(dnaStringSet4[[3]])
FamilyDNA4 <- as.data.frame(dnaStringSet4[[4]])
FamilyDNA5 <- as.data.frame(dnaStringSet4[[5]])
FamilyDNA6 <- as.data.frame(dnaStringSet4[[6]])
FamilyDNA7 <- as.data.frame(dnaStringSet4[[7]])
FamilyDNA8 <- as.data.frame(dnaStringSet4[[8]])
FamilyDNA9 <- as.data.frame(dnaStringSet4[[9]])
FamilyDNA10 <- as.data.frame(dnaStringSet4[[10]])
FamilyDNA11 <- as.data.frame(dnaStringSet4[[11]])
FamilyDNA12 <- as.data.frame(dnaStringSet4[[12]])
FamilyDNA13 <- as.data.frame(dnaStringSet4[[13]])
FamilyDNA14 <- as.data.frame(dnaStringSet4[[14]])
FamilyDNA15 <- as.data.frame(dnaStringSet4[[15]])
FamilyDNA16 <- as.data.frame(dnaStringSet4[[16]])
FamilyDNA <- rbind(FamilyDNA1, FamilyDNA2, FamilyDNA3, FamilyDNA4, FamilyDNA5, FamilyDNA6, FamilyDNA7, FamilyDNA8, FamilyDNA9, FamilyDNA10, FamilyDNA11, FamilyDNA12, FamilyDNA13, FamilyDNA14, FamilyDNA15, FamilyDNA16)

############################
#Cam's refactoring of above change to dataframes

FamilyDNA = lapply(dnaStringSet4, as.data.frame)


###########################



#Add the bin_uri
FamilyDNA$bin_uri <- row.names(FamilyDNA)
#Merge with the information for dfAllSeq
dfFamilyDNA <- merge(FamilyDNA, dfAllSeq, by.x = "bin_uri", by.y = "bin_uri", all.x = TRUE)
#Rename the coloumn with your aligned sequences
colnames(dfFamilyDNA)[2] <- "FinalSequences"

#Remove the reference sequences
referencefind1 <- which(dfFamilyDNA$bin_uri == "reference")
dfFamilyDNA <- dfFamilyDNA[-referencefind1, ]
referencefind2 <- which(dfFamilyDNA$bin_uri == "reference1")
dfFamilyDNA <- dfFamilyDNA[-referencefind2, ]
referencefind3 <- which(dfFamilyDNA$bin_uri == "reference2")
dfFamilyDNA <- dfFamilyDNA[-referencefind3, ]
referencefind4 <- which(dfFamilyDNA$bin_uri == "reference3")
dfFamilyDNA <- dfFamilyDNA[-referencefind4, ]
referencefind5 <- which(dfFamilyDNA$bin_uri == "reference4")
dfFamilyDNA <- dfFamilyDNA[-referencefind5, ]
referencefind6 <- which(dfFamilyDNA$bin_uri == "reference5")
dfFamilyDNA <- dfFamilyDNA[-referencefind6, ]
referencefind7 <- which(dfFamilyDNA$bin_uri == "reference6")
dfFamilyDNA <- dfFamilyDNA[-referencefind7, ]
referencefind8 <- which(dfFamilyDNA$bin_uri == "reference7")
dfFamilyDNA <- dfFamilyDNA[-referencefind8, ]
referencefind9 <- which(dfFamilyDNA$bin_uri == "reference8")
dfFamilyDNA <- dfFamilyDNA[-referencefind9, ]
referencefind10 <- which(dfFamilyDNA$bin_uri == "reference9")
dfFamilyDNA <- dfFamilyDNA[-referencefind10, ]
referencefind11 <- which(dfFamilyDNA$bin_uri == "reference10")
dfFamilyDNA <- dfFamilyDNA[-referencefind11, ]
referencefind12 <- which(dfFamilyDNA$bin_uri == "reference11")
dfFamilyDNA <- dfFamilyDNA[-referencefind12, ]
referencefind13 <- which(dfFamilyDNA$bin_uri == "reference12")
dfFamilyDNA <- dfFamilyDNA[-referencefind13, ]
referencefind14 <- which(dfFamilyDNA$bin_uri == "reference13")
dfFamilyDNA <- dfFamilyDNA[-referencefind14, ]
referencefind15 <- which(dfFamilyDNA$bin_uri == "reference14")
dfFamilyDNA <- dfFamilyDNA[-referencefind15, ]
referencefind16 <- which(dfFamilyDNA$bin_uri == "reference15")
dfFamilyDNA <- dfFamilyDNA[-referencefind16, ]


###################################################
# Cam's refactoring of above

#this function isn't really necessary, but if you start repetitive strings with
#more then 15 instances then the need to automate arises. Plus it helps us learn about
#iteration
get_reference_names = function(top_ref_num = 15){
	#the highest number used in the references is passed in, default is 15
	ref_names = c()	
	prefix = "reference"
	for(i in 1:top_ref_num){
		new_str = paste(prefix, as.character(i), sep='')
		ref_names = c(ref_names, new_str)
	}
	return(ref_names)
}

reference_names = get_reference_names()

#simple demonstration of the %in% pattern
x = c('reference1', 'bill', 'george', 'reference2')
x[x %in% reference_names]
x[!x %in% reference_names]

#I'm not famility with the 'which' function, but this is a bit simpler 
dfFamilyDNA = dfFamilyDNA[!dfFamilyDNA$bin_uri %in% reference_names , ]

###################################################


#Pull names from dataframe
familyList <- lapply(unique(dfFamilyDNA$family_name), 
                     function(x) dfFamilyDNA[dfFamilyDNA$family_name == x, ])
#Create new dnaStringSet
dnaStringSet5 <- sapply(familyList, function(x) DNAStringSet(x$FinalSequences))
#Pull bin names from list
binNames <- sapply(familyList, function(x)(x$bin_uri))
#Name the stringsets
for(i in seq(from = 1, to = length(dnaStringSet5), by = 1)) {
  names(dnaStringSet5[[i]]) <- binNames[[i]]
}

#Save family as a fasta file
#For file names make sure to list each family name
familyFileNames <- list("Carabidae", "Curculionidae", "Dytiscidae", "Coccinellidae", "Leiodidae", "Chrysomelidae", "Staphylinidae", "Buprestidae", "Hydrophilidae", "Haliplidae", "Cantharidae", "Gyrinidae", "Elateridae", "Cryptophagidae", "Scirtidae", "Latridiidae")
familyFileNames <- foreach(i=1:length(familyFileNames)) %do% 
  paste("Alignment", familyFileNames[[i]], ".fas", sep="")
#Send to your desired working directory
foreach(i=1:length(dnaStringSet5)) %do% writeXStringSet(dnaStringSet5[[i]], file=familyFileNames[[i]], format="fasta")

#Create a separate function for each family to generate files, must be in same working directory you wrote to
#Double check that these are named with the proper family 
phyDatDyt <- read.phyDat("AlignmentDytiscidae.fas", format="fasta", type="DNA")
phyDatCara <- read.phyDat("AlignmentCarabidae.fas", format="fasta", type="DNA")
phyDatCur <- read.phyDat("AlignmentCurculionidae.fas", format="fasta", type="DNA")
phyDatCocc <- read.phyDat("AlignmentCoccinellidae.fas", format="fasta", type="DNA")
phyDatLei <- read.phyDat("AlignmentLeiodidae.fas", format="fasta", type="DNA")
phyDatChry <- read.phyDat("AlignmentChrysomelidae.fas", format="fasta", type="DNA")
phyDatBup <- read.phyDat("AlignmentBuprestidae.fas", format="fasta", type="DNA")
phyDatHyd <- read.phyDat("AlignmentHydrophilidae.fas", format="fasta", type="DNA")
phyDatHal <- read.phyDat("AlignmentHaliplidae.fas", format="fasta", type="DNA") 
phyDatCan <- read.phyDat("AlignmentCantharidae.fas", format="fasta", type="DNA")
phyDatGyr <- read.phyDat("AlignmentGyrinidae.fas", format="fasta", type="DNA")
phyDatEla <- read.phyDat("AlignmentElateridae.fas", format="fasta", type="DNA")
phyDatCryp <- read.phyDat("AlignmentCryptophagidae.fas", format="fasta", type="DNA")
phyDatSci <- read.phyDat("AlignmentScirtidae.fas", format="fasta", type="DNA")
phyDatLat <- read.phyDat("AlignmentLatridiidae.fas", format="fasta", type="DNA")
#phyDatStap <- read.phyDat("AlignmentStaphylinidae.fas", format="fasta", type="DNA")
#Had an error with this family, so it was excluded from the rest of the analysis. Plan to revisit at a later date.

#Create distance matrix
dm1 <- dist.ml(phyDatDyt)
dm2 <- dist.ml(phyDatCara)
dm3 <- dist.ml(phyDatCur)
dm4 <- dist.ml(phyDatCocc)
dm5 <- dist.ml(phyDatLei)
dm6 <- dist.ml(phyDatChry)
dm7 <- dist.ml(phyDatBup)
dm8 <- dist.ml(phyDatHyd)
dm9 <- dist.ml(phyDatHal)
dm10 <- dist.ml(phyDatCan)
dm11 <- dist.ml(phyDatGyr)
dm12 <- dist.ml(phyDatEla)
dm13 <- dist.ml(phyDatCryp)
dm14 <- dist.ml(phyDatSci)
dm15 <- dist.ml(phyDatLat)

#Create tree
tree1 <- fastme.bal(dm1)
tree2 <- fastme.bal(dm2)
tree3 <- fastme.bal(dm3)
tree4 <- fastme.bal(dm4)
tree5 <- fastme.bal(dm5)
tree6 <- fastme.bal(dm6)
tree7 <- fastme.bal(dm7)
tree8 <- fastme.bal(dm8)
tree9 <- fastme.bal(dm9)
tree10 <- fastme.bal(dm10)
tree11 <- fastme.bal(dm11)
tree12 <- fastme.bal(dm12)
tree13 <- fastme.bal(dm13)
tree14 <- fastme.bal(dm14)
tree15 <- fastme.bal(dm15)

#Run model test
mt1 <- modelTest(phyDatDyt)
mt2 <- modelTest(phyDatCara)
mt3 <- modelTest(phyDatCur)
mt4 <- modelTest(phyDatCocc)
mt5 <- modelTest(phyDatLei)
mt6 <- modelTest(phyDatChry)
mt7 <- modelTest(phyDatBup)
mt8 <- modelTest(phyDatHyd)
mt9 <- modelTest(phyDatHal)
mt10 <- modelTest(phyDatCan)
mt11 <- modelTest(phyDatGyr)
mt12 <- modelTest(phyDatEla)
mt13 <- modelTest(phyDatCryp)
mt14 <- modelTest(phyDatSci)
mt15 <- modelTest(phyDatLat)
#Create environemnt
env1 <-attr(mt1, "env")
env2 <- attr(mt2, "env")
env3 <- attr(mt3, "env")
env4 <- attr(mt4, "env")
env5 <- attr(mt5, "env")
env6 <- attr(mt6, "env")
env7 <- attr(mt7, "env")
env8 <- attr(mt8, "env")
env9 <- attr(mt9, "env")
env10 <- attr(mt10, "env")
env11 <- attr(mt11, "env")
env12 <- attr(mt12, "env")
env13 <- attr(mt13, "env")
env14 <- attr(mt14, "env")
env15 <- attr(mt15, "env")
#Find parameters 
fit1 <- eval(get("HKY+G+I", env1), env1)
fit2 <- eval(get("HKY+G+I", env2), env2)
fit3 <- eval(get("HKY+G+I", env3), env3)
fit4 <- eval(get("HKY+G+I", env4), env4)
fit5 <- eval(get("HKY+G+I", env5), env5)
fit6 <- eval(get("HKY+G+I", env6), env6)
fit7 <- eval(get("HKY+G+I", env7), env7)
fit8 <- eval(get("HKY+G+I", env8), env8)
fit9 <- eval(get("GTR+G+I", env9), env9)
fit10 <- eval(get("HKY+G+I", env10), env10)
fit11 <- eval(get("GTR+G+I", env11), env11)
fit12 <- eval(get("HKY+G+I", env12), env12)
fit13 <- eval(get("HKY+G+I", env13), env13)
fit14 <- eval(get("HKY+G+I", env14), env14)
fit15 <- eval(get("GTR+G+I", env15), env15)
#Compute likelihood 
ML_Dyt <- pml(tree1, phyDatDyt, k=4, inv= 0.4681883)
ML_Car <- pml(tree2, phyDatCara, k=4, inv= 0.5308776)
ML_Cur <- pml(tree3, phyDatCur, k=4, inv= 0.5482425)
ML_Cocc <- pml(tree4, phyDatCocc, k=4, inv= 0.5581981)
ML_Lei <- pml(tree5, phyDatLei, k=4, inv= 0.5581981)
ML_Chry <- pml(tree6, phyDatChry, k=4, inv= 0.528835)
ML_Bup <- pml(tree7, phyDatBup, k=4, inv= 0.5104581)
ML_Hyd <- pml(tree8, phyDatHyd, k=4, inv= 0.5886127)
ML_Hal <- pml(tree9, phyDatHal, k=4, inv= 0.6695774)
ML_Can <- pml(tree10, phyDatCan, k=4, inv= 0.5145766)
ML_Gyr <- pml(tree11, phyDatGyr, k=4, inv= 0.6361499)
ML_Ela <- pml(tree12, phyDatEla, k=4, inv= 0.5619085)
ML_Cryp <- pml(tree13, phyDatCryp, k=4, inv= 0.4869501)
ML_Sci <- pml(tree14, phyDatSci, k=4, inv= 0.5677839)
ML_Lat <- pml(tree15, phyDatLat, k=4, inv= 0.5619141)
#Compute likelihood and optimize parameters 
#Change model based on results of model test 
ML_Dyt <- optim.pml(ML_Dyt, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Car <- optim.pml(ML_Car, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Cur <- optim.pml(ML_Cur, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Cocc <- optim.pml(ML_Cocc, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Lei <- optim.pml(ML_Lei, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Chry <- optim.pml(ML_Chry, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Bup <- optim.pml(ML_Bup, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Hyd <- optim.pml(ML_Hyd, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Hal <- optim.pml(ML_Hal, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "GTR")
ML_Can <- optim.pml(ML_Can, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Gyr <- optim.pml(ML_Gyr, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "GTR")
ML_Ela <- optim.pml(ML_Ela, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Cryp <- optim.pml(ML_Cryp, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Sci <- optim.pml(ML_Sci, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Lat <- optim.pml(ML_Lat, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "GTR")

#Create seperate variable for tree 
ML_Tree_Dyt <- ML_Dyt$tree
ML_Tree_Car <- ML_Car$tree
ML_Tree_Cur <- ML_Cur$tree
ML_Tree_Cocc <- ML_Cocc$tree
ML_Tree_Lei <- ML_Lei$tree
ML_Tree_Chry <- ML_Chry$tree
ML_Tree_Bup <- ML_Bup$tree
ML_Tree_Hyd <- ML_Hyd$tree
ML_Tree_Hal <- ML_Hal$tree
ML_Tree_Can <- ML_Can$tree
ML_Tree_Gyr <- ML_Gyr$tree
ML_Tree_Ela <- ML_Ela$tree
ML_Tree_Cryp <- ML_Cryp$tree
ML_Tree_Sci <- ML_Sci$tree
ML_Tree_Lat <- ML_Lat$tree
#Plot trees
plot(ML_Tree_Dyt)
plot(ML_Tree_Car)
plot(ML_Tree_Cur)
plot(ML_Tree_Cocc)
plot(ML_Tree_Lei)
plot(ML_Tree_Chry)
plot(ML_Tree_Bup)
plot(ML_Tree_Hyd)
plot(ML_Tree_Hal)
plot(ML_Tree_Can)
plot(ML_Tree_Gyr)
plot(ML_Tree_Ela)
plot(ML_Tree_Cryp)
plot(ML_Tree_Sci)
plot(ML_Tree_Lat)


#Remove unneeded variables 
rm(env1, env2, env3, env4, env5, env6, env7, env8, env9, env10, env11, env12, env13, env14, env15, fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, fit11, fit12, fit13, fit14, fit15, mt1, mt2, mt3, mt4, mt5, mt6, mt7, mt8, mt9, mt10, mt11, mt12, mt13, mt14, mt15, dm1, dm2, dm3, dm4, dm5, dm6, dm7, dm8, dm9, dm10, dm11, dm12, dm13, dm14, dm15, dm16, tree1, tree2, tree3, tree4, tree5, tree6, tree7, tree8, tree9, tree10, tree11, tree12, tree13, tree14, tree15, binNames, familyFileNames, familyList, familySequenceNames, referencefind1, referencefind2, referencefind3, referencefind4, referencefind5, referencefind6, referencefind7, referencefind8, referencefind9, referencefind10, referencefind11, referencefind12, referencefind13, referencefind14, referencefind15, referencefind16)


############################
#Cam's refactoring of above change to dataframes

#note I haven't tested all this block of code yet, need Sam's original files, 
#we can go through it together and look for any mistakes I've made

#going to load all of the files into a list of objects
#this avoids having to repeat the code for each individually and is a lot safer
#because we reduce the likelihood of copy/paste errors in the code construction
list_of_files = c("AlignmentDytiscidae.fas", "AlignmentCarabidae.fas",
	"AlignmentCurculionidae.fas","AlignmentCoccinellidae.fas","AlignmentLeiodidae.fas", 
	"AlignmentChrysomelidae.fas","AlignmentBuprestidae.fas","AlignmentHydrophilidae.fas", 
	"AlignmentHaliplidae.fas", "AlignmentCantharidae.fas", "AlignmentGyrinidae.fas",
	"AlignmentElateridae.fas","AlignmentCryptophagidae.fas", "AlignmentScirtidae.fas", 
	"AlignmentLatridiidae.fas")


#this is a motif for the application of a function to a list of things
# we can talk about this in person/draw it out but its really handy once you learn the pattern

#note if any of the intermediate steps here are not needed, this could be further reduced by
#writing a function that conducts two or more of the functions below in one go
#i.e. what you've done above with RefSeqTrim
phylo_dat = lapply(list_of_files, function(x){
	read.phyDat(x, format="fasta", type="DNA")
	})

dm_outs = lapply(phylo_dat, function(x){
	dist.ml(x)
	})

tree_outs = lapply(dm_outs, function(x){
	fastme.bal(x)
	})

model_tests = lapply(phylo_dat, function(x){
	modelTest(x)
	})

env_out = lapply(model_tests, function(x){
	attr(x)
	})

#look at the get thing here and check this is okay
model_fit = lapply(env_out, function(x){
	eval(get("HKY+G+I", x),x)
	})

# two argument function so use mapply and a mapper function
# this is done so we can hardcode the k and inv paramaters
# if there were just the two variables, could use only mapply
# note the variable names aren't the best here
pml_wrapper = function(tree, phylo, inv){
	return(pml(tree, phylo, k=4, inv))
}

ml_out = mapply(pml_wrapper, tree_outs, phylo_dat, inv_are_hiding)

ml_out = lapply(ml_out, function(x){
	optim.pml(x, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
	})

ml_trees = = lapply(ml_out, function(x){
	x['tree',]
	})

#here is a single use example of the original that I was using to write the functions above
phyDatDyt <- read.phyDat("AlignmentDytiscidae.fas", format="fasta", type="DNA")
dm1 <- dist.ml(phyDatDyt)
tree1 <- fastme.bal(dm1)
mt1 <- modelTest(phyDatDyt)
env1 <-attr(mt1, "env")
fit1 <- eval(get("HKY+G+I", env1), env1)
ML_Dyt <- pml(tree1, phyDatDyt, k=4, inv= 0.4681883)
ML_Dyt <- optim.pml(ML_Dyt, optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = "HKY")
ML_Tree_Dyt <- ML_Dyt$tree
plot(ML_Tree_Dyt)


#you can take the patterns we've talked about above and apply them to Part 5 below


###########################


#Part 5: NTI and NRI----

#Create a filter for BINs found in Churchill
ChurchillFilter <- which(dfAllSeq$bin_uri %in% dfOrder_Churchill$bin_uri)
#Create a filter for the BINs not found in Churchill 
NotChurchillFilter <- which(!(dfAllSeq$bin_uri %in% dfOrder_Churchill$bin_uri))
#Apply the filters 
dfFilter_Churchill <- dfAllSeq[ChurchillFilter, ]
dfFilter_NotChurchill <- dfAllSeq[NotChurchillFilter, ]
#Change to data table and set to 1 if present in Churchill and 0 if not in Churchill
dfFilter_Churchill <- as.data.table(dfFilter_Churchill)
dfFilter_Churchill <- dfFilter_Churchill[, churchill := 1]
dfFilter_NotChurchill <- as.data.table(dfFilter_NotChurchill)
dfFilter_NotChurchill <- dfFilter_NotChurchill[, churchill := 0]
#Bind the new data frames to taxalistcomplete 
dfAllSeq <- rbind(dfFilter_Churchill, dfFilter_NotChurchill)

#Create a presence absence matrix for bin_uri in Churchill
#Create new data frame  
dfAllSeq2 <- dfAllSeq [, c("bin_uri", "churchill", "family_name")]
#Split into family dataframes and remove the family column 
dfAllSeq2 <- split(dfAllSeq2, list(dfAllSeq$family_name))
dfBuprestidae <- as.data.frame(dfAllSeq2[1])
dfBuprestidae <- dfBuprestidae[, c("Buprestidae.bin_uri", "Buprestidae.churchill")]
dfCantharidae <- as.data.frame(dfAllSeq2[2])
dfCantharidae <- dfCantharidae[, c("Cantharidae.bin_uri", "Cantharidae.churchill")]
dfCarabidae <- as.data.frame(dfAllSeq2[3])
dfCarabidae <- dfCarabidae[, c("Carabidae.bin_uri", "Carabidae.churchill")]
dfChrysomelidae <- as.data.frame(dfAllSeq2[4])
dfChrysomelidae <- dfChrysomelidae[, c("Chrysomelidae.bin_uri", "Chrysomelidae.churchill")]
dfCoccinellidae <- as.data.frame(dfAllSeq2[5])
dfCoccinellidae <- dfCoccinellidae[, c("Coccinellidae.bin_uri", "Coccinellidae.churchill")]
dfCryptophagidae <- as.data.frame(dfAllSeq2[6])
dfCryptophagidae <- dfCryptophagidae[, c("Cryptophagidae.bin_uri", "Cryptophagidae.churchill")]
dfCurculionidae <- as.data.frame(dfAllSeq2[7])
dfCurculionidae <- dfCurculionidae[, c("Curculionidae.bin_uri", "Curculionidae.churchill")]
dfDytiscidae <- as.data.frame(dfAllSeq2[8])
dfDytiscidae <- dfDytiscidae[, c("Dytiscidae.bin_uri", "Dytiscidae.churchill")]
dfElateridae <- as.data.frame(dfAllSeq2[9])
dfElateridae <- dfElateridae[, c("Elateridae.bin_uri", "Elateridae.churchill")]
dfGyrinidae <- as.data.frame(dfAllSeq2[10])
dfGyrinidae <- dfGyrinidae[, c("Gyrinidae.bin_uri", "Gyrinidae.churchill")]
dfHaliplidae <- as.data.frame(dfAllSeq2[11])
dfHaliplidae <- dfHaliplidae[, c("Haliplidae.bin_uri", "Haliplidae.churchill")]
dfHydrophilidae <- as.data.frame(dfAllSeq2[12])
dfHydrophilidae <- dfHydrophilidae[, c("Hydrophilidae.bin_uri", "Hydrophilidae.churchill")]
dfLatridiidae <- as.data.frame(dfAllSeq2[13])
dfLatridiidae <- dfLatridiidae[, c("Latridiidae.bin_uri", "Latridiidae.churchill")]
dfLeiodidae <- as.data.frame(dfAllSeq2[14])
dfLeiodidae <- dfLeiodidae[, c("Leiodidae.bin_uri", "Leiodidae.churchill")]
dfScirtidae <- as.data.frame(dfAllSeq2[15])
dfScirtidae <- dfScirtidae[, c("Scirtidae.bin_uri", "Scirtidae.churchill")]

#Create family matrices
Family_matrix1 <- melt(dfBuprestidae, id.var="Buprestidae.churchill")
Family_matrix2 <- melt(dfCantharidae, id.var="Cantharidae.churchill")
Family_matrix3 <- melt(dfCarabidae, id.var="Carabidae.churchill")
Family_matrix4 <- melt(dfChrysomelidae, id.var="Chrysomelidae.churchill")
Family_matrix5 <- melt(dfCoccinellidae, id.var="Coccinellidae.churchill")
Family_matrix6 <- melt(dfCryptophagidae, id.var="Cryptophagidae.churchill")
Family_matrix7 <- melt(dfCurculionidae, id.var="Curculionidae.churchill")
Family_matrix8 <- melt(dfDytiscidae, id.var="Dytiscidae.churchill")
Family_matrix9 <- melt(dfElateridae, id.var="Elateridae.churchill")
Family_matrix10 <- melt(dfGyrinidae, id.var="Gyrinidae.churchill")
Family_matrix11 <- melt(dfHaliplidae, id.var="Haliplidae.churchill")
Family_matrix12 <- melt(dfHydrophilidae, id.var="Hydrophilidae.churchill")
Family_matrix13 <- melt(dfLatridiidae, id.var="Latridiidae.churchill")
Family_matrix14 <- melt(dfLeiodidae, id.var="Leiodidae.churchill")
Family_matrix15 <- melt(dfScirtidae, id.var="Scirtidae.churchill")

Family_matrix1 <- with(Family_matrix1, table(Buprestidae.churchill, value))
Family_matrix2 <- with(Family_matrix2, table(Cantharidae.churchill, value))
Family_matrix3 <- with(Family_matrix3, table(Carabidae.churchill, value))
Family_matrix4 <- with(Family_matrix4, table(Chrysomelidae.churchill, value))
Family_matrix5 <- with(Family_matrix5, table(Coccinellidae.churchill, value))
Family_matrix6 <- with(Family_matrix6, table(Cryptophagidae.churchill, value))
Family_matrix7 <- with(Family_matrix7, table(Curculionidae.churchill, value))
Family_matrix8 <- with(Family_matrix8, table(Dytiscidae.churchill, value))
Family_matrix9 <- with(Family_matrix9, table(Elateridae.churchill, value))
Family_matrix10 <- with(Family_matrix10, table(Gyrinidae.churchill, value))
Family_matrix11 <- with(Family_matrix11, table(Haliplidae.churchill, value))
Family_matrix12 <- with(Family_matrix12, table(Hydrophilidae.churchill, value))
Family_matrix13 <- with(Family_matrix13, table(Latridiidae.churchill, value))
Family_matrix14 <- with(Family_matrix14, table(Leiodidae.churchill, value))
Family_matrix15 <- with(Family_matrix15, table(Scirtidae.churchill, value))

Family_matrix1[Family_matrix1 > 1] <- 1
Family_matrix2[Family_matrix2 > 1] <- 1
Family_matrix3[Family_matrix3 > 1] <- 1
Family_matrix4[Family_matrix4 > 1] <- 1
Family_matrix5[Family_matrix5 > 1] <- 1
Family_matrix6[Family_matrix6 > 1] <- 1
Family_matrix7[Family_matrix7 > 1] <- 1
Family_matrix8[Family_matrix8 > 1] <- 1
Family_matrix9[Family_matrix9 > 1] <- 1
Family_matrix10[Family_matrix10 > 1] <- 1
Family_matrix11[Family_matrix11 > 1] <- 1
Family_matrix12[Family_matrix12 > 1] <- 1
Family_matrix13[Family_matrix13 > 1] <- 1
Family_matrix14[Family_matrix14 > 1] <- 1
Family_matrix15[Family_matrix15 > 1] <- 1

Family_matrix1 <- unclass(Family_matrix1)
Family_matrix2 <- unclass(Family_matrix2)
Family_matrix3 <- unclass(Family_matrix3)
Family_matrix4 <- unclass(Family_matrix4)
Family_matrix5 <- unclass(Family_matrix5)
Family_matrix6 <- unclass(Family_matrix6)
Family_matrix7 <- unclass(Family_matrix7)
Family_matrix8 <- unclass(Family_matrix8)
Family_matrix9 <- unclass(Family_matrix9)
Family_matrix10 <- unclass(Family_matrix10)
Family_matrix11 <- unclass(Family_matrix11)
Family_matrix12 <- unclass(Family_matrix12)
Family_matrix13 <- unclass(Family_matrix13)
Family_matrix14 <- unclass(Family_matrix14)
Family_matrix15 <- unclass(Family_matrix15)

#Calculate net relatedness index (NRI) and nearest taxon index (NTI) using ML Tree
#Ensure ML tree is in correct format 
phy.dist1 <- cophenetic(ML_Tree_Bup)
phy.dist2 <- cophenetic(ML_Tree_Can)
phy.dist3 <- cophenetic(ML_Tree_Car)
phy.dist4 <- cophenetic(ML_Tree_Chry)
phy.dist5 <- cophenetic(ML_Tree_Cocc)
phy.dist6 <- cophenetic(ML_Tree_Cryp)
phy.dist7 <- cophenetic(ML_Tree_Cur)
phy.dist8 <- cophenetic(ML_Tree_Dyt)
phy.dist9 <- cophenetic(ML_Tree_Ela)
phy.dist10 <- cophenetic(ML_Tree_Gyr)
phy.dist11 <- cophenetic(ML_Tree_Hal)
phy.dist12 <- cophenetic(ML_Tree_Hyd)
phy.dist13 <- cophenetic(ML_Tree_Lat)
phy.dist14 <- cophenetic(ML_Tree_Lei)
phy.dist15 <- cophenetic(ML_Tree_Sci)

#Calculate NRI
ses.mpd.result_ML_Bup <- ses.mpd(Family_matrix1, phy.dist1, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Can <- ses.mpd(Family_matrix2, phy.dist2, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Car <- ses.mpd(Family_matrix3, phy.dist3, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Chry <- ses.mpd(Family_matrix4, phy.dist4, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Cocc <- ses.mpd(Family_matrix5, phy.dist5, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Cryp <- ses.mpd(Family_matrix6, phy.dist6, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Cur <- ses.mpd(Family_matrix7, phy.dist7, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Dyt <- ses.mpd(Family_matrix8, phy.dist8, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Ela <- ses.mpd(Family_matrix9, phy.dist9, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Gyr <- ses.mpd(Family_matrix10, phy.dist10, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Hal <- ses.mpd(Family_matrix11, phy.dist11, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Hyd <- ses.mpd(Family_matrix12, phy.dist12, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Lat <- ses.mpd(Family_matrix13, phy.dist13, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Lei <- ses.mpd(Family_matrix14, phy.dist14, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mpd.result_ML_Sci <- ses.mpd(Family_matrix15, phy.dist15, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)

#Calculate NTI
ses.mntd.result_ML_Bup <- ses.mntd(Family_matrix1, phy.dist1, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Can <- ses.mntd(Family_matrix2, phy.dist2, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Car <- ses.mntd(Family_matrix3, phy.dist3, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Chry <- ses.mntd(Family_matrix4, phy.dist4, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Cocc <- ses.mntd(Family_matrix5, phy.dist5, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Cryp <- ses.mntd(Family_matrix6, phy.dist6, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Cur <- ses.mntd(Family_matrix7, phy.dist7, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Dyt <- ses.mntd(Family_matrix8, phy.dist8, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Ela <- ses.mntd(Family_matrix9, phy.dist9, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Gyr <- ses.mntd(Family_matrix10, phy.dist10, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Hal <- ses.mntd(Family_matrix11, phy.dist11, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Hyd <- ses.mntd(Family_matrix12, phy.dist12, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Lat <- ses.mntd(Family_matrix13, phy.dist13, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Lei <- ses.mntd(Family_matrix14, phy.dist14, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
ses.mntd.result_ML_Sci <- ses.mntd(Family_matrix15, phy.dist15, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)

#Remove unneeded variables 
rm(env, Family_phyDat, fit, mt, phy.dist, dm, dfFilter_Churchill, dfFilter_NotChurchill, ChurchillFilter, NotChurchillFilter, dfAllseq2)

#Part 6: Trait Analysis: ANOVA----

#Read in character matrix
Coleoptera_Matrix_NRI <- read_csv(file="C:/Users/sammi/Dropbox/Sam Majoros/R Code/Coleoptera_Matrix_NRI.csv")
#Run ANOVAs for both traits 
Coleoptera_ANOVA_NRI_Feeding <- aov(structure ~ adult_diet, data = Coleoptera_Matrix_NRI)
Coleoptera_ANOVA_NRI_Habitat <- aov(structure ~ habitat, data = Coleoptera_Matrix_NRI)
#Get ANOVA summary 
summary(Coleoptera_ANOVA_NRI_Habitat)
summary(Coleoptera_ANOVA_NRI_Feeding)

#Repeat for NTI values
#Read in matrix
Coleoptera_Matrix_NTI <- read.csv(file="C:/Users/sammi/Dropbox/Sam Majoros/R Code/Coleoptera_Matrix_NTI.csv")
#Run ANOVAs for both traits 
Coleoptera_ANOVA_NTI_Habitat <- aov(structure ~ habitat, data = Coleoptera_Matrix_NTI)
Coleoptera_ANOVA_NTI_Feeding <- aov(structure ~ adult_diet, data = Coleoptera_Matrix_NTI)
#Get ANOVA summary 
summary(Coleoptera_ANOVA_NTI_Habitat)
summary(Coleoptera_ANOVA_NTI_Feeding)

#Part 7: Trait Analysis: PGLS----

#Read in matrix
PGLSdata_NRI <- read.csv("Coleoptera_Matrix_NRI.csv")
#Read in tree
PGLStree <- read.nexus("PGLS_tree_Coleoptera")
#Set branch lengths to one
PGLStree$edge.length <- replicate((length(PGLStree$edge[, 1])), 1)
PGGLStree <- force.ultrametric(PGLStree, method="extend")
#Set the row names to family names 
PGLSdata_NRI <- PGLSdata_NRI %>%
  column_to_rownames(var = 'family_name')
#Make sure tree and dataframe are in the same order
PGLSdata_NRI <- PGLSdata_NRI[match(PGLStree$tip.label, rownames(PGLSdata_NRI)), ]
#Run PGLS analysis 
pglsModel_NRI1 <- gls(structure ~ habitat, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NRI, method = "ML")
pglsModel_NRI2 <- gls(structure ~ adult_diet, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NRI, method = "ML")
#Get PGLS summary 
summary(pglsModel_NRI1)
summary(pglsModel_NRI2)

#Create boxplots for traits vs. clustering matrix
plot1 <- boxplot(PGLSdata_NRI$structure ~ PGLSdata_NRI$habitat )
plot2 <- boxplot(PGLSdata_NRI$structure ~ PGLSdata_NRI$adult_diet )

#Repeat for NTI
PGLSdata_NTI <- read.csv("Coleoptera_Matrix_NTI.csv")
#Set row names to family names 
PGLSdata_NTI <- PGLSdata_NTI %>%
  column_to_rownames(var = 'family_name')
#Make sure tree and dataframe are in the same order
PGLSdata_NTI <- PGLSdata_NTI[match(PGLStree$tip.label, rownames(PGLSdata_NTI)), ]
#Run PGLS analysis 
pglsModel_NTI1 <- gls(structure ~ habitat, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NTI, method = "ML")
pglsModel_NTI2 <- gls(structure ~ adult_diet, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NTI, method = "ML")
#Get PGLS summary 
summary(pglsModel_NTI1)
summary(pglsModel_NTI2)




