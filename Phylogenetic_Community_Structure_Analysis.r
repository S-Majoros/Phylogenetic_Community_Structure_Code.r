#Pipeline for phylogenetic community structure analysis
#Centroid and alignment code adapted from Matthew Orton https://github.com/m-orton/Evolutionary-Rates-Analysis-Pipeline/blob/master/EvolutionaryComparisonPipelineSmallTaxa.R
#Sequence trimming and outlier removal functions adapted from Jacqueline May https://github.com/jmay29/phylo/blob/master/refSeqTrim.R
#Script by Sam Majoros. Email: smajoros@uoguelph.ca
#Thank you to Cameron Nugent, Jacqueline May, and Matthew Orton for help with the code design.
#Thank you to Sarah Adamowicz for help designing the project and coauthoring the paper. 
#Thank you to Alex Smith and Kamil Chatila-Amos for helpful project disscussions. 

#Part 1: Inputting and Filtering Data ----

#This sections includes steps to load the needed packages, load the data into R, and to filter the data.
#Packages needed for this analysis
#If you do not already have these packages, uncomment the code and install.
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
#install.packages("BiocManager")
#BiocManager::install(c("Biostrings", "muscle"))
library(Biostrings)
library(muscle)
#install.packages("phangorn")
library(phangorn)
#install.packages("picante")
library(picante)
#install.packages("data.table")
library(data.table)
#install.packages("phytools")
library(phytools)

#Upload order data into R
#Uncomment the following code to download data directly from BOLD, specifying the required geographical locations. The data I used was retrieved on June 19th, 2019. 
#dfOrder <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Coleoptera&geo=Alaska|Canada&format=tsv")
#Write file to hard disk
#write_tsv(dfOrder, "Coleoptera_download_June19")
#Read in saved order data
dfOrder <- read_tsv("Coleoptera_download_June19")

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
#Remove redundant "BOLD:" component from each BIN uri in the BIN column
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
#First convert to datatable
dfOrder_Churchill <- as.data.table(dfOrder_Churchill)
#create datatable showing number of sequences per family
total_species_number <- dfOrder_Churchill[ , .(.N),by=.(family_name)]
#create datatable showing the number of BINs per family
number_of_unique_species <- dfOrder_Churchill[ , .(number_of_species=length(unique(bin_uri))), by=family_name]
#convert to dataframe
number_of_unique_species <- as.data.frame(number_of_unique_species)
#Filter down to families with 3 or more species
number_of_unique_species <- filter(number_of_unique_species, number_of_unique_species$number_of_species > 2)

#Create filter to filter down the order to families with three or more species in the subset
dfOrder_filter <- which(dfOrder$family_name %in% number_of_unique_species$family_name)
#Apply filter
dfOrder <- dfOrder[dfOrder_filter, ]
                    
#Here I save the filtered data as a file. This data is provided on Github as the original download is too large. 
#To read in the pre filtered data, uncomment the following code. 
#dfOrder <- read_tsv("Coleoptera_download_post_filtering")                    

#Remove unneeded variables
rm(number_of_unique_species, total_species_number, SubsetFilter_Churchill, containLatLon, endNGap, internalNGap, nGapCheck, split, startNGap, dfOrder_filter, i)

#Part 2: Choosing a Centroid----

#In this section we find a centroid sequence for each BIN present in the Canada and Alaska dataset. This includes running a multiple sequence alignment and computing a distance matrix. 
#Create smaller dataframe with needed info
dfBinList <- (dfOrder[, c("processid", "bin_uri", "nucleotides")])
#Create groupings by BIN, each with different bin_uri
binList <- lapply(unique(dfOrder$bin_uri), function(x) dfOrder[dfOrder$bin_uri==x, ])
#Find the number of processids in each bin
binSize <- sapply(binList, function(x)length(x$processid))
#Create new data frame with bin_uri and bin size
dfOrder_bins <- data.frame(binSize)
dfOrder_bins$bin_uri <- c(unique(dfOrder$bin_uri))
#Merge dfBinList and dfOrder_bins
dfBinList <- merge(dfBinList, dfOrder_bins, by.x="bin_uri", by.y="bin_uri")
#Reorder dfFamily_bins by bin_uri
dfOrder_bins <- dfOrder_bins[order(dfOrder_bins$bin_uri), ]

#Find BINs with more than one member
largeBin <- which(dfBinList$binSize > 1)
#Create dataframe with only bins with more than one member
if (length(largeBin) > 0) {
  dfCentroid <- dfBinList[largeBin, ]
}

#Subset dfOrder_bins down to number of BINs in dfOrder
dfOrder_bins <- subset(dfOrder_bins, dfOrder_bins$bin_uri %in% dfCentroid$bin_uri)

#Find number of unique BINs in dfCentroid
binNumberCentroid <- unique(dfCentroid$bin_uri)
binNumberCentroid <- length(binNumberCentroid)

#Create dataframe with BINs with only one sequence
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

#Here I saved dfAllSeq to the file. This Centroid step takes a long time to complete, so the data post centroid is provided on Github
#If you want to use the post centroid data, uncomment and run the following code.
#dfAllSeq <- read_tsv("Post_Centroid")                     
                        
#Remove unneeded dataframes and variables
rm(alignment1, binList, dfBinList, dfCentroid, dfOrder_bins, dfNonCentroid, dnaBINCentroid, dnaStringSet1, geneticDistanceCentroid, largeBinList, largeBinProcessid, binNumberCentroid, binSize,  centroidSeq, i, largeBin)

#Part 3: Alignment----

#This sections includes code to trim the sequences and run a mutliple sequence alignment. The first alignment is a preliminary alignment that allows the sequences to be trimmed. The second is the final alignment once the outgroups have been added. Outgroups for each family were chosen from a different suborder of Coleoptera.
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
  #Find stop and start positions in reference
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
#Visualize the values in the distance matrix using a histogram
hist(distanceMatrix)

#Using upper threshold of IQR to detect outliers
lowerQuantile <- quantile(distanceMatrix)[2]
upperQuantile <- quantile(distanceMatrix)[4]
iqr <-upperQuantile - lowerQuantile
#Set threshold to 1.5. In order to only remove extreme outliers this can be change to 3.
upperThreshold <- (iqr*1.5) + upperQuantile
#Remove 0 values
distanceMatrix[distanceMatrix==0] <- NA
#Convert to data table
dfOutliers <- as.data.table(distanceMatrix, keep.rownames = T)
#Change the "rn" column to bin_uri
setnames(dfOutliers, "rn", "bin_uri")
#Identify divergent BINs
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

#name nucleotide column and set as character
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

#Create new dataframes that sort families by suborder
dfAdephaga <- rbind(taxalistcomplete[[1]], taxalistcomplete[[3]], taxalistcomplete[[10]], taxalistcomplete[[12]])

dfPolyphaga <- rbind(taxalistcomplete[[2]], taxalistcomplete[[4]], taxalistcomplete[[5]], taxalistcomplete[[6]], taxalistcomplete[[7]],taxalistcomplete[[8]],taxalistcomplete[[9]],taxalistcomplete[[11]],taxalistcomplete[[13]],taxalistcomplete[[14]],taxalistcomplete[[15]],taxalistcomplete[[16]])

#Randomly select outgroup sequences
dfAdephaga_Outgroup <- sample_n(dfAdephaga, 1)
dfPolyphaga_outgroup <- sample_n(dfPolyphaga, 1)    
                           
#Put outgroup dataframes into a list
Suborders <- list(dfPolyphaga_outgroup, dfAdephaga_Outgroup)                  

dataframeNums <- c(1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2)
#Creating a vector of the outgroup names.
outgroupNames <- c("Outgroup_Car", "Outgroup_Cur", "Outgroup_Dyt", "Outgroup_Cocc", "Outgroup_Lei", "Outgroup_Chry", "Outgroup_Staph", "Outgroup_Bup", "Outgroup_Hydro", "Outgroup_Hal", "Outgroup_Can", "Outgroup_Gyr", "Outgroup_Ela", "Outgroup_Crypt", "Outgroup_Scirt", "Outgroup_Lat")

 l_outgroups <- list()

#Use a for loop to sample the sequences from taxalistcomplete.
for (i in 1:length(dataframeNums)) {
#Take the ith dataframe number.
  number <- dataframeNums[[i]]
#Sample from the corresponding dataframe in taxalistcomplete.
 outgroupSeqs <- Suborders[[number]]
#Append outgroupSeqs to the list we created earlier.
l_outgroups[[i]] <- outgroupSeqs

}
names(l_outgroups) <- outgroupNames

taxalistcompleteWO <- list()
#Now let's add the outgroup sequences to the corresponding dataframes taxalistcomplete.
for (i in 1:length(taxalistcomplete)) {
  
  #Take the ith outgroupseqs.
  outgroupSeqs <- l_outgroups[[i]]
  #Take the ith taxa from taxalistcomplete.
  taxa <- taxalistcomplete[[i]]
  #Combine them using rbind
  dfCombined <- rbind(taxa, outgroupSeqs)
  #Append to our list.
  taxalistcompleteWO[[i]] <- dfCombined
  
}                          
      
#Trim sequences again now that outgroups have been added
taxalistcompleteWO <- lapply(taxalistcompleteWO, RefSeqTrim)                         

#Extract sequences and bin_uri
familyBin <- foreach(i=1:length(taxalistcompleteWO)) %do% taxalistcompleteWO[[i]]$bin_uri
familySequences <- foreach(i=1:length(taxalistcompleteWO)) %do% taxalistcompleteWO[[i]]$nucleotides
familySequenceNames <- familyBin

#Take reference sequences
alignmentref <- as.character(dfRefSeq$nucleotides)
dfRefSeq$reference <- "reference"
#Name reference as a reference
alignmentRefNames <- dfRefSeq$reference
#Merge reference with other sequences
alignmentSequencesPlusRef <- foreach(i=1:length(taxalistcompleteWO)) %do%
  append(familySequences[[i]], alignmentref[[1]])

#Merge names together
alignmentNames <- foreach(i=1:length(taxalistcompleteWO)) %do%
  append(familySequenceNames[[i]], alignmentRefNames[[1]])

#Convert sequences to DNAStringSet format
dnaStringSet3 <- foreach(i=1:length(alignmentSequencesPlusRef)) %do%
  DNAStringSet(alignmentSequencesPlusRef[[i]])

  #Name each sequence
  for(i in 1:16){
    names(dnaStringSet3[[i]]) = alignmentNames[[i]]
  }

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

#Part 4: Create Maximum Likelihood tree----

#In this section, the data is saved to a fasta file for each family, distance matrices and neighbour joining trees are found, the best model and paramters are found for each family and maximum likelihood trees are generated. A bootstrap analysis is perfromed for each family and maxCladeCred consensus trees were found.
#Create function to convert DNAStringSets to dataframes
dna_string_to_df = function(dna_string_set){
  out_df = as.data.frame(dna_string_set[[1]])
  for(i in 2:length(dna_string_set)){
    new_df = as.data.frame(dna_string_set[[i]])
    out_df = rbind(out_df, new_df)
  }
  return(out_df)
}
#convert stringsets to dataframes
FamilyDNA = dna_string_to_df(dnaStringSet4)
                           
#Add the bin_uri
FamilyDNA$bin_uri <- row.names(FamilyDNA)
#Bind list of dataframes into one dataframe
dfTaxaWO <- bind_rows(taxalistcompleteWO, .id = "column_label")
# Add column_label info to dfAllSeq.
dfAllSeqWO <- merge(dfAllSeq, dfTaxaWO[, c(1:2)], by = "bin_uri", all.y = T)
# Merge with the information for dfAllSeqWO.
dfFamilyDNA <- merge(FamilyDNA, dfAllSeqWO, by = "bin_uri",  all.x = TRUE)
#Rename the coloumn with your aligned sequences
colnames(dfFamilyDNA)[2] <- "FinalSequences"

#create function to get reference names
get_reference_names <- function(top_ref_num = 15){
  ref_names <- c("reference")
  prefix <- "reference"
  for(i in 1:top_ref_num){
    new_str <- paste(prefix, as.character(i), sep='')
    ref_names <- c(ref_names, new_str)
  }
  return(ref_names)
}
#create list of reference names
reference_names = get_reference_names()
#remove reference from dataframe
dfFamilyDNA <- dfFamilyDNA[!dfFamilyDNA$bin_uri %in% reference_names , ]

#Groups by column label
familyList <- lapply(na.omit(unique(dfFamilyDNA$column_label)), 
                     function(x) dfFamilyDNA[dfFamilyDNA$column_label == x, ])                           

 #Create function to remove NAs                        
 RemoveAllNAs <- function(dfOriginal) {
  
 #Function for removing rows that consist of entirely NA values.

 #Reset the row numbers.
  row.names(dfOriginal) <- NULL
  # First check if there are any rows with all values missing.
  test <- apply(dfOriginal, 1, function(x) all(is.na(x)))
  
  if (any(test) == TRUE) {
    
    #Identify the missing rows.
    missing <- which(test == TRUE)
    #Remove them from dfOriginal.
    dfOriginal <- dfOriginal[-missing, ]
    
  } 
  
  return(dfOriginal)
  
}

#Cleaning up familyList.
familyList <- lapply(familyList, RemoveAllNAs)
                
familyList <- lapply(familyList, function(x) {
  
  noBin <- is.na(x$bin_uri)
  x <- x[!noBin, ]
  
})               
                           
#Create new dnaStringSet
dnaStringSet5 <- sapply(familyList, function(x) DNAStringSet(x$FinalSequences))
#Pull BIN names from list
binNames <- sapply(familyList, function(x)(x$bin_uri))
#Name the stringsets
for(i in seq(from = 1, to = length(dnaStringSet5), by = 1)) {
  names(dnaStringSet5[[i]]) <- binNames[[i]]
}

#Save family as a fasta file
#For file names make sure to list each family name
familyFileNames <- list("Carabidae", "Curculionidae", "Dytiscidae", "Coccinellidae", "Leiodidae", "Chrysomelidae", "Staphylinidae", "Buprestidae", "Hydrophilidae", "Haliplidae", "Cantharidae", "Gyrinidae", "Elateridae", "Cryptophagidae", "Scirtidae", "Latridiidae")
#Add alignment and .fas to each family name
familyFileNames <- foreach(i=1:length(familyFileNames)) %do%
  paste("Alignment", familyFileNames[[i]], ".fas", sep="")
#Send to your desired working directory
foreach(i=1:length(dnaStringSet5)) %do% writeXStringSet(dnaStringSet5[[i]], file=familyFileNames[[i]], format="fasta")

#create a list of alignment files
#Calling the alignments in alphabetical order allows for easier analysis during the NRI/NTI step
list_of_files <- c("AlignmentBuprestidae.fas", "AlignmentCantharidae.fas", "AlignmentCarabidae.fas",
                   "AlignmentChrysomelidae.fas", "AlignmentCoccinellidae.fas", "AlignmentCryptophagidae.fas",
                   "AlignmentCurculionidae.fas", "AlignmentDytiscidae.fas", "AlignmentElateridae.fas",
                   "AlignmentGyrinidae.fas","AlignmentHaliplidae.fas", "AlignmentHydrophilidae.fas",
                   "AlignmentLatridiidae.fas", "AlignmentLeiodidae.fas", "AlignmentScirtidae.fas",
                   "AlignmentStaphylinidae.fas")

#read the alignments into phyDat format
phylo_dat <- lapply(list_of_files, function(x){
  read.phyDat(x, format="fasta", type="DNA")
})

#create distance matrices
dm <- lapply(phylo_dat, function(x){
  dist.ml(x)
})

#creating NJ tree
tree <- lapply(dm, function(x){
  NJ(x)
})

#run model tests
model_tests <- lapply(phylo_dat, function(x){
  modelTest(x)
})

#create environments
env <- lapply(model_tests, function(x){
  attr(x, "env")
})

#create function to find best model for each family
get_best_model = function(model_df){
  best_model = model_df['Model'][model_df['BIC'] == min(model_df['BIC']) ]
  return(best_model)
}

#create a vector containing the best models
list_of_models = unlist(lapply(model_tests, function(x){
  get_best_model(x)
}))


#get parameters for each model
model_fit <- lapply(env, function(x){
  eval(get(list_of_models, x),x)
})

#create vector containing inv values
inv_values <- lapply(1:length(model_fit), function(i){
  model_fit[[i]]$inv
})

#compute likelihood
ml_out = lapply(1:length(tree), function(i){
  pml(tree[[i]], phylo_dat[[i]], k=4, inv = inv_values[[i]])
})

#drop the suffix from each of the model names
new_list_of_models = unlist(lapply(list_of_models , function(x){unlist(strsplit(x, "\\+"))[[1]]}))

#compute likelihood and optimize parameters
ml_families = lapply(1:length(ml_out), function(i){
  optim.pml(ml_out[[i]], optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = new_list_of_models[[i]])
})

#Do a bootstrap analysis for each family
bs <- lapply(1:length(ml_families), function(i){
  bootstrap.pml(ml_families[[i]], bs=1000, optNni=TRUE, multicore = FALSE)
})

#Step to create a maxCladeCred consensus tree 
bsTree <- lapply(1:length(bs), function(i){
  maxCladeCred(bs[[i]])
})
                   
#Reorder outgroup dataframe list so it is in the same order as the trees
l_outgroups <- l_outgroups[order(names(l_outgroups))]

#Root the trees using the outgroups
Tree_Rooted <- lapply(1:length(bsTree), function(i){
  root(bsTree[[3]], outgroup = l_outgroups[[3]]$bin_uri, resolve.root = TRUE)
})

#Remove the outgroup
Tree_Final <- lapply(1:length(Tree_Rooted), function(i){
  drop.tip(Tree_Rooted[[i]], l_outgroups[[i]]$bin_uri)
})


#Remove unneeded variables
rm(env, tree, model_tests, model_fit, dm, binNames, familyFileNames, familyList, dfFamilyDNA, dnaStringSet4, dnaStringSet5, FamilyDNA, ml_families, ml_out)

#Part 5: NTI and NRI----

#In this sections, family matrices are created and net relatedness indces(NRI) and nearest taxon indces (NTI) are calculated. NRI and NTI calculate the pairwise distance between two species and use this to estimate the community relatedness. NRI averages the evolutionary distances between all pairs of tips in the community, while NTI takes only the distances between nearest neighbours. A presence/absence matrix is created for each family in order to see how thier distribution relates to thier phylogenetic relationships.            
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
#Split into family dataframes
dfAllSeq2 <- split(dfAllSeq2, list(dfAllSeq$family_name))
#Remove the family name column
dfAllSeq2 <- lapply(dfAllSeq2, function(x){
  x[,-3]
})

#Create family matrices
#create family matrices
Family_matrices <- lapply(dfAllSeq2, function(x){
  melt(x, id.var="churchill")
})
Family_matrices <- lapply(Family_matrices, as.data.frame)
Family_matrices <- lapply(Family_matrices, function(x){
  with(x, table(churchill, value))
})
Family_matrices <- lapply(Family_matrices, unclass)

#Calculate net relatedness index (NRI) and nearest taxon index (NTI) using ML Tree
#Ensure ML tree is in correct format
phy.dist <- lapply(Tree_Final, cophenetic)

#Calculate NRI
NRI_Results = lapply(1:length(phy.dist), function(i){
  ses.mpd(Family_matrices[[i]], phy.dist[[i]], null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
})

#Calculate NTI
NTI_Results = lapply(1:length(phy.dist), function(i){
  ses.mntd(Family_matrices[[i]], phy.dist[[i]], null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
})

#Remove unneeded variables
rm(phy.dist, dfFilter_Churchill, dfFilter_NotChurchill, ChurchillFilter, NotChurchillFilter, dfAllSeq2, Family_matrices)

#Part 6: Trait Analysis: ANOVA----

#In this section, character matrices are loaded in and ANOVAs are performed. The ANOVA will tell us whether the phylogenetic community structure is related to the biological traits of the families. A P-value of > 0.5 indicates that families with that trait are significantly phylogenetically clustered or overdispersed. 
#Read in character matrix
Coleoptera_Matrix_NRI <- read_csv(file="C:/Users/Sam/Documents/Coleoptera_Matrix_NRI.csv")
#Run ANOVAs for all traits 
Coleoptera_ANOVA_NRI_Feeding_Adult <- aov(structure ~ adult_diet, data = Coleoptera_Matrix_NRI)
Coleoptera_ANOVA_NRI_Feeding_Larval <- aov(structure ~ larval_diet, data = Coleoptera_Matrix_NRI)
Coleoptera_ANOVA_NRI_Habitat <- aov(structure ~ habitat, data = Coleoptera_Matrix_NRI)
#Get ANOVA summary 
summary(Coleoptera_ANOVA_NRI_Habitat)
summary(Coleoptera_ANOVA_NRI_Feeding_Adult)
summary(Coleoptera_ANOVA_NRI_Feeding_Larval)

#Repeat for NTI values
#Read in matrix
Coleoptera_Matrix_NTI <- read.csv(file="C:/Users/Sam/Documents/Coleoptera_Matrix_NTI.csv")
#Run ANOVAs for both traits 
Coleoptera_ANOVA_NTI_Habitat <- aov(structure ~ habitat, data = Coleoptera_Matrix_NTI)
Coleoptera_ANOVA_NTI_Feeding_Adult <- aov(structure ~ adult_diet, data = Coleoptera_Matrix_NTI)
Coleoptera_ANOVA_NTI_Feeding_Larval <- aov(structure ~ larval_diet, data = Coleoptera_Matrix_NTI)
#Get ANOVA summary 
summary(Coleoptera_ANOVA_NTI_Habitat)
summary(Coleoptera_ANOVA_NTI_Feeding_Adult)
summary(Coleoptera_ANOVA_NTI_Feeding_Larval)

#Part 7: Trait Analysis: PGLS----

#In this section, the phylogenetic tree created based on the literature is read in and phylogenetic generalized least squares analyses are performed. 
#The PGLS will tell us whether the phylogenetic community structure is related to the biological traits of the families, while taking into account the phylogenetic relationships. A P-value of > 0.5 indicates that families with that trait are significantly phylogenetically clustered or overdispersed.
#Read in matrix
PGLSdata_NRI <- read.csv("Coleoptera_Matrix_NRI.csv")
#Read in tree
PGLStree <- read.nexus("Coleoptera Tree 2020")
#Set branch lengths to one
PGLStree$edge.length <- replicate((length(PGLStree$edge[, 1])), 1)
PGLStree <- force.ultrametric(PGLStree, method="extend")
#Set the row names to family names 
PGLSdata_NRI <- PGLSdata_NRI %>%
  column_to_rownames(var = 'family_name')
#Make sure tree and dataframe are in the same order
PGLSdata_NRI <- PGLSdata_NRI[match(PGLStree$tip.label, rownames(PGLSdata_NRI)), ]
#Run PGLS analysis 
pglsModel_NRI_Habitat <- gls(structure ~ habitat, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NRI, method = "ML")
pglsModel_NRI_Adult_Diet <- gls(structure ~ adult_diet, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NRI, method = "ML")
pglsModel_NRI_Larval_Diet <- gls(structure ~ larval_diet, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NRI, method = "ML")
#Get PGLS summary 
summary(pglsModel_NRI_Habitat)
summary(pglsModel_NRI_Adult_Diet)
summary(pglsModel_NRI_Larval_Diet)

#Create boxplots for traits vs. clustering matrix
plot_Habitat_NRI <- boxplot(PGLSdata_NRI$structure ~ PGLSdata_NRI$habitat)
plot_Adult_Diet_NRI <- boxplot(PGLSdata_NRI$structure ~ PGLSdata_NRI$adult_diet)
plot_Larval_Diet_NRI <- boxplot(PGLSdata_NRI$structure ~ PGLSdata_NRI$larval_diet)

#Repeat for NTI
PGLSdata_NTI <- read.csv("Coleoptera_Matrix_NTI.csv")
#Set row names to family names 
PGLSdata_NTI <- PGLSdata_NTI %>%
  column_to_rownames(var = 'family_name')
#Make sure tree and dataframe are in the same order
PGLSdata_NTI <- PGLSdata_NTI[match(PGLStree$tip.label, rownames(PGLSdata_NTI)), ]
#Run PGLS analysis 
pglsModel_NTI_Habitat <- gls(structure ~ habitat, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NTI, method = "ML")
pglsModel_NTI_Adult_Diet <- gls(structure ~ adult_diet, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NTI, method = "ML")
pglsModel_NTI_Larval_Diet <- gls(structure ~ larval_diet, correlation = corBrownian(phy = PGLStree), data = PGLSdata_NTI, method = "ML")
#Get PGLS summary 
summary(pglsModel_NTI_Habitat)
summary(pglsModel_NTI_Adult_Diet)
summary(pglsModel_NTI_Larval_Diet)

#Create boxplots for traits vs. clustering matrix
plot_Habitat_NTI <- boxplot(PGLSdata_NTI$structure ~ PGLSdata_NTI$habitat)
plot_Adult_Diet_NTI <- boxplot(PGLSdata_NTI$structure ~ PGLSdata_NTI$adult_diet)
plot_Larval_Diet_NTI <- boxplot(PGLSdata_NTI$structure ~ PGLSdata_NTI$larval_diet)

