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
dfOrder <- read_tsv("Coleoptera_download_Oct26")

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

#Part 4: Create Maximum Liklihood tree----

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
#Merge with the information for dfAllSeq
#This step does not work with the new code. Issue with differing number of rows.
dfFamilyDNA <- merge(FamilyDNA, dfAllSeq, by.x = "bin_uri", by.y = "bin_uri", all.x = TRUE)
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

################################################################
#for cam's computer only:
list_of_files <- c("../data/AlignmentDytiscidae.fas", "../data/AlignmentCarabidae.fas",
                   "../data/AlignmentCurculionidae.fas","../data/AlignmentCoccinellidae.fas","../data/AlignmentLeiodidae.fas",
                   "../data/AlignmentChrysomelidae.fas","../data/AlignmentBuprestidae.fas","../data/AlignmentHydrophilidae.fas",
                   "../data/AlignmentHaliplidae.fas", "../data/AlignmentCantharidae.fas", "../data/AlignmentGyrinidae.fas",
                   "../data/AlignmentElateridae.fas","../data/AlignmentCryptophagidae.fas", "../data/AlignmentScirtidae.fas",
                   "../data/AlignmentLatridiidae.fas")
###################################################################################

#create a list of alignment files
list_of_files <- c("AlignmentBuprestidae.fas", "AlignmentCantharidae.fas", "AlignmentCarabidae.fas",
                   "AlignmentChrysomelidae.fas", "AlignmentCoccinellidae.fas", "AlignmentCryptophagidae.fas",
                   "AlignmentCurculionidae.fas", "AlignmentDytiscidae.fas", "AlignmentElateridae.fas",
                   "AlignmentGyrinidae.fas","AlignmentHaliplidae.fas", "AlignmentHydrophilidae.fas",
                   "AlignmentLatridiidae.fas", "AlignmentLeiodidae.fas", "AlignmentScirtidae.fas")

#read the alignments into phyDat format
phylo_dat <- lapply(list_of_files, function(x){
  read.phyDat(x, format="fasta", type="DNA")
})

#create distance matrices
dm <- lapply(phylo_dat, function(x){
  dist.ml(x)
})

#Create a tree for each family
tree <- lapply(dm, function(x){
  fastme.bal(x)
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

#I added a line to the function above and it removes everything after the '+' in the model name
new_list_of_models == model_list
#note these don't match!
#possible reasons: a.error in the model_list vectore
#                   b. The code to generate list_of_models is picking the minimum incorrectly, Not BIC?
#                   c. The models are conducting some permutation and the results are different each time
#                       in this case we may need to set a random seed for the script:
#                       setting seed: http://rfunction.com/archives/62

#compute likelihood and optimize parameters
ml_families = lapply(1:length(ml_out), function(i){
  optim.pml(ml_out[[i]], optNni = TRUE, optGamma = TRUE, optInv = TRUE, model = model_list[[i]])
})

#Create seperate variable for trees
ML_Trees <- lapply(ml_families, function(x){
  x$tree})

#Remove unneeded variables
rm(env, tree, model_tests, model_fit, dm, binNames, familyFileNames, familyList, dfFamilyDNA, dnaStringSet4, dnaStringSet5, FamilyDNA, ml_families, ml_out)

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
#Split into family dataframes
dfAllSeq2 <- split(dfAllSeq2, list(dfAllSeq$family_name))
#Removing Staphylinidae until issue is sorted out
dfAllSeq2[16] <- NULL
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
phy.dist <- lapply(ML_Trees, cophenetic)

#Calculate NRI
NRI_Results = lapply(1:length(phy.dist), function(i){
  ses.mpd(Family_matrices[[i]], phy.dist[[i]], null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
})

#Calculate NTI
NTI_Results = lapply(1:length(phy.dist), function(i){
  ses.mntd(Family_matrices[[i]], phy.dist[[i]], null.model = "taxa.labels", abundance.weighted = FALSE, runs = 1000)
})

#Remove unneeded variables
rm(Family_phyDat, phy.dist, dfFilter_Churchill, dfFilter_NotChurchill, ChurchillFilter, NotChurchillFilter, dfAllSeq2, Family_matrices)

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
