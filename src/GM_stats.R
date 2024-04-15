#load required Packages
library(ggplot2)
library(geomorph)
library(SlicerMorphR)
library(tidyverse)
library(gtools)


#Set a reproducibility seed
set.seed(2145)

# Slicermorph import and cross validation ---------------------------------


#this section brings in the slicermorph landmarks and configures that analysis for use in R
##Unless noted this code is the same as that provided in the slicermorph tutorials here:https://github.com/SlicerMorph/Tutorials/blob/main/GPA_3/parser_and_sample_R_analysis.md
##This code is replicated so that all the same basic sanity checks are run to ensure consistent results.
nparser <- function(logfile=NULL, forceLPS=FALSE) {
  cut = function(x) return(strsplit(x, "=")[[1]][2])
  temp = unlist(lapply(X = readLines(logfile,w=F), FUN = cut))
  log = list()
  log$input.path = "./data/LMs"
  log$output.path = "./data"
  log$files = mixedsort(list.files(log$input.path))
  log$format = temp[6]
  log$no.LM = as.numeric(temp[7])
  if (is.na(temp[8])) {
    log$skipped = FALSE
  } else {
    log$skipped = TRUE
    log$skipped.LM = unlist(strsplit(temp[8], ","))
  }
  if (is.na(temp[15])) {
    log$semi = FALSE
  } else {
    log$semi = TRUE
    log$semiLMs = unlist(strsplit(temp[15], ","))
  }
  if (as.logical(temp[9])) {
    log$scale = TRUE
  } else { 
    log$scale = FALSE
  }
  log$MeanShape = temp[10]
  log$eigenvalues = temp[11]
  log$eigenvectors = temp[12]
  log$OutputData = temp[13]
  log$pcScores = temp[14]
  log$ID = gsub(log$format, "", fixed = T, log$files)
  log$ID = gsub(".json", "", fixed = T, log$ID)
  
  if (!log$skipped) {
    log$LM = array(dim = c(log$no.LM, 3, length(log$files)), 
                   dimnames = list(1:log$no.LM, c("x", "y", "z"), log$ID))
    if (log$format == ".fcsv") {
      for (i in 1:length(log$files)) log$LM[, , i] = read.markups.fcsv(paste(log$input.path, 
                                                                             log$files[i], sep = "/"), forceLPS)
    } else { for (i in 1:length(log$files)) log$LM[, , i] = read.markups.json(paste(log$input.path, 
                                                                                    log$files[i], sep = "/"))
    }
  } else {
    log$LM = array(dim = c(log$no.LM - length(log$skipped.LM), 
                           3, length(log$files)), dimnames = list((1:log$no.LM)[-as.numeric(log$skipped.LM)], 
                                                                  c("x", "y", "z"), log$ID))
    if (log$format == ".fcsv") {
      for (i in 1:length(log$files)) log$LM[, , i] = read.markups.fcsv(paste(log$input.path, 
                                                                             log$files[i], sep = "/"), forceLPS)[-c(as.numeric(log$skipped.LM)), ]
    } else { for (i in 1:length(log$files)) log$LM[, , i] = read.markups.json(paste(log$input.path, 
                                                                                    log$files[i], sep = "/"))[-c(as.numeric(log$skipped.LM)), ]
    }
  }
  return(log)
}


#This code parses he slicermorph GPA output into R readable data
SM.log.file <- "./data/analysis.log"
SM.log <- nparser(SM.log.file, forceLPS = TRUE)

SM.output <- read.csv(file = paste(SM.log$output.path, 
                                   SM.log$OutputData, 
                                   sep = "/"))

##The next few lines of code are new additions that bring in the classification data
## and reorganizes the slicermorph data to align with the order that data appears in

###Read in the "Classifier" sheet within the Data folder
classifier <- as.data.frame(read.csv("./data/Classifier_sheet.csv", header = T))
classifier$Inst_code <- paste(classifier$Institution,classifier$Specimen,sep = "_")
clsf_reord <- classifier[,c("Inst_code","Nose","Bitework","UKC.breeding.Standard","AKC.breeding.standard")]
clsf_reord <- clsf_reord[order(clsf_reord$Inst_code),]

SM.reorder <- SM.output[order(SM.output$Sample_name),]
##After running these steps the columns SM.reorder and clsf_reorder should list specimens in the same order

#Slicermorph Tutorial Code resumes

SlicerMorph.PCs <- read.table(file = paste(SM.log$output.path, SM.log$pcScores, sep="/"), sep = ",", header = TRUE, row.names = 1)
SlicerMorph.repc<-SlicerMorph.PCs[order(rownames(SlicerMorph.PCs)),]

## Adding to match breed names with slicer data
breed_list <- read.csv("./data/Breed_list.csv")

SlicerMorph.repc$Breed <- rep(NA, nrow(SlicerMorph.repc))

for (j in 1:nrow(SlicerMorph.repc)){
  specimen_name <- unlist(str_split(rownames(SlicerMorph.repc)[j], "_"))
  which_row <- which(breed_list$Institution == specimen_name[1] & 
                       breed_list$Specimen == specimen_name[2])
  if(length(which_row) == 0) next
  SlicerMorph.repc$Breed[j] <- breed_list$Breed[which_row]
}

### NOTE FOR LINDSAY: You should only need to run the code from the tutorial this far to get the data genrated for the figures.
### You can skip th rest of this section


PD <- SM.output[,2]
if (!SM.log$skipped) {
  no.LM <- SM.log$no.LM 
} else {
  no.LM <- SM.log$no.LM - length(SM.log$skipped.LM)
}
PD
no.LM

Coords <- arrayspecs(SM.output[, -c(1:3)], p = no.LM, k = 3)

dimnames(Coords) <- list(1:no.LM, c("x","y","z"),SM.log$ID)

gdf <- geomorph.data.frame(Size = SM.output$centeroid, Coords = Coords)
fit.slicermorph <- procD.lm(Coords ~ Size, data = gdf, print.progress = FALSE)

gpa <- gpagen(SM.log$LM)
pca <- gm.prcomp(gpa$coords)
geomorph.PCs <- pca$x

gdf2 <- geomorph.data.frame(size = gpa$Csize, coords = gpa$coords)
fit.rawcoords <- procD.lm(coords ~ size, data = gdf2, print.progress = FALSE)

pd <- function(M, A) sqrt(sum(rowSums((M-A)^2)))
geomorph.PD <- NULL
for (i in 1:length(SM.log$files)) 
  geomorph.PD[i] <- pd(gpa$consensus, gpa$coords[,, i])

plot(gpa$Csize, SM.output$centeroid, 
     pch = 20, ylab = 'SlicerMorph', 
     xlab = 'geomorph', main = "Centroid Size")
cor(gpa$Csize, SM.output$centeroid)

plot(geomorph.PD, SM.output[,2], 
     pch = 20, ylab = 'SlicerMorph', 
     xlab = 'geomorph', main = "Procrustes Distance")
cor(geomorph.PD, SM.output[,2])

plot(geomorph.PCs[,1], SlicerMorph.PCs[,1], 
    pch = 20, ylab = 'SlicerMorph', 
    xlab = 'geomorph', main = "PC1 scores")

plot(geomorph.PCs[,2], SlicerMorph.PCs[,2], 
     pch = 20, ylab = 'SlicerMorph', 
     xlab = 'geomorph', main = "PC2 scores")

for (i in 1:10) print(cor(SlicerMorph.PCs[, i], 
                          geomorph.PCs[, i]))

summary(fit.slicermorph)
summary(fit.rawcoords)
#End of Slicermorph tutorial code and sanity check


# KMeans analysis ---------------------------------------------------------

#First make a function to run the K-means for a specified number of clusters on the first two PC axes
Kanalyses<-function(clusters) {
  kmeans(SlicerMorph.repc[,1:2],clusters)
}

#Next inspect different k values to find the optimal number of clusters
#First make a storage variable
km.var<-c()

#next a for loop that pulls the total sum of the squares for each k means cluster size 
for (k in 1:15) {
  holding<-Kanalyses(k)
  km.var[k]<-holding$tot.withinss
}

#now a simple plot to inspect those results
k.clusts<-data.frame(x=1:15,y=km.var)
ggplot(k.clusts,aes(x=x,y=y))+geom_line(color="#23669D")+geom_point()+labs(x="number of clusters",y="Total Within Sum of Squares ")

#now view the kmeans for whichever number of clusters is most appropriate
km<-Kanalyses(3)
