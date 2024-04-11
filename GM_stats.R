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
  log$input.path = "./Data/LMs"
  log$output.path = "./Data"
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
SM.log.file <- ".\\Data\\analysis.log"
SM.log <- nparser(SM.log.file, forceLPS = TRUE)

SM.output <- read.csv(file = paste(SM.log$output.path, 
                                   SM.log$OutputData, 
                                   sep = "/"))

##The next few lines of code are new additions that bring in the classification data
## and reorganizes the slicermorph data to align with the order that data appears in

###Read in the "Classifier" sheet within the Data folder
classifier <- as.data.frame(read.csv(file = file.choose(), header = T))
classifier$Inst_code <- paste(classifier$Institution,classifier$Specimen,sep = "_")
clsf_reord <- classifier[,c("Inst_code","Nose","Bitework","UKC.breeding.Standard","AKC.breeding.standard")]
clsf_reord <- clsf_reord[order(clsf_reord$Inst_code),]

SM.reorder <- SM.output[order(SM.output$Sample_name),]
##After running these steps the columns SM.reorder and clsf_reorder should list specimens in the same order

#Slicermorph Tutorial Code resumes

SlicerMorph.PCs <- read.table(file = paste(SM.log$output.path, SM.log$pcScores, sep="/"), sep = ",", header = TRUE, row.names = 1)
SlicerMorph.repc<-SlicerMorph.PCs[order(rownames(SlicerMorph.PCs)),]
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

##First generate some generic plotting functions that will be useful for the remaining plots
plotbuild<-function(datagroup,colgroup){
  ggplot(SlicerMorph.repc, aes(x=SlicerMorph.repc[,1],y=SlicerMorph.repc[,2],color = colgroup, label = rownames(SlicerMorph.repc)))
}
Pointadd<-function(datagroup){
  geom_point(SlicerMorph.repc, mapping = aes(x=PC.1,y=PC.2, shape = datagroup))
}
plotshapes<-scale_shape_manual(values=c(5,1,19,15,17,19,15,17,19,15,17))
plotcolors<-scale_color_viridis_d(option = "plasma")

#now build the plot
plotbuild(factor(km$cluster),factor(km$cluster))+Pointadd(factor(km$cluster))+
  plotshapes+plotcolors+
  stat_ellipse(geom="polygon",level = .95,fill=1,alpha=.1)+
  labs(title="k-means optimal clustering",x="PC 1 (Var=50.2%)",y="PC 2 (Var=13.3%)",color="kmeans",shape="kmeans")

# grouping plots ----------------------------------------------------------
#First added your ordered classifiers onto the PCA data
SlicerMorph.repc$UKC<-clsf_reord$UKC.breeding.Standard
SlicerMorph.repc$AKC<-clsf_reord$AKC.breeding.standard
SlicerMorph.repc$Nosework<-clsf_reord$Nose
SlicerMorph.repc$Bitework<-clsf_reord$Bitework

#Now to ease viewing we'll factor them into broad categories
##UKC first
SlicerMorph.repc$UKC.shapefact<-factor(SlicerMorph.repc$UKC, levels = c("FOX","NATURAL","Herding","Terrier","Companion","Gun dog","Guardian","MIXED","Sighthound","Northern Breed","Scenthound"))
SlicerMorph.repc$UKC.fact<-factor(SlicerMorph.repc$UKC)
levels(SlicerMorph.repc$UKC.fact)<-c("companion","FOX","Working","Hunting","Herding","companion","NATURAL","Working","Hunting","Hunting","Terrier")
SlicerMorph.repc$UKC.fact<-factor(SlicerMorph.repc$UKC.fact,levels = c("FOX","NATURAL","Herding","Terrier","Working","Hunting","companion"))

##AKC next
SlicerMorph.repc$AKC.shapefact<-factor(SlicerMorph.repc$AKC, levels = c("FOX","NATURAL","Herding","Terrier","Toy","Sporting","Working","Non-sporting","Hound","MIXED"))
SlicerMorph.repc$AKC.fact<-factor(SlicerMorph.repc$AKC)
levels(SlicerMorph.repc$AKC.fact)<-c("FOX","Herding","Hunting","companion","NATURAL","companion","Hunting","Terrier","companion","Working")
SlicerMorph.repc$UKC.fact<-factor(SlicerMorph.repc$AKC.fact,levels = c("FOX","NATURAL","Herding","Terrier","Working","Hunting","companion"))

#now a whole bunch of plotting functionsto visualize the results
#First set which group you want plots of from the options: "UKC" "AKC" "Nosework" "Bitework"
grptoggle<-"Bitework"
coltoggle<-paste(grptoggle,".fact",sep = "")
shapetoggle<-paste(grptoggle,".shapefact",sep = "")
#this plots the chosen data with a 95% confidence ellipse
plotbuild(SlicerMorph.repc[,grptoggle],SlicerMorph.repc[,coltoggle])+
  stat_ellipse(geom="polygon",level = .95,fill=1,alpha=.1, aes(group=SlicerMorph.repc[,shapetoggle]))+
  Pointadd(SlicerMorph.repc[,shapetoggle])+
  plotcolors+plotshapes+
  labs(title = paste(grptoggle,"grouping w/ 95% conf. interval"),x="PC 1 (Var=50.2%)",y="PC 2 (Var=13.3%)",shape="groups")
#ONLY USED FOR BITE AND NOSE
plotbuild(SlicerMorph.repc[,grptoggle],SlicerMorph.repc[,grptoggle])+
  stat_ellipse(geom="polygon",level = .95,fill=1,alpha=.1)+
  Pointadd(SlicerMorph.repc[,grptoggle])+
  plotcolors+plotshapes+
  labs(title = paste(grptoggle,"grouping w/ 95% conf. interval"),x="PC 1 (Var=50.2%)",y="PC 2 (Var=13.3%)",shape="groups")
#this plots without the ellipses
plotbuild(SlicerMorph.repc[,grptoggle],SlicerMorph.repc[,coltoggle])+
  Pointadd(SlicerMorph.repc[,grptoggle])+
  plotcolors+plotshapes+
  labs(title=paste(grptoggle,"grouping"),x="PC 1 (Var=50.2%)",y="PC 2 (Var=13.3%)",color="groups",shape="groups")

#and this is with the ellipses of the k-means data
plotbuild(SlicerMorph.repc[,grptoggle])+
  stat_ellipse(aes(group=factor(km$cluster)),geom="polygon",level = .95,fill=1,alpha=.1)+
  Pointadd(SlicerMorph.repc[,grptoggle])+
  plotcolors+plotshapes+
  labs(title="groupings with k-means cluster 95% confidence intervals",x="PC 1 (Var=50.2%)",y="PC 2 (Var=13.3%)",color="AKC groups",shape="AKC groups")
