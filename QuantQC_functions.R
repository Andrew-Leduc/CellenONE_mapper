
# Read in raw data
#ev_SCP <- read.delim('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-leduc.an@husky.neu.edu/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users/Oliver/report.tsv')
#ev_QE  <- read.delim('/Users/andrewleduc/Desktop/nPOP_Paper/dat/raw_data/plexDIA/report_plexDIA_mel_nPOP.tsv')

#install.packages('gghighlight')
# Functions

library(BiocManager)
library(psych)
library(dplyr)
library(reshape2)
library(ggplot2)
library(diann)
library(ggpubr)
library(patchwork)
library(genefilter)
library(matrixStats)
library(gghighlight)
library(viridis)
library(ggpointdensity)
library(sva)
library(yaImpute)
library(tidyr)
#library(data.table)


my_col3 <- c("purple2","black")
cv<-function(x){
  
  sd(x, na.rm=T) / mean(x, na.rm=T)
  
}

filt.mat.cr<-function(mat, pct.r,pct.c){
  
  kc<-c()
  for(k in 1:ncol(mat)){
    
    pct.na<-length(which(is.na(mat[,k]))) / length(mat[,k])
    if(pct.na <= pct.c){ kc<-c(kc,k)}
    #print(pct.na)
    
    
  }
  
  mat<-mat[,kc]
  
  kr<-c()
  for(k in 1:nrow(mat)){
    
    pct.na<-length(which(is.na(mat[k,]))) / length(mat[k,])
    if(pct.na <= pct.r){ kr<-c(kr,k)}
    #print(pct.na)
    
    
  }
  
  mat<-mat[kr,]
  
  
  
  return(mat)
  
}

hknn<-function(dat, k){
  
  # Create a copy of the data, NA values to be filled in later
  dat.imp<-dat
  
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat<-as.matrix( dist(t(dat)) )
  #dist.mat<- 1-as.matrix(cor((dat), use="pairwise.complete.obs"))
  
  #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))
  
  # Column names of the similarity matrix, same as data matrix
  cnames<-colnames(dist.mat)
  
  # For each column in the data... 
  for(X in cnames){
    
    # Find the distances of all other columns to that column 
    distances<-dist.mat[, X]
    
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered<-distances[order(distances, decreasing = F)]
    
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered<-dat[ , names(distances.ordered ) ]
    
    # Take the values in the column of interest
    vec<-dat[, X]
    
    # Which entries are missing and need to be imputed...
    na.index<-which( is.na(vec) )
    
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      
      # Find the most similar columns that have a non-NA value in this row
      closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )
      
      #print(length(closest.columns))
      
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
        
      }
      
      
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i]<-mean( dat[ i, closest.columns ] )
        
      }
      
      
    }
    
    # Populate a the matrix with the new, imputed values
    dat.imp[,X]<-vec
    
  }
  
  return(dat.imp)
  
}


Normalize <- function(evnew,log = 'no'){
  evnew <- as.matrix(evnew)
  evnew[evnew==0] <- NA
  for(i in 1:ncol(evnew)){
    evnew[,i] <- evnew[,i]/median(evnew[,i],na.rm = T)
  }
  for(i in 1:nrow(evnew)){
    evnew[i,] <- evnew[i,]/mean(evnew[i,],na.rm = T)
  }
  if(log == 'yes'){
    evnew <- log2(evnew)
  }
  
  return(evnew)
}
Normalize_saad <- function(dat,log = 'no'){
  dat <- as.matrix(dat)
  dat[dat==0] <- NA

  refVec <- rowMedians(x = dat, cols = 1:ncol(dat), na.rm = T)
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k] * median(refVec/dat[,k], na.rm = T)
  }
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]/mean(dat[k,], na.rm = T)
  }
  if(log == 'yes'){
    dat <- log2(dat)
  }
  return(dat)
}

Normalize_saad_log <- function(dat){
  refVec <- rowMedians(x = dat, cols = 1:ncol(dat), na.rm = T)
  for(k in 1:ncol(dat)){
    dat[,k]<-dat[,k] + median(refVec - dat[,k], na.rm = T)
  }
  for(k in 1:nrow(dat)){
    dat[k,]<-dat[k,]-mean(dat[k,], na.rm = T)
  }
  return(dat)
}


CVs <- function(mat,thresh,cell_id){
  abs_mat <- melt(mat,id = c('Protein.Group','seqcharge'))
  rownames(mat) <- mat$seqcharge
  mat_cv <- mat
  mat_cv$Protein.Group <- NULL
  mat_cv$seqcharge <- NULL
  mat_cv <- as.matrix(mat_cv)     
  mat_norm <- Normalize_saad(mat_cv)
  evnew_CVs <- as.data.frame(mat_norm)
  evnew_CVs$Protein <- mat$Protein.Group
  evnew_CVs$pep <- rownames(evnew_CVs)
  evnew_CVs <- melt(evnew_CVs, id = c('Protein','pep'))
  
  suppressMessages({
  abs_mat <- abs_mat %>%
    dplyr::group_by(Protein.Group, variable) %>%
    dplyr::summarise(cvq = mean(value,na.rm=T))
  
  CV_mat <- evnew_CVs %>%
    dplyr::group_by(Protein, variable) %>%
    dplyr::summarise(cvq = cv(value),ct = sum(is.na(value)==F))
  })
  CV_mat2 <- CV_mat
  CV_mat2$intense <- abs_mat$cvq
  CV_mat_plot <- CV_mat2 %>% filter(ct > 2)
  CV_mat_plot <- CV_mat_plot %>% filter(intense != Inf)
  CV_mat_plot <- CV_mat_plot %>% filter(intense != -Inf)
  CV_mat_plot <- CV_mat_plot %>% filter(is.nan(intense) == F)
  CVs_vs_intensity <- ggplot(CV_mat_plot,aes(x=log10(intense),y=cvq)) + 
    geom_pointdensity() + scale_color_viridis()+
    ylim(c(0,2))+xlim(c(min(log10(CV_mat_plot$intense)),max(log10(CV_mat_plot$intense))))+
    theme_bw() +ylab('Protein CVs')+xlab('log10(Protein intensity)')+
    font("xylab", size=30) +
    font("xy.text", size=20) + rremove('legend')
  
  CV_mat <- CV_mat %>%
    dplyr::group_by( variable) %>%
    dplyr::summarise(cvq = median(cvq,na.rm = T))
  
  pos <- cell_id %>% filter(sample != 'neg')
  
  CV_mat$value <- NA
  CV_mat$value[CV_mat$variable %in% pos$ID] <- 'cell'
  CV_mat$value[!CV_mat$variable %in% pos$ID] <- 'neg'
  
  CV_mat_pos <- CV_mat %>% filter(value == 'cell')
  CV_mat_neg <- CV_mat %>% filter(value == 'neg')
    
  CV_plot <- ggplot(data=CV_mat, aes(x=cvq,fill=value)) + geom_density( alpha=0.5,adjust=1.5) + theme_pubr() +
    scale_fill_manual(values=my_col3) + 
    xlab("CV") + ylab("Fraction of cells") + rremove("y.ticks") + rremove("y.text") +
    font("xylab", size=30) +
    font("x.text", size=20) +
    coord_cartesian(xlim=c(.1,.8))+
    annotate("text", x=0.2, y= 14, label=paste0(sum(CV_mat_pos$cvq < thresh)," cells"), size=10, color=my_col3[c(1)])+
    annotate("text", x=0.64, y= 12, label=paste0(sum(CV_mat_neg$cvq > thresh)," Ctr -"), size=10, color=my_col3[c(2)])+
    annotate("text", x=0.63, y= 14, label=paste0(sum(CV_mat_pos$cvq > thresh)," cells"), size=10, color=my_col3[c(1)])+
    annotate("text", x=0.2, y= 12, label=paste0(sum(CV_mat_neg$cvq < thresh)," Ctr -"), size=10, color=my_col3[c(2)])+
    #annotate("text", x=0.25, y= 3, label="Macrophage-like", size=6) +
    rremove("legend") + geom_vline(xintercept=thresh, lty=2, size=2, color="gray50") + theme(plot.margin = margin(1, 1, 0, 1, "cm"))
  
  CV_mat <- CV_mat %>% filter(value == 'cell')
  return(list(CVs_vs_intensity,CV_plot,CV_mat))
  
}
Pep_cor <- function(mat,prot_list){
  #mat = evnew_norm
  #prot_list <- protein_list
  mat <- as.matrix(mat)
  count = 0
  mat_stor = matrix(data = NA,nrow = 10000,ncol = 2)
  mat_stor_fake = matrix(data = NA,nrow = 10000,ncol = 2)
  upl <- unique(prot_list)
  for(i in upl){
    mat_p1 <- mat[which(prot_list == i),]
    if(is.null(nrow(mat_p1)) == FALSE ){
      obs_mat <- psych::pairwiseCount(t(mat_p1))
      cor_mat <- cor(t(mat_p1),use = 'pairwise.complete.obs')
      cor_mat[obs_mat < 4] <- NA
      cor_mat <- cor_mat[lower.tri(cor_mat)]
      if(is.na(median(cor_mat,na.rm = T)) == F){
        count = count + 1
        mat_stor[count,1] <- i
        mat_stor[count,2] <- median(cor_mat,na.rm = T)
      }
    }
  }
  prot_list <- sample(prot_list)  
  for(i in upl){
    mat_p1 <- mat[which(prot_list == i),]
    if(is.null(nrow(mat_p1)) == FALSE ){
      obs_mat <- psych::pairwiseCount(t(mat_p1))
      cor_mat <- cor(t(mat_p1),use = 'pairwise.complete.obs')
      cor_mat[obs_mat < 4] <- NA
      cor_mat <- cor_mat[lower.tri(cor_mat)]
      if(is.na(median(cor_mat,na.rm = T)) == F){
        count = count + 1
        mat_stor_fake[count,1] <- i
        mat_stor_fake[count,2] <- median(cor_mat,na.rm = T)
      }
    }
  }
  
  mat_stor <- as.data.frame(mat_stor)
  colnames(mat_stor) <- c('Protein','Cor')
  mat_stor <- mat_stor %>% filter(is.na(Cor) == F)
  mat_stor$Cor <- as.numeric(mat_stor$Cor)
  mat_stor$Cor <- mat_stor$Cor - median(as.numeric(mat_stor_fake[,2]),na.rm = T)
  pep_cor <- ggplot(mat_stor,aes(x = Cor)) + geom_histogram() + ylab('# of Proteins')+xlab('Median Cor between peptides mapping to a protein')
  return(mat_stor)
}

analyzeCellenONE <- function(allDays,pickupPath2,pickupPath1,labelPath){
    
    allDays[grepl("Transmission",allDays$X),]$X <- NA
    allDays <- allDays %>% fill(2:7, .direction = "up") %>% drop_na(XPos)
    
    #### Labelling and Pickup Field Files
    ## labelling file
    
    con_lab <-file(labelPath)
    lines_lab<- readLines(con_lab)
    close(con_lab)
    slines_lab <- strsplit(lines_lab,"\t")
    colCount_lab <- max(unlist(lapply(slines_lab, length)))
    
    label <- read.table(labelPath, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)
    
    fieldRaw_label <- label[grepl("\\[\\d",label$position),]$position
    fieldBrack_label <- gsub("\\[|\\]", "", fieldRaw_label)
    fieldComma_label <- strsplit(fieldBrack_label, "," )
    fieldNum_label <- unlist(lapply(fieldComma_label, `[[`, 1))
    label$field[grepl("\\[\\d",label$position)] <- fieldNum_label
    label <- label %>% fill(field, .direction = "down")
    label <-  label[!label$well=="",]
    label$field <- as.numeric(label$field) + 1
    
    labelxyPos <- strsplit(label$position, "\\/")
    label$yPos <- unlist(lapply(labelxyPos, '[[', 1))
    label$xPos <- unlist(lapply(labelxyPos, '[[', 2))
    
    ## sample pickup file
    
    con_pickup <-file(pickupPath1)
    lines_pickup<- readLines(con_pickup)
    close(con_pickup)
    slines_pickup <- strsplit(lines_pickup,"\t")
    colCount_pickup <- max(unlist(lapply(slines_pickup, length)))
    
    pickup <- read.table(pickupPath1, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 22)
    
    
    
    con_pickup <-file(pickupPath2)
    lines_pickup<- readLines(con_pickup)
    close(con_pickup)
    slines_pickup <- strsplit(lines_pickup,"\t")
    colCount_pickup <- max(unlist(lapply(slines_pickup, length)))
    pickup2 <- read.table(pickupPath2, sep="\t",fill=TRUE,header = F,col.names=c("position", "well","volume","field"), quote = "", skip = 27)
    
    pickup <- rbind(pickup,pickup2)
    
    fieldRaw_pickup <- pickup[grepl("\\[\\d",pickup$position),]$position
    fieldBrack_pickup <- gsub("\\[|\\]", "",fieldRaw_pickup)
    fieldComma_pickup <- strsplit(fieldBrack_pickup, "," )
    fieldNum_pickup <- unlist(lapply(fieldComma_pickup, `[[`, 1))
    pickup$field[grepl("\\[\\d",pickup$position)] <- fieldNum_pickup
    pickup <- pickup %>% fill(field, .direction = "down")
    pickup <-  pickup[!pickup$well=="",]
    pickup$field <- as.numeric(pickup$field) + 1
    
    
    
    pickupxyPos <- strsplit(pickup$position, "\\/")
    pickup$yPos <- unlist(lapply(pickupxyPos, '[[', 1))
    pickup$xPos <- unlist(lapply(pickupxyPos, '[[', 2))
    
    
    ### making sure that label "slide" refers to fields or actual slides.
    ## There is a sequence of 108 x,y and there are 108 spots per slide
    # order is y/x positions
    LabelxyPos <- strsplit(label$position, "\\/")
    label$yPos <- unlist(lapply(LabelxyPos, '[[', 1))
    label$xPos <- unlist(lapply(LabelxyPos, '[[', 2))
    label$well <- substring(label$well, 3)
    label$well <- as.numeric(gsub(",","",label$well))
    matchTMTSCP <- paste0("TMT", 1:length(unique(label$well)))
    
    for (i in 1:length(matchTMTSCP)) {
      
      label[which(label$well == i),]$well <- matchTMTSCP[i]
      
    }
    
    label$well <- substring(label$well, 4)
    
    #### Trying to map label to cell
    ###  lets try and keep all three important in one
    allDays <- transform(allDays, xyf = paste0(XPos, YPos, Field))
    label <- transform(label, xyf = paste0(xPos, yPos, field))
    
    ### Isolation and Label merged
    
    isoLab <- allDays %>% group_by(Target) %>% left_join(label, by = 'xyf')
    labelCount <- isoLab %>% group_by(well) %>% dplyr::summarize(count=n())
    
    ## trying to get pickup working in the same way
    # fixing pickup file
    pickup$target <- 1
    
    pickup$well <- substring(pickup$well, 2)
    pickup$well <- gsub(",","",pickup$well)
    
    slidesUsedForPrep <- unique(isoLab$Target)
    
    #pickup <- pickup[which(pickup$target %in% slidesUsedForPrep),]
    
    ### clustering together pickup points with sample points using ANN
    ## super convolutedly/idiotically
    
    isoLab$ann <- NA
    isoLab$pickupX <- NA
    isoLab$pickupY <- NA
    
    ann_123 <- ann(ref = as.matrix(unique(pickup[(-which(pickup$field == 4)) , c("xPos","yPos")])),  target = as.matrix(isoLab[(-which(isoLab$Field == 4)) , c("xPos","yPos")]), k=1)
    
    ann_4 <- ann(ref = as.matrix(unique(pickup[(which(pickup$field == 4)) , c("xPos","yPos")])),  target = as.matrix(isoLab[(which(isoLab$Field == 4)) , c("xPos","yPos")]), k=1)
    
    isoLab[(-which(isoLab$Field == 4)),]$ann <-  ann_123$knnIndexDist[,1]
    isoLab[(which(isoLab$Field == 4)),]$ann <-  ann_4$knnIndexDist[,1]
    
    ## split - combine
    isoLab_123 <- isoLab[-which(isoLab$Field == 4),]
    isoLab_4 <- isoLab[which(isoLab$Field == 4),]
    
    
    notFieldFourPickUnique <- unique(pickup[-(which(pickup$field == 4)) , c("xPos","yPos")])
    fieldFourPickUnique <- unique(pickup[(which(pickup$field == 4)) , c("xPos","yPos")])
    
    isoLab_123$pickupX <- notFieldFourPickUnique[isoLab_123$ann,]$xPos
    isoLab_123$pickupY <- notFieldFourPickUnique[isoLab_123$ann,]$yPos
    
    isoLab_4$pickupX <- fieldFourPickUnique[isoLab_4$ann,]$xPos
    isoLab_4$pickupY <- fieldFourPickUnique[isoLab_4$ann,]$yPos
    isoLab_bound <- rbind(isoLab_123, isoLab_4)
    
    
    ### Merge pickup and isoLab
    isoLab_bound <- transform(isoLab_bound, xyft = paste0(pickupX, pickupY, Field, Target))
    pickup <-  transform(pickup, xyft = paste0(xPos, yPos, field, target))
    
    intersect(pickup$xyft,isoLab_bound$xyft)
    
    isoLab_final <- isoLab_bound %>% left_join(pickup, by = 'xyft')
    wellCount <- isoLab_final %>% group_by(well.y) %>% dplyr::summarize(count=n())
    
    
    
    ### Clean up to yield final dataframe
    cellenOne_data <- data.frame (sample = isoLab_final$condition, isoTime = isoLab_final$Time, diameter = isoLab_final$Diameter, elongation = isoLab_final$Elongation, slide = isoLab_final$Target, field = isoLab_final$Field, dropXPos = isoLab_final$XPos, dropYPos = isoLab_final$YPos, label = isoLab_final$well.x, pickupXPos = isoLab_final$pickupX, pickupYPos = isoLab_final$pickupY, injectWell = isoLab_final$well.y)
    
    
    ##*sigh* not done yet
    
    cellenOne_data$wellAlph <- substring(cellenOne_data$injectWell, 1, 1)
    cellenOne_data$wellNum <- substring(cellenOne_data$injectWell, 2)
    
    cellenOne_data <- cellenOne_data %>% arrange(wellAlph, as.numeric(wellNum), as.numeric(label))
    cellenOne_data <- cellenOne_data %>% dplyr::select(dropYPos,dropXPos,sample,diameter,elongation,field,label,injectWell)
    return(cellenOne_data)
}




