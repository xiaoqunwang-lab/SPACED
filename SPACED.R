#This method requires the following parameters as input:
#                 1.degs -> differential expressed genes
#                 2.ident_use -> cell identity class
#                 3.celltype ->  cell identity


#function prop_calc calculates probability of gene expression in a given cell type
prop_calc <- function(degs,ident_use,celltype){
  library(Seurat)
  prop <- matrix(0, nrow(degs), 1)
  dim(prop)
  pos <- grep(celltype[1],ex_neuron@meta.data[ident_use][[1]])
  length(pos)
  
  for(i in 1:nrow(degs)){
    for(j in 1:ncol(prop)){
      prop[i,j] <-  length(grep('TRUE',as.matrix(ex_neuron@assays$RNA@counts[degs$gene[i],pos]>0)))/length(pos)
  }
    
    }
  for(i in 2:length(celltype)){
    tmp <- matrix(0, nrow(degs), 1)
    dim(tmp)
    pos <- grep(celltype[i],ex_neuron@meta.data[ident_use][[1]])
    length(pos)
    
    for(i in 1:nrow(degs)){
      for(j in 1:ncol(tmp)){
        tmp[i,j] <-  length(grep('TRUE',as.matrix(ex_neuron@assays$RNA@counts[degs$gene[i],pos]>0)))/length(pos)  
      }
    } 
    prop <- cbind(prop,tmp)
    print(celltype)
  }
  colnames(prop)<- celltype
  rownames(prop)<- degs$gene
  return(prop)
}

#test <- prop_calc(ex_neuron_marker[1:20,],'merge_id_ed1',unique(ex_neuron$merge_id_ed1))

makeprobsvec<-function(p){ 
  phat<-p/sum(p) 
  phat[is.na(phat)] = 0 
  phat 
} 

shannon.entropy <- function(p) { 
  if (min(p) < 0 || sum(p) <=0) 
    return(Inf) 
  p.norm<-p[p>0]/sum(p) 

  -sum( log2(p.norm)*p.norm) 
} 


JSdistVec <- function(p, q){ 
  JSdiv <- shannon.entropy((p + q)/2) - (shannon.entropy(p) +  
                                           shannon.entropy(q)) * 0.5 
  JSdiv[is.infinite(JSdiv)] <- 1 
  JSdiv[JSdiv < 0] <- 0 
  JSdist <- sqrt(JSdiv) 
  JSdist 
} 

specificity_scorer = function(normpropmat){ 
  marker_gene_specificities <- lapply(1:ncol(normpropmat), function(cell_type_i){ 
    perfect_specificity <- rep(0.0, ncol(normpropmat)) 
    perfect_specificity[cell_type_i] <- 1.0 
    apply(normpropmat, 1, function(x) {  
      if (sum(x) > 0) 1 - JSdistVec(makeprobsvec(x), perfect_specificity) 
      else 0 
    }) 
  }) 
  return(do.call(cbind, marker_gene_specificities)) 
} 

#function specifcity_score calculates the specificity 
#of gene expression among given cell types 
specificity_score <- function(prop){
  deg_specificities = specificity_scorer(prop)
  colnames(deg_specificities)<- colnames(prop)
  
  deg_specificities <- as.data.frame(deg_specificities)
  
  return(deg_specificities)
}

#function top_n_spec_deg calculates the top n cell type
#specific genes
top_n_spec_deg <- function(deg_specificities,n){
  library(reshape2)
  library(dplyr)
  deg_specificities$gene = rownames(deg_specificities)
  deg_specificities = melt(deg_specificities,id='gene')
  deg_specificities = deg_specificities[order(deg_specificities$value,decreasing = T),]
  deg_specificities = arrange(deg_specificities,value,.by_group = T)
  deg_specificities = deg_specificities%>%group_by(variable)%>%arrange(desc(value),.by_group = T)
  topn = deg_specificities%>%group_by(variable)%>%top_n(n,value)
  return(topn)
}

#function SPACED calculates regional specificity for a given subtype
SPACED <- function(insitu_intensity,celltype){
  library(matrixStats)
  insitu_intensity_norm <- matrix(0, nrow(insitu_intensity), ncol(insitu_intensity))
  for(j in 1:ncol(insitu_intensity)){
    min_tmp =min(insitu_intensity[,j])
    max_tmp = max(insitu_intensity[,j])
    for(i in 1:nrow(insitu_intensity)){
      insitu_intensity_norm[i,j]<- ((insitu_intensity[i,j]-min_tmp)/(max_tmp-min_tmp))+1
    }
  }
  for(i in 1:ncol(insitu_intensity_norm)){
    insitu_intensity_norm[,i]<- insitu_intensity_norm[,i]*colSds(as.matrix(insitu_intensity))[i]
  }
  #print(insitu_intensity_norm)
  
  colnames(insitu_intensity_norm)<- colnames(insitu_intensity)
  rownames(insitu_intensity_norm) <- rownames(insitu_intensity)
  insitu_intensity_norm <- data.frame(SuG=mean(insitu_intensity_norm[1,]),OP=mean(insitu_intensity_norm[2,]),
                                      InG=mean(insitu_intensity_norm[3,]),DPG=mean(insitu_intensity_norm[4,]))
  rownames(insitu_intensity_norm)<-celltype
  insitu_intensity_norm_binary <- matrix(0,nrow(insitu_intensity_norm),ncol(insitu_intensity_norm))
  print(insitu_intensity_norm_binary)
  
  for(i in 1:nrow(insitu_intensity_norm)){
    insitu_intensity_norm_binary[i,as.numeric(which.max(insitu_intensity_norm[i,]))]<-insitu_intensity_norm[i,which.max(insitu_intensity_norm[i,])]
  }
  rownames(insitu_intensity_norm_binary)<- rownames(insitu_intensity_norm)
  colnames(insitu_intensity_norm_binary)<- colnames(insitu_intensity_norm)
  return(insitu_intensity_norm)
}

#function SPACED_binary calculates binarized regional specificity for a given subtype
#required input: 1.in situ signal intensity accessed from Allen Institute
#                (https://mouse.brain-map.org/search/index)
#                2.cell type, for which regional specificity should be calculated
SPACED_binary <- function(insitu_intensity,celltype){
  library(matrixStats)
  insitu_intensity_norm <- matrix(0, nrow(insitu_intensity), ncol(insitu_intensity))
  for(j in 1:ncol(insitu_intensity)){
    min_tmp =min(insitu_intensity[,j])
    max_tmp = max(insitu_intensity[,j])
    for(i in 1:nrow(insitu_intensity)){
      insitu_intensity_norm[i,j]<- ((insitu_intensity[i,j]-min_tmp)/(max_tmp-min_tmp))+1
    }
  }
  for(i in 1:ncol(insitu_intensity_norm)){
    insitu_intensity_norm[,i]<- insitu_intensity_norm[,i]*colSds(as.matrix(insitu_intensity))[i]
  }
  #print(insitu_intensity_norm)
  
  colnames(insitu_intensity_norm)<- colnames(insitu_intensity)
  rownames(insitu_intensity_norm) <- rownames(insitu_intensity)
  insitu_intensity_norm <- data.frame(SuG=mean(insitu_intensity_norm[1,]),OP=mean(insitu_intensity_norm[2,]),
                                          InG=mean(insitu_intensity_norm[3,]),DPG=mean(insitu_intensity_norm[4,]))
  rownames(insitu_intensity_norm)<-celltype
  insitu_intensity_norm_binary <- matrix(0,nrow(insitu_intensity_norm),ncol(insitu_intensity_norm))
  #print(insitu_intensity_norm_binary)
  
  for(i in 1:nrow(insitu_intensity_norm)){
    insitu_intensity_norm_binary[i,as.numeric(which.max(insitu_intensity_norm[i,]))]<-insitu_intensity_norm[i,which.max(insitu_intensity_norm[i,])]
  }
  rownames(insitu_intensity_norm_binary)<- rownames(insitu_intensity_norm)
  colnames(insitu_intensity_norm_binary)<- colnames(insitu_intensity_norm)
  return(insitu_intensity_norm_binary)
}


