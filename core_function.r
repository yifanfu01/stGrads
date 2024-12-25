###DisPlotSpatial 


FindPrimSpot <- function(seurat.obj,celltype_prop,celltype=NULL,probs=0.9){
  st.p3 <- seurat.obj
  if(is.null(celltype))
    celltype=celltype_prop
  
  st.p3@meta.data[,celltype] <- ifelse(st.p3@meta.data[,celltype_prop]>quantile(st.p3@meta.data[,celltype_prop],probs=probs),celltype,'Others')
  # p3 <- SpatialDimPlot(st.p3,group.by = celltype,pt.size.factor = 4e3)
  # p3
  return(st.p3)
}


###CalcNearDis 


CalcNearDis <- function(seurat.obj,celltype,pheno_choose=NULL){
  
  st.p3 <- seurat.obj
  names(st.p3@images) <- 'image1'
  spot_coordinates <- st.p3@images$image1@coordinates
  if(is.null(pheno_choose))
    pheno_choose=celltype
  
  ref_spot = spot_coordinates[c(rownames(st.p3@meta.data[st.p3@meta.data[,celltype]==pheno_choose,])),] 
  query_spot = spot_coordinates[c(rownames(st.p3@meta.data[st.p3@meta.data[,celltype]=='Others',])),] 
  
  
  ref_spot <- ref_spot[ref_spot$tissue==1,]
  query_spot <- query_spot[query_spot$tissue==1,]

  
  calculate_distance <- function(query_col,query_row, ref_col,ref_row) {
    sqrt((query_col - ref_col)^2 + (query_row - ref_row)^2)
  }
  
  find_nearest_ref <- function(query_col,query_row, ref_spot) {
    min.dis=999
    nearest_ref=''
    for( i in 1:length(rownames(ref_spot))){
      distances=calculate_distance(query_col = query_col,query_row = query_row, 
                                   ref_col = ref_spot[i, 'col'],
                                   ref_row = ref_spot[i,'row']
      )
      nearest_ref=ifelse(min.dis>distances,rownames(ref_spot)[i],nearest_ref)
      min.dis=ifelse(min.dis>distances,distances,min.dis)
      # temp.df <- data.frame(ref=NULL,dis=NULL)
      # temp.df <- rbind(temp.df,data.frame(rownames(ref_spot)[i],distances))
      # print(temp.df)
    }
    
    return(c(nearest_ref, min.dis))
  }
  

  nearest_ref_info <- data.frame(query_barcode='',distance='',nearst_barcode='')
  for( i in 1:length(rownames(query_spot))){
    res=find_nearest_ref(query_col = query_spot[i,c('col')],
                         query_row = query_spot[i,c('row')],
                         ref_spot =ref_spot )
    nearest_ref_info <- rbind(nearest_ref_info,data.frame(query_barcode=rownames(query_spot)[i],
                                                          distance=res[2],
                                                          nearst_barcode=res[1]))
  }
  nearest_ref_info <- nearest_ref_info[-1,]
  return(nearest_ref_info)
  
}



###PlotNearDis 


PlotNearDis <- function(seurat.obj,nearest_ref_info,color=NULL,max.dis=20,image.alpha = 0,pt.size.factor = 4e3){
  st.p3 <- seurat.obj
  st.p3$distance=0
  st.p3@meta.data[nearest_ref_info$query_barcode,'distance']=as.numeric(nearest_ref_info$distance)
  st.p3$distance <- as.numeric(st.p3$distance)
  
  if(is.null(color))
    color=c('darkblue','yellow')
  
  
  p4 <- SpatialFeaturePlot(subset(st.p3,distance>0),features = 'distance',max.cutoff = max.dis,
                           image.alpha = image.alpha,pt.size.factor = 4e3)+
    scale_fill_gradient(low = color[1],high = color[2])
  return(p4)
  
}






###PlotDisExpr 

PlotDisExpr <- function(seurat.obj,nearest_ref_info,ref_col,layer='data',assay='SCT',gene,log.trans=F){
  
  st.p3 <- seurat.obj

  st.p3$distance=0
  st.p3@meta.data[nearest_ref_info$query_barcode,'distance']=as.numeric(nearest_ref_info$distance)
  st.p3$distance <- as.numeric(st.p3$distance)
  df <- st.p3@meta.data
  df$target <- GetAssayData(st.p3,assay = assay,layer = layer)[gene,]
  
  # df$dis
  # df$target
  # p6 <- ggplot(df,aes(x=dis,y=target))+geom_point(aes(alpha=I(1/10),color=Duct2_MUC),size=0.75)+
  #   geom_smooth(method = 'lm',se = T)+labs(x="Distance",y=paste0('Relative Expression of ',gene),title = 'Dist-Expr Relationship')+
  #   ggpubr::theme_pubr()
  if(log.trans){
    df$distance <- log(df$distance+1)
  }

  
  p6 <- ggpubr::ggscatter(df,x='distance',y='target',conf.int = T,color =ref_col ,
                          xlab = 'Distance',ylab = paste0('Relative Expression of ',gene),
                          cor.coef = T,add = 'reg.line',size = 0.75,alpha=0.2,
                          add.params = list(color='darkred'))
  return(p6)
}


###PlotDisProp 

PlotDisProp <- function(seurat.obj,nearest_ref_info,ref_col,layer='data',assay='SCT',size=0.8,alpha=0.5,prop_cell,show.gene=NULL,log.trans=F){
  
  st.p3 <- seurat.obj
  
  st.p3$distance=0
  st.p3@meta.data[nearest_ref_info$query_barcode,'distance']=as.numeric(nearest_ref_info$distance)
  st.p3$distance <- as.numeric(st.p3$distance)
  df <- st.p3@meta.data
  if(log.trans){
    df$distance <- log(df$distance+1)
  }
  
  if(is.null(show.gene)){
    p6 <- ggpubr::ggscatter(df,x='distance',y=prop_cell,conf.int = T,color =ref_col ,
                            xlab = 'Distance',ylab = paste0('Proportion of ',prop_cell),cor.method = 'pearson',
                            cor.coef = T,add = 'reg.line',size = size,alpha=alpha,
                            add.params = list(color='darkred'))
  }
  else{
    
    df$target_use_temp <- GetAssayData(st.p3,assay = assay,layer = layer)[show.gene,]
    colnames(df) <- gsub('target_use_temp',show.gene,colnames(df) )
    p6 <- ggpubr::ggscatter(df,x='distance',y=prop_cell,conf.int = T,shape =ref_col ,color=show.gene,
                            xlab = 'Distance',ylab = paste0('Proportion of ',prop_cell),cor.method = 'pearson',
                            cor.coef = T,add = 'reg.line',size = size,alpha=alpha,
                            add.params = list(color='darkred'))+scale_color_gradient(low = 'lightgrey',high = 'darkblue')
  }

  
  
  
  return(p6)
}




