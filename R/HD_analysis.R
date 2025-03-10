

SeuLocation_StConvert <- function(seurat_obj,assay='Spatial',
                                  file.st,
                                  ref.name){
  library(DropletUtils)
  library(tibble)
  library(Seurat)
  library(jsonlite)
  library(stringr)
  
  st.r <- seurat_obj
  st.r
  st.spatial <- GetAssayData(st.r,assay = 'Spatial')
  st.use <- CreateSeuratObject(st.spatial,assay = 'Spatial')
  st.use@meta.data <- st.r@meta.data
  st.use@images <- st.r@images
  names(st.use@images) <- 'slice1'
  #names(st.use@images)
  st.use <- SCTransform(st.use, assay = "Spatial")  
  
  # colnames(st.use)
  # colnames(st.use) <- gsub('sample_12L_gongmo:','S12L-',colnames(st.use))
  
  
  
  DropletUtils::write10xCounts(x = GetAssayData(st.r,assay = assay,layer = 'counts'), 
                               #barcodes = str_split(colnames(st.use@assays$Spatial$counts),'_',simplify = TRUE)[,3] ,
                               gene.symbol = rownames(GetAssayData(st.r,assay = assay,layer = 'counts')), 
                               gene.type = "Gene Expression", 
                               overwrite = TRUE, type = 'HDF5', 
                               genome = "CRCh38", version = "3", 
                               chemistry = "Single Cell 3' v1", 
                               library.ids = "slide_1", ### library.ids is key, and use"adata_vis = sc.read_visium(spatial,library_id='slice1')" in py
                               path = paste0(file.st,'/filtered_feature_bc_matrix.h5'))
  
  
  low_factor <- st.use@images$slice1@scale.factors$lowres
  hi_factor <- st.use@images$slice1@scale.factors$hires
  
  
  coordinates <- st.use@images$slice1@coordinates
  coordinates <- rownames_to_column(coordinates)
  # coordinates[,1] <- str_split(coordinates[,1],'_',simplify = TRUE)[,2] #可以用这行修改barcodes前缀
  write.table(coordinates,file =paste0 (file.st,'/spatial/tissue_positions_list.csv'),
              quote = F,row.names = F,col.names = F,sep=',')
  
  scale_factors <- data.frame(spot_diameter_fullres=st.use@images$slice1@scale.factors$spot,
                              tissue_hires_scalef=st.use@images$slice1@scale.factors$hires,
                              fiducial_diameter_fullres=st.use@images$slice1@scale.factors$fiducial,
                              tissue_lowres_scalef=st.use@images$slice1@scale.factors$lowres)
  
  jsondata <- toJSON(scale_factors)
  cat(str_remove_all(jsondata, '[\\[\\]]'), 
      file = paste0(file.st,'/spatial/scalefactors_json.json'), 
      fill = FALSE, labels = NULL, append = FALSE)
}


FindPrimSpot <- function(seurat.obj,celltype_prop,celltype=NULL,probs=0.9){
  st.p3 <- seurat.obj
  if(is.null(celltype))
    celltype=celltype_prop
  
  st.p3@meta.data[,celltype] <- ifelse(st.p3@meta.data[,celltype_prop]>
                                         quantile(st.p3@meta.data[,celltype_prop],probs=probs),
                                       celltype,'Others')

  return(st.p3)
}




calculate_distance <- function(query_col,query_row, ref_col,ref_row) {
  sqrt((query_col - ref_col)^2 + (query_row - ref_row)^2)
}




find_nearest_ref <- function(query_col,query_row, ref_spot,max.r=50,model='Linear') {
  if(is.null(max.r))
  {max.r=999}
  min.dis=max.r+1
  nearest_ref='None'
  ref_spot.filt <- ref_spot[abs(ref_spot$row-query_row)<=max.r & abs(ref_spot$col-query_col)<=max.r,]
  strength.sum=0
  strength=0
  if(length(rownames(ref_spot.filt))!=0){
    
    for( i in 1:length(rownames(ref_spot.filt))){
      distances=calculate_distance(query_col = query_col,query_row = query_row, 
                                   ref_col = ref_spot.filt[i, 'col'],
                                   ref_row = ref_spot.filt[i,'row']
      )
      if(model=='Linear')
        strength.sum=strength.sum+1/distances ###Linear
      if(model=='Exp')
        strength.sum=strength.sum+exp((-1)*distances) ###Exp
      if(model=='Lg')
        strength.sum=strength.sum+1/log(distances+1)  ###Lg
      
      nearest_ref=ifelse(min.dis>distances,rownames(ref_spot.filt)[i],
                         ifelse(min.dis==distances,paste0(nearest_ref,',',rownames(ref_spot.filt)[i]),nearest_ref))
      min.dis=ifelse(min.dis>distances,distances,min.dis)
      
      # temp.df <- data.frame(ref=NULL,dis=NULL)
      # temp.df <- rbind(temp.df,data.frame(rownames(ref_spot)[i],distances))
      # print(temp.df)
    }
    strength=ifelse(min.dis<=max.r,length(str_split(nearest_ref,pattern = ',')[[1]]),0)
  }
  
  return(c(nearest_ref, min.dis,strength,strength.sum))
}



find_nearest_ref_HD <- function(query_col,query_row, ref_spot,max.r=50,model='Linear') {
  library(stringr)
  if(is.null(max.r))
  {max.r=999}
  min.dis=max.r+1
  nearest_ref='None'
  ref_spot.filt <- ref_spot[abs(ref_spot[,2]-query_row)<=max.r & abs(ref_spot[,1]-query_col)<=max.r,]
  strength.sum=0
  strength=0
  if(length(rownames(ref_spot.filt))!=0){
    
    for( i in 1:length(rownames(ref_spot.filt))){
      distances=calculate_distance(query_col = query_col,query_row = query_row, 
                                   ref_col = ref_spot.filt[i, 'x'],
                                   ref_row = ref_spot.filt[i,'y']
      )
      if(model=='Linear')
        strength.sum=strength.sum+1/distances ###Linear
      if(model=='Exp')
        strength.sum=strength.sum+exp((-1)*distances) ###Exp
      if(model=='Lg')
        strength.sum=strength.sum+1/log(distances+1)  ###Lg
      
      nearest_ref=ifelse(min.dis>distances,rownames(ref_spot.filt)[i],
                         ifelse(min.dis==distances,paste0(nearest_ref,',',rownames(ref_spot.filt)[i]),nearest_ref))
      min.dis=ifelse(min.dis>distances,distances,min.dis)
      
      # temp.df <- data.frame(ref=NULL,dis=NULL)
      # temp.df <- rbind(temp.df,data.frame(rownames(ref_spot)[i],distances))
      # print(temp.df)
    }
    strength=ifelse(min.dis<=max.r,length(str_split(nearest_ref,pattern = ',')[[1]]),0)
  }
  
  return(c(nearest_ref, min.dis,strength,strength.sum))
}


###CalcNearDis 计算最近距离矩阵
#seurat.obj 输入对象，只可以包含一张图。
#celltype是对应的细胞标签列名,且该列默认为同样的
#pheno_choose, 若上一列表，列名内选择不相同，则pheno_choose为选择的阳性，默认NULL，与celltype相同
#max.r 是传入find_nearest_ref的max.r 通常表示距离，默认Visium=10
#calc.strength 是否计算距离源强度，默认false
#model是强度衰减模型
#返回一个数据框

CalcNearDis <- function(seurat.obj,celltype,pheno_choose=NULL,calc.strength=FALSE,model='Linear',max.r=10){
  
  st.p3 <- seurat.obj
  names(st.p3@images) <- 'image1'
  spot_coordinates <- st.p3@images$image1@coordinates
  if(is.null(pheno_choose))
    pheno_choose=celltype
  
  ref_spot = spot_coordinates[c(rownames(st.p3@meta.data[st.p3@meta.data[,celltype]==pheno_choose,])),] 
  query_spot = spot_coordinates[c(rownames(st.p3@meta.data[st.p3@meta.data[,celltype]=='Others',])),] 
  
  
  ref_spot <- ref_spot[ref_spot$tissue==1,]
  query_spot <- query_spot[query_spot$tissue==1,]
  
  
  
  
  # find_nearest_ref(query_col = 43,query_row = 3,ref_spot = ref_spot)
  # #"PDACP_8_GTCTATCTGAGTTTCT-1" "12"  
  # SpatialDimPlot(st.p3,cells.highlight =list (query='PDACP_8_AAACTGCTGGCTCCAA-1',
  #                                             ref_nearst='PDACP_8_TACGACTGCCTCTTAG-1',
  #                                             all_other_ref=rownames(ref_spot)[rownames(ref_spot)!='PDACP_8_TACGACTGCCTCTTAG-1']
  # ),
  # cols.highlight = c('darkred','darkgreen','darkblue','white')
  # ,pt.size.factor = 4e3)
  if(!calc.strength){
    nearest_ref_info <- data.frame(query_barcode='',distance='',nearst_barcode='')
    for( i in 1:length(rownames(query_spot))){
      res=find_nearest_ref(query_col = query_spot[i,c('col')],
                           query_row = query_spot[i,c('row')],
                           ref_spot =ref_spot ,max.r = max.r)
      nearest_ref_info <- rbind(nearest_ref_info,data.frame(query_barcode=rownames(query_spot)[i],
                                                            distance=res[2],
                                                            nearst_barcode=res[1]))
      
    }
    nearest_ref_info <- nearest_ref_info[-1,]
  }
  else{
    nearest_ref_info <- data.frame(query_barcode='',distance='',nearst_barcode='',strength=0,strength.sum=0)
    for( i in 1:length(rownames(query_spot))){
      res=find_nearest_ref(query_col = query_spot[i,c('col')],
                           query_row = query_spot[i,c('row')],
                           ref_spot =ref_spot ,max.r = max.r,model=model)
      nearest_ref_info <- rbind(nearest_ref_info,data.frame(query_barcode=rownames(query_spot)[i],
                                                            distance=res[2],strength=res[3],strength.sum=res[4],
                                                            nearst_barcode=res[1]))
      
    }
    nearest_ref_info$strength <- as.numeric(nearest_ref_info$strength)
    nearest_ref_info$strength.sum <- as.numeric(nearest_ref_info$strength.sum)
    nearest_ref_info <- nearest_ref_info[-1,]
  }
  
  
  return(nearest_ref_info)
  
}


#HD函数可以不使用Find来定义primary spot，而是直接指定特定的细胞
#resolution 默认8um，用来归一化物理距离
#query_celltype是特定的查询亚群,NULL则为所有的亚群

CalcNearDis_HD <- function(seurat.obj,celltype,query_celltype=NULL,
                           pheno_choose=NULL,calc.strength=FALSE,model='Linear',
                           resolution=8,max.r=10){
  
  sce1.cp <- seurat.obj
  names(sce1.cp@images) <- 'image1'
  spot_coordinates <-sce1.cp@images[["image1"]]@boundaries[["centroids"]]@coords
  spot_coordinates <- spot_coordinates/sce1.cp@images[["image1"]]@scale.factors$spot
  spot_coordinates <- spot_coordinates*resolution
  
  if(!identical(length(spot_coordinates[,1]),length(colnames(sce1.cp))))
    print("!warning: Barcode corrected error!")
  rownames(spot_coordinates) <- colnames(sce1.cp)
  
  # x y
  if(is.null(pheno_choose))
    pheno_choose=celltype
  
  ref_spot = spot_coordinates[c(rownames(sce1.cp@meta.data[ !is.na(sce1.cp@meta.data[,celltype]) &
                                                              sce1.cp@meta.data[,celltype]==pheno_choose,])),] 
  print(paste0('Ref_spots numbers is:',length(rownames(ref_spot))))
  if(is.null(query_celltype)){
    query_spot =  spot_coordinates[c(rownames(sce1.cp@meta.data[ !is.na(sce1.cp@meta.data[,celltype]) &
                                                                   sce1.cp@meta.data[,celltype]!=pheno_choose,])),] 
  }
  else{
    query_spot =  spot_coordinates[c(rownames(sce1.cp@meta.data[ !is.na(sce1.cp@meta.data[,celltype]) &
                                                                   sce1.cp@meta.data[,celltype]==query_celltype,])),] 
  }
  
  
  
  
  print(paste0('Query_spots numbers is:',length(rownames(query_spot))))
  # ref_spot <- ref_spot[ref_spot$tissue==1,]
  # query_spot <- query_spot[query_spot$tissue==1,]
  
  
  # x(col) y (row)
  
  #s_008um_00017_00821-1 4971.317 10048.944
  
  # find_nearest_ref_HD(query_col =4971.317 ,query_row = 10048.944,ref_spot = ref_spot,max.r = 1000)
  # #"s_008um_00032_00806-1" "155.076633447315"      "1"                     "0.314818629754032"    
  # SpatialDimPlot((sce1.cp),cells.highlight =list (query='s_008um_00017_00821-1',
  #                                             ref_nearst='s_008um_00032_00806-1',
  #                                             all_other_ref=rownames(ref_spot)[rownames(ref_spot)!='s_008um_00032_00806-1']
  # ),
  # cols.highlight = c('darkblue','darkred','darkgreen','white'))
  if(!calc.strength){
    nearest_ref_info <- data.frame(query_barcode='',distance='',nearst_barcode='')
    for( i in 1:length(rownames(query_spot))){
      if(i %% 10000 ==1)
      {t.start=Sys.time()}
      res=find_nearest_ref_HD(query_col = query_spot[i,c('x')],
                              query_row = query_spot[i,c('y')],
                              ref_spot =ref_spot ,max.r = max.r)
      nearest_ref_info <- rbind(nearest_ref_info,data.frame(query_barcode=rownames(query_spot)[i],
                                                            distance=res[2],
                                                            nearst_barcode=res[1]))
      
      if(i %% 10000==0){
        t.end=Sys.time()
        print(paste0('proceed ',i,' bins. Use ',round((t.end-t.start),2),' second'))
        
      }
      
      
    }
    nearest_ref_info <- nearest_ref_info[-1,]
  }
  else{
    nearest_ref_info <- data.frame(query_barcode='',distance='',nearst_barcode='',strength=0,strength.sum=0)
    for( i in 1:length(rownames(query_spot))){
      if(i %% 10000 ==1)
      {t.start=Sys.time()}
      res=find_nearest_ref(query_col = query_spot[i,c('x')],
                           query_row = query_spot[i,c('y')],
                           ref_spot =ref_spot ,max.r = max.r,model=model)
      nearest_ref_info <- rbind(nearest_ref_info,data.frame(query_barcode=rownames(query_spot)[i],
                                                            distance=res[2],strength=res[3],strength.sum=res[4],
                                                            nearst_barcode=res[1]))
      if(i %% 10000==0){
        t.end=Sys.time()
        print(paste0('proceed ',i,' bins. Use ',round((t.end-t.start),2),' second'))
        
      }
      
    }
    nearest_ref_info$strength <- as.numeric(nearest_ref_info$strength)
    nearest_ref_info$strength.sum <- as.numeric(nearest_ref_info$strength.sum)
    nearest_ref_info <- nearest_ref_info[-1,]
  }
  
  
  return(nearest_ref_info)
  
}



###PlotNearDis 绘制Spot-最近距离
#seurat.obj 输入对象，只可以包含一张图。
#nearest_ref_info,是CalcNearDis的输出结果
#color是颜色向量，默认darkblue=yellow
#max.dis，最远距离颜色阈值，默认20(spot)
#shape是圆形21（默认），HD可以改为22方块



PlotNearDis <- function(seurat.obj,nearest_ref_info,color=NULL,shape=21,max.dis=20,image.alpha = 0,
                        pt.size.factor = 4e3){
  library(ggpubr)
  st.p3 <- seurat.obj
  st.p3$distance=0
  st.p3@meta.data[nearest_ref_info$query_barcode,'distance']=as.numeric(nearest_ref_info$distance)
  st.p3$distance <- as.numeric(st.p3$distance)
  
  if(is.null(color))
    color=c('darkblue','yellow')
  
  
  p4 <- SpatialFeaturePlot(subset(st.p3,distance>0),features = 'distance',max.cutoff = max.dis,shape=shape,
                           image.alpha = image.alpha,pt.size.factor = pt.size.factor)+
    scale_fill_gradient(low = color[1],high = color[2])
  return(p4)
  
}



###PlotStrengthDis 绘制Spot-最近距离
#seurat.obj 输入对象，只可以包含一张图。
#nearest_ref_info,是CalcNearDis的输出结果,且calc.strength=T
#color是颜色向量，默认darkblue=yellow
#



PlotStrengthDis <- function(seurat.obj,nearest_ref_info,color=NULL,image.alpha = 0,pt.size.factor = 4e3){
  st.p3 <- seurat.obj
  st.p3$distance=0
  st.p3@meta.data[nearest_ref_info$query_barcode,'distance']=as.numeric(nearest_ref_info$distance)
  st.p3$distance <- as.numeric(st.p3$distance)
  
  st.p3$strength.sum=0
  st.p3@meta.data[nearest_ref_info$query_barcode,'strength.sum']=as.numeric(nearest_ref_info$strength.sum)
  #st.p3$strength.sum <- as.numeric(st.p3$strength.sum)
  
  if(is.null(color))
    color=c('lightgrey','darkgreen')
  
  
  p4 <- SpatialFeaturePlot(subset(st.p3,distance>0),features = 'strength.sum',
                           image.alpha = image.alpha,pt.size.factor = pt.size.factor)+
    scale_fill_gradient(low = color[1],high = color[2])
  return(p4)
  
}



###PlotStrengthExpr 绘制Spot-强度、表达量
#seurat.obj 输入对象，只可以包含一张图。
#nearest_ref_info,是CalcNearDis的输出结果,且calc.strength=T
#color是颜色向量，默认darkblue=yellow
#gene是表达量
#



PlotStrengthExpr <- function(seurat.obj,nearest_ref_info,color=NULL,image.alpha = 0,pt.size.factor = 4e3,layer='data',assay='SCT',gene){
  st.p3 <- seurat.obj
  st.p3$distance=0
  st.p3@meta.data[nearest_ref_info$query_barcode,'distance']=as.numeric(nearest_ref_info$distance)
  st.p3$distance <- as.numeric(st.p3$distance)
  
  st.p3$strength.sum=0
  st.p3@meta.data[nearest_ref_info$query_barcode,'strength.sum']=as.numeric(nearest_ref_info$strength.sum)
  #st.p3$strength.sum=1-st.p3$strength.sum/max(st.p3$strength.sum)
  
  #st.p3$strength.sum <- as.numeric(st.p3$strength.sum)
  
  
  
  
  if(is.null(color))
    color='darkgreen'
  
  
  df <- st.p3@meta.data
  df$target <- GetAssayData(st.p3,assay = assay,layer = layer)[gene,]
  p6 <- ggpubr::ggscatter(df,x='strength.sum',y='target',conf.int = T,color =color ,
                          xlab = 'Strength.sum',ylab = paste0('Relative Expression of ',gene),
                          cor.coef = T,add = 'reg.line',size = 0.75,alpha=0.2,
                          add.params = list(color='darkred'))
  return(p6)
  
  
  
}





###PlotDisExpr 绘制距离与表达量
#log.trans 是指距离是否需要log校正，默认F
#ref_col是距离ref的标签列名
#nearest_ref_info是前面CalcNearDis的结果
#filt_far 删去far距离（不含）以上的点，默认NULL不删除
#filt_zero 删去核心点，仅HD
PlotDisExpr <- function(seurat.obj,nearest_ref_info,ref_col,layer='data',assay='SCT',gene,log.trans=F,filt_zero=F,filt_far=NULL){
  
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
  if(!is.null(filt_far)){
    df <- df[df$distance<=filt_far,]
  }
  if(filt_zero){
    df <- df[df$distance>0,]
  }
  
  
  p6 <- ggpubr::ggscatter(df,x='distance',y='target',conf.int = T,color =ref_col ,
                          xlab = 'Distance',ylab = paste0('Relative Expression of ',gene),
                          cor.coef = T,add = 'reg.line',size = 0.75,alpha=0.2,
                          add.params = list(color='darkred'))
  return(p6)
}


###PlotDisProp 绘制距离与比例
#log.trans 是指距离是否需要log校正，默认F
#ref_col是距离ref的标签列名
#nearest_ref_info是前面CalcNearDis的结果
#show.gene是是否需要用点表示表达量，默认NULL,否则请传入基因名

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






