library(Seurat)
library(stringr)
set.seed(1113)
# image_hires <- Read10X_Image(image.name = 'tissue_hires_image.png',assay='Spatial',
#                              image.dir = '~/yifanfu/Rserver/24sweet/Matrix/sp/forR/treat/outs/binned_outputs/square_050um/spatial/')



####Visium Demo
st1.t1 <- readRDS('~/yifanfu/Arx_data/Public/24NatST/st1_demo_stGrads.rds')
st1.t1 <- UpdateSeuratObject(st1.t1)
st.p1.ann <- read.csv('/public/home/qiujiangdong/yifanfu/ijupyter/PDAC/HT/p1.cell2location.csv',row.names = 1)
#factor=4e3
st1.t1@images$IU_PDA_T1@scale.factors
st1.t1@meta.data <- cbind(st1.t1@meta.data,st.p1.ann)
colnames(st1.t1@meta.data) <- gsub('q05cell_abundance_w_sf_','',colnames(st1.t1@meta.data))


st1.t1 <- FindPrimSpot(
  seurat.obj   = st1.t1,
  celltype_prop = 'MDSC',
  celltype      = 'MDSC_spot',
  probs         = 0.9
)
table(st1.t1$MDSC_spot)

SpatialDimPlot(st1.t1,group.by = 'MDSC_spot',
               pt.size.factor = 4e3)




df.res <- CalcNearDis(
  st1.t1,
  celltype      = 'MDSC_spot',
  pheno_choose  = NULL,       # NOT RECOMMENDED now
  calc.strength = T,
  model         = 'Linear',   # 'Linear' | 'Log' | 'Exp'
  max.r         = 10          # for Visium you can use 4; here 10 shows a wider range
)

# df.res is a data.frame summarizing nearest-ring relationships (and strength if requested)
head(df.res)

p1 <- PlotNearDis(
  st1.t1,
  nearest_ref_info = df.res,
  color            = c('darkred','gray'),
  shape            = 21,
  max.dis          = 6,
  image.alpha      = 0.5,
  pt.size.factor   = 4e3
)+ggtitle('Distance')

p1 <- PlotStrengthDis(
  st1.t1,
  df.res,
  img.use = 'hires',
  color        = c('gray','darkred'),
  image.alpha  = 0.5,
  pt.size.factor = 4e3
)+ggtitle('Strength_Linear model')
p2 <- PlotStrengthDis(
  st1.t1,
  df.res,
  color        = c('gray','darkred'),
  image.alpha  = 0.5,
  pt.size.factor = 4e3
)+ggtitle('Strength_Log model')
p3 <- PlotStrengthDis(
  st1.t1,
  df.res,
  color        = c('gray','darkred'),
  image.alpha  = 0.5,
  pt.size.factor = 4e3
)+ggtitle('Strength_Exp model')



p1|p2|p3

PlotStrengthExpr(
  st1.t1,
  df.res,
  image.alpha    = 0.5,
  pt.size.factor = 4e3,
  layer          = 'data',
  assay          = 'SCT',
  gene           = 'COL1A1'
)



PlotDisExpr(
  st1.t1,
  df.res,
  ref_col   = c('darkred'),
  gene      = 'COL1A1',
  log.trans = FALSE
)


PlotDisProp(
  st1.t1,
  df.res,
  alpha = 0.5,
  ref_col    = c('darkred'),
  prop_cell  = 'PSC'     # name of the proportion column for the target cell type
)



#############

st.zl <- Load10X_Spatial('~/yifanfu/Rserver/24sweet/Matrix/sp/forR/treat/outs/binned_outputs/square_050um/' )
st.zl <- subset(st.zl,nCount_Spatial>0)

SpatialDimPlot(st.zl,shape=22,crop = F)
st.zl <- SCTransform(st.zl,assay = 'Spatial')

#read in cell2loc results
{
  comp <- read.csv('~/yifanfu/Rserver/24sweet/ZL2_component.csv',row.names = 1)
  colnames(comp) <- gsub('q05cell_abundance_w_sf_','',colnames(comp))

  integrated_compositions <- comp
  if (max(integrated_compositions) <= 1){integrated_compositions = integrated_compositions * 20}
  #这一行用来处理非cell2loc的数据
  colnames(integrated_compositions)[colnames(integrated_compositions)=='KC_Diff']='dKC'
  colnames(integrated_compositions)[colnames(integrated_compositions)=='CD8Tcell']='CD8T'
  colnames(integrated_compositions)[colnames(integrated_compositions)=='KC_UnDiff']='udKC'
  colnames(integrated_compositions)[colnames(integrated_compositions)=='Vellus.hair.follicle']='VHF'
  colnames(integrated_compositions)[colnames(integrated_compositions)=='LangerhansCell']='Langerhans'
  colnames(integrated_compositions)[colnames(integrated_compositions)=='MastCell']='Mast'
  colnames(integrated_compositions)[colnames(integrated_compositions)=='Endothelial']='Endo'
  colnames(integrated_compositions)[colnames(integrated_compositions)=='CD4Tcell']='CD4T'

}
st.zl <- AddMetaData(st.zl,integrated_compositions)

colnames(st.zl) %in% rownames(integrated_compositions)


SpatialFeaturePlot(st.zl,features = 'dKC',alpha=c(0.4,1),shape=22,crop = F)

#Step1
##Find the primary infiltration spots of dKC
#celltype_prop is given proportion of specific celltype
#celltype  store classification
#probs is cutoff q90

st.zl <- FindPrimSpot(seurat.obj = st.zl,celltype_prop = 'dKC',
                       celltype = 'dKC_spot',probs = 0.95)
table(st.zl$dKC_spot)
# dKC_spot   Others
# 182     3443


st.zl <- FindPrimSpot(seurat.obj = st.zl,celltype_prop = 'Neu_4',
                      celltype = 'Neu_4_spot',probs = 0.95)
table(st.zl$Neu_4_spot)
# Neu_4_spot     Others
# 182       3443


SpatialDimPlot(st.zl,group.by = 'dKC_spot',alpha=c(0.8),shape=22,crop = F)
SpatialDimPlot(st.zl,group.by = 'Neu_4_spot',alpha=c(0.8),shape=22,crop = F)



#Step2
##Calculate Nearest spots,return a dataframe
#celltype is classification from FindPrimSpot()
# pheno_choose is a list of spots, NOT RECOMMENDED now.
#model is "Linear", "Log", or "Exp"
#max.r is limitation of SPOTS, def=10, for Visium can use 4 (55um * 4 ~ 200um )
df.res <- CalcNearDis_HD(st.zl,celltype = 'dKC_spot',pheno_choose = 'Neu_4_spot',resolution = 50,
                      calc.strength = F,model = 'Linear',max.r = 200)





#####HD demo######

set.seed(1113)
library(Seurat)
library(ggplot2)
library(ggsci)
library(stringr)

library(stGrads)
setwd("~/yifanfu/Rserver/BH")
options(future.globals.maxSize = 10*1024^3)

# ---- Load Visium HD Seurat object ----
st.hd <- readRDS("HD_demo_1.rds")
SpatialDimPlot(
  st.hd,
  group.by       = "first_type",
  image.scale    = "hires",
  crop           = TRUE,
  shape          = 22,  # square spots for HD
  pt.size.factor = 80 * st.hd@images[["slice1"]]@scale.factors[["hires"]],
  alpha          = 0.75
) +
  scale_fill_manual(values = c("#FFF68F","#CDBE70","#CD2626","#008B8B","#8B3A3A"))



dis_mtx.Duct2_myCAF <- CalcNearDis_HD(
  st.hd,
  celltype       = "first_type",
  query_celltype = "Duct2",
  pheno_choose   = "myCAF",
  max.r          = 100,
  resolution     = 8
)

dis_mtx.Duct2_Mp <- CalcNearDis_HD(
  st.hd,
  celltype       = "first_type",
  query_celltype = "Duct2",
  pheno_choose   = "Macrophage",
  max.r          = 100,
  resolution     = 8
)


dis_mtx.Duct2_myCAF$ref <- "myCAF"
dis_mtx.Duct2_Mp$ref    <- "Macrophage"
dff <- rbind(dis_mtx.Duct2_myCAF, dis_mtx.Duct2_Mp)
dff$distance <- as.numeric(dff$distance)

ggplot(dff) +
  geom_density(aes(x = distance, colour = ref), linewidth = 1) +
  theme_bw() +
  labs(title = "Duct2 → myCAF / Macrophage", x = "Distance (μm)", y = "Density") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#CD2626","#008B8B"))

stg_mtx.Duct2_myCAF <- CalcNearDis_HD(
  st.hd,
  celltype       = "first_type",
  query_celltype = "Duct2",
  pheno_choose   = "myCAF",
  max.r          = 100,
  resolution     = 8,
  calc.strength  = TRUE,
  model          = "Linear"
)


stg_mtx.Duct2_Mp <- CalcNearDis_HD(
  st.hd,
  celltype       = "first_type",
  query_celltype = "Duct2",
  pheno_choose   = "Macrophage",
  max.r          = 100,
  resolution     = 8,
  calc.strength  = TRUE,
  model          = "Linear"
)

stg_mtx.Duct2_myCAF$ref <- "myCAF"
stg_mtx.Duct2_Mp$ref    <- "Macrophage"
dff.stg <- rbind(stg_mtx.Duct2_myCAF, stg_mtx.Duct2_Mp)

# Normalize summed strength for visualization
dff.stg$strength.sum <- dff.stg$strength.sum / max(dff.stg$strength.sum, na.rm = TRUE)

ggplot(dff.stg) +
  geom_density(aes(x = strength.sum, colour = ref), linewidth = 2) +
  theme_bw() +
  labs(title = "Duct2 → myCAF / Macrophage (r = 100 μm)", x = "Relative Strength", y = "Density") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#CD2626","#008B8B"))
st.hd <- UpdateSeuratObject(st.hd)

p1 <- PlotNearDis(
  st.hd,
  nearest_ref_info = dis_mtx.Duct2_Mp,
  max.dis       = 100,
  color         = c("darkred","gray"),
  shape         = 22,
  #img.use       = "hires",
  image.alpha   = 1,
  pt.size.factor = 80 * st.hd@images[["slice1"]]@scale.factors[["hires"]]
) + labs(title = "Duct2 → Macrophage (r = 100 μm)")

p2 <- PlotStrengthDis(
  st.hd,
  nearest_ref_info = stg_mtx.Duct2_Mp,
  #max.dis       = 100,
  color         = c("gray","darkgreen"),
  shape         = 22,
  img.use       = "hires",
  image.alpha   = 1,
  pt.size.factor = 80 * st.hd@images[["slice1"]]@scale.factors[["hires"]]
) + labs(title = "Duct2 → Macrophage (r = 100 μm)")



p1|p2



st.hd <- NormalizeData(st.hd)

DRGs <- CalcDRG(
  st.hd,
  nearest_ref_info = dis_mtx.Duct2_Mp,
  layer      = "data",
  assay      = "Spatial",
  pthresh    = 1,
  log.trans  = FALSE,
  filt_root  = TRUE,
  filt_far   = 100
)
DRGs <- lm.df
head(DRGs)
# rho     pvalue
# NOC2L   -0.09532062 0.46493640
# KLHL17  -0.05884603 0.65236098
# PLEKHN1  0.08831333 0.49853035
# HES4    -0.23190598 0.07211730
# ISG15   -0.31018989 0.01497929
# AGRN    -0.12869089 0.32294101
PlotDRGs(DRGs)


p1 <- PlotDisExpr(
  seurat.obj       = st.hd,
  nearest_ref_info = dis_mtx.Duct2_Mp,
  ref_col   = "distance",
  assay     = "Spatial",
  gene      = "SPARC",
  col       = "darkred",
  log.trans = FALSE,
  filt_root = TRUE,
  filt_far  = 100
)

p2 <- PlotDisExpr(
  seurat.obj       = st.hd,
  nearest_ref_info = dis_mtx.Duct2_Mp,
  ref_col   = "distance",
  assay     = "Spatial",
  gene      = "APOL1",
  col       = "darkgreen",
  log.trans = FALSE,
  filt_root = TRUE,
  filt_far  = 100
)


p1|p2
sessionInfo()
