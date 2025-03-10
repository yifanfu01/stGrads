library(stGrad)


#Step1
##Find the primary infiltration spots of neutrophils
#celltype_prop is given proportion of specific celltype
#celltype  store classification
#probs is cutoff q90

st1.t1 <- FindPrimSpot(seurat.obj = st1.t1,celltype_prop = 'MDSC',celltype = 'MDSC_spot',probs = 0.9)
table(st1.t1$MDSC_spot)
# MDSC_spot    Others 
# 353      3177 


#Step2 
##Calculate Nearest spots,return a dataframe
#celltype is classification from FindPrimSpot()
# pheno_choose is a list of spots, NOT RECOMMENDED now.
#model is "Linear", "Log", or "Exp"
#max.r is limitation of SPOTS, def=10, for Visium can use 4 (55um * 4 ~ 200um )
df.res <- CalcNearDis(st1.t1,celltype = 'MDSC_spot',pheno_choose = NULL,
                      calc.strength = T,model = 'Linear',max.r = 10)

#Step3
##Visualization-basic
PlotNearDis(st1.t1,nearest_ref_info = df.res,color = c('darkred','gray'),shape = 21,
            max.dis = 6,image.alpha = 0.5,pt.size.factor = 4e3)

#Step3
##Visualization-strength
PlotStrengthDis(st1.t1,df.res,color =c('gray','darkred'),
                image.alpha = 0.5,pt.size.factor = 4e3 )


#Step4
##Correlation of Expr-Strength
PlotStrengthExpr(st1.t1,df.res,image.alpha = 0.5,pt.size.factor = 4e3,
                 layer = 'data',assay = 'SCT',gene = 'COL1A1')


#Step4
##Correalation of Expr-Dis
PlotDisExpr(st1.t1,df.res,ref_col = c('darkred'),gene = 'COL1A1',log.trans = F)



#Step4
##Correalation of Proportion-Dis
PlotDisProp(st1.t1,df.res,ref_col = c('darkred'),prop_cell = 'PSC')



