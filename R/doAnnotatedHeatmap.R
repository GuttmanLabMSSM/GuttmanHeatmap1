doAnnotatedHeatmap <-function(
  grafname, # grafname
  wdt, # width
  hgt, # height
  setname, # description of the gene set
  mat, # matrix of expression
  column_names_max_height = unit(8, "cm"),
  annot.col, # annotation data
  colsplitx=NULL, # splitting column
  column_ha = NULL,
  cfx = NULL, # lgfch of contrasts of interest
  fdx  = NULL, # fdr or pvalue for the contrasts of interest
  tabNames = NULL # example of 3 contrasts
){
  
  ncontrasts = dim(cfx)[2]
  
  mat = mat[,rownames(annot.col)]
  
  scale_rows<-function (x) 
  {
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m)/s)
  }
  
  mat=scale_rows(mat)
  
  
  source('../Rfunctions//doTab.R')
  Tab = doTab(fdrsx = fdx , coefsx = cfx)
  
  
  require(circlize)
  require(ComplexHeatmap)
  col_fun = colorRamp2(c(-1.5,0,1.5), c("blue", "white", "red"))
  
  row_dend = as.dendrogram(hclust(dist(mat)))
  
  ht = ComplexHeatmap::Heatmap(mat,
                               name = "Z-score", 
                               # row_km = 3,
                               column_split = colsplitx ,
                               column_title_gp = gpar(fontsize =20),
                               col = col_fun,
                               top_annotation = column_ha ,
                               cluster_rows = T,
                               cluster_columns = F,  
                               show_row_names = F,
                               show_column_names = F,
                               column_order = (rownames(annot.col)),
                               use_raster = TRUE)
  
  Tab[,-1]=sapply(Tab[,-1],formatTab)
  
  
  
  ha1 = rowAnnotation(Gene = anno_text(gt_render(Tab[,1], align_widths = TRUE),location=0,just='left' ,
                                       gp = gpar(fill ='white', border= 'white', col = "black", fontface = 2,  fontsize = 12,
                                                 fontfamily = "mono")))
  
  
  ha2 = rowAnnotation( 'Tab2' = anno_text( gt_render((Tab[,2]), align_widths = TRUE),
                                           gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,
                                                     fontfamily = "mono"),
                                           just = "left",show_name = FALSE))
  
  ha3 = rowAnnotation( 'Tab3' = anno_text( gt_render((Tab[,3]), align_widths = TRUE),
                                           gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,
                                                     fontfamily = "mono"),
                                           just = "left",show_name = FALSE))
  
  
  if(ncontrasts >2){
  ha4 = rowAnnotation( 'Tab4' = anno_text( gt_render((Tab[,4]), align_widths = TRUE),
                                           gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,
                                                     fontfamily = "mono"),
                                           just = "left",show_name = FALSE))
  }
  
  if(ncontrasts >3){
  ha5 = rowAnnotation( 'Tab5' = anno_text( gt_render((Tab[,5]), align_widths = TRUE),
                                           gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,
                                                     fontfamily = "mono"),
                                           just = "left",show_name = FALSE))
  }
  
  if(ncontrasts>4){
    ha6 = rowAnnotation( 'Tab6' = anno_text( gt_render((Tab[,6]), align_widths = TRUE),
                                             gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,
                                                       fontfamily = "mono"),
                                             just = "left",show_name = FALSE))
  }
  
  if(ncontrasts>5){
    ha7 = rowAnnotation( 'Tab7' = anno_text( gt_render((Tab[,7]), align_widths = TRUE),
                                             gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,
                                                       fontfamily = "mono"),
                                             just = "left",show_name = FALSE))
  }
  

  
  if(ncontrasts==2){
    htlist = ht + ha1 + ha2 + ha3
  }else if(ncontrasts==3){
    htlist = ht + ha1 + ha2 + ha3 + ha4
  }else if(ncontrasts==4){
    htlist = ht + ha1 + ha2 + ha3 + ha4 + ha5
  }else if(ncontrasts==5){
    htlist = ht + ha1 + ha2 + ha3 + ha4 + ha5 + ha6
  }else if(ncontrasts==6){
    htlist = ht + ha1 + ha2 + ha3 + ha4 + ha5 + ha6 + ha7
  }
  
  
  # write the file to a pdf file
  pdf(file = grafname,width=wdt,height=hgt)
  draw(htlist,padding = unit(c(8, 8, 12, 12), "mm") ,  column_title = setname ,column_title_side = 'bottom')
  
  
  decorate_annotation('Tab2', {
    grid.text(tabNames[1], y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
              gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=8))
  })
  decorate_annotation('Tab3', {
    grid.text(tabNames[2], y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
              gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=8))
  })
  
  if(ncontrasts>2){
  decorate_annotation('Tab4', {
    grid.text(tabNames[3], y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
              gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=8))
  })}
  
  if(ncontrasts>3){
  decorate_annotation('Tab5', {
    grid.text(tabNames[4], y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
              gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=8))
  })}
  
  if(ncontrasts>4){
    decorate_annotation('Tab6', {
      grid.text(tabNames[5], y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
                gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=8))
    })}
  
  if(ncontrasts>5){
    decorate_annotation('Tab7', {
      grid.text(tabNames[6], y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
                gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=8))
    })}
  
  if(ncontrasts>6){
    decorate_annotation('Tab8', {
      grid.text(tabNames[7], y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
                gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=8))
    })}
  
  
  dev.off()
  
}
