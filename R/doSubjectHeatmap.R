doSubjectHeatmap <-function(
  grafname, # grafname
  wdt, # width
  hgt, # height
  setname, # description of the gene set
  matx, # matrix of expression
  column_names_max_height = unit(8, "cm"),
  annot.col, # annotation data
  colsplitx, # splitting column
  column_ha = NULL,
  annot_bottom = NULL,
  clust = FALSE, #cluster columns
  rown=F #show rownames
){
  
  mat = as.matrix(matx)
  
  
  scale_rows<-function (x) 
  {
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m)/s)
  }
  
  mat=scale_rows(mat)
  
  require(circlize)
  require(ComplexHeatmap)
  col_fun = colorRamp2(c(-1.5,0,1.5), c("blue", "white", "red"))
  
  row_dend = as.dendrogram(hclust(dist(mat)))
  
  ht = ComplexHeatmap::Heatmap(mat,
                               name = "Z-score", 
                               # row_km = 3,
                               column_split = colsplitx ,
                               col = col_fun,
                               top_annotation = column_ha,
                               bottom_annotation = annot_bottom,
                               column_title_gp = gpar(fontsize =20),
                               cluster_rows = T,
                               cluster_columns = clust,  
                               show_row_names = rown,
                               show_column_names = F,
                               column_order = (rownames(annot.col)),
                               use_raster = TRUE)
  
  # write the file to a pdf file
  pdf(file = grafname,width=wdt,height=hgt)
  htlist = ht
  draw(htlist,padding = unit(c(8, 8, 12, 12), "mm") ,  column_title = setname ,column_title_side = 'bottom')
  dev.off()
  
}
