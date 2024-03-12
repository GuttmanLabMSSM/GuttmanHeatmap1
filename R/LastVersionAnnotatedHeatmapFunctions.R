
doTab<-function(fdrsx,coefsx){
  require(stringr)
  require(weights)
  
  p.cuts.stars=c(0.001,0.01, 0.05, 0.1,1)
  p.symbols.stars= c('***',"**", "*", "+","  ")
  
  coefsx[coefsx==0]<-1e-10
  
  fdrsx[] <- sapply(fdrsx,starmaker, p.levels =p.cuts.stars, symbols = p.symbols.stars)
  
  transformfch <- function(lgfch) {
    fch <- sign(lgfch) * (2^abs(lgfch))
    return(fch)
  }
  
  
  
  padding<-function(x,n){
    #print(x)
    m=ceiling(log10(abs(x)+.001))
    x<-as.character(x)
    
    initspace=n-m
     if (substr(x,1,1)!='-'){
       initspace<-initspace+1
     }
    if (!str_detect(x,'\\.')){
      x<-paste(x,'.',sep='')
    }
    while (nchar(x)<n-initspace+4){
      x<-paste(x,'0',sep='')
    }
    initstr=paste0("<span style='color:white'>",str_c(rep('\\-',initspace),collapse=''))
    x<-paste(initstr,"</span>",x,sep='')
    return(x)
  }


  coefsx[]<-sapply(coefsx,transformfch)
  coefsx<-round(coefsx,2)
  
  ss<-rownames(coefsx)
  maxss <- max(sapply(ss, nchar))
  adspaceText <- function(x, n) {paste0(x,"<span style='font-size:12pt;color:black'>",' ',str_c(rep("\\-",(n-nchar(x))),collapse = ''),"</span>","\\:")}
  
  ss <- sapply(ss, adspaceText, maxss+1)
  
  maxdigits=max(ceiling(log10(max(abs(coefsx)))),3)
  Tab <- data.frame(Symbol = ss)
  for (i in c(1:ncol(coefsx))) {
    Tab[colnames(coefsx)[i]] <- 
      paste(sapply(coefsx[, i],padding,maxdigits), 
            fdrsx[,i], sep = "")
  }
  rownames(Tab)= ss
  colnames(Tab)= c('Symbol', colnames(coefsx)) 
  
  return(Tab)
}


doAnnotatedHeatmap <-function(
    grafname, # grafname
    mat, # matrix of expression
    annot.col, # annotation data
    colsplitx=NULL, # splitting column
    column_ha = NULL,
    BigTab = NULL,
    wdt=0, # width
    hgt=0, # height
    setname='', # description of the gene set
    column_names_max_height = unit(8, "cm"),
    BigTabcols=NULL, #contrasts from bigtab to include in annotations. Defaults to all
    colnamesubs=c()

    ){
  if (!is.null(BigTabcols)){
  allBigTabcols<-c()
  for (i in BigTabcols){

    
    for (j in colnames(BigTab)){
      if (strsplit(j,'_')[[1]][2]==i){
        allBigTabcols<-append(allBigTabcols,j)
      }
    }
  
    }
  BigTab<-BigTab[allBigTabcols]
  }
  
  cfxfunc<-function(x){substr(x,1,5)=='lgFCH'}
  fdxfunc<-function(x){substr(x,1,5)=='pvals'}
  cfx=BigTab[,colnames(BigTab)[unlist(lapply(colnames(BigTab),cfxfunc))]]
  fdx=BigTab[,colnames(BigTab)[unlist(lapply(colnames(BigTab),fdxfunc))]]
  
  ncontrasts = dim(cfx)[2]
  
  if (wdt==0){wdt=dim(mat)[2]*2+dim(BigTab)[2]/10}
  if (hgt==0){hgt=dim(mat)[1]/4+2}
  
  mat = mat[rownames(BigTab),rownames(annot.col)]
  
  scale_rows<-function (x) 
  {
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m)/s)
  }
  
  mat=scale_rows(mat)
  
  
  
  Tab = doTab(fdrsx = fdx , coefsx = cfx)
  
  
  require(circlize)
  require(ComplexHeatmap)
  col_fun = colorRamp2(c(-2,0,2), c("blue", "white", "red"))
  
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
  
  
  ha_names = rowAnnotation(Gene = anno_text(gt_render(Tab[,1], align_widths = TRUE),location=0,just='left' ,
                                       gp = gpar(fill ='white', border= 'white', col = "black", fontface = 2,  fontsize = 12,
                                                 fontfamily = "mono")))
  htlist=ht+ha_names

  for (i in 1:ncontrasts){
    
    thisha=rowAnnotation( x=anno_text( gt_render((Tab[,i+1]), align_widths = FALSE),
                                               gp = gpar(box_col = "white", box_lwd =1, fontface=2, fontsize = 12,
                                                         fontfamily = "mono"),
                                               just = "left",show_name = FALSE))
    names(thisha)<-paste('Tab',i,sep='')
    htlist=htlist+thisha
  }

  
  # write the file to a pdf file
  pdf(file = grafname,width=wdt,height=hgt)
  draw(htlist,padding = unit(c(8, 8, 12, 12), "mm") ,  column_title = setname ,column_title_side = 'bottom')
  
  if (length(colnamesubs)>0){
    colnamesubmat<-matrix(colnamesubs,ncol=2,byrow = TRUE)
    for (i in 1:dim(colnamesubmat)[1]){
      for (j in 1:dim(Tab)[2]){
        colnames(Tab)[j]<-gsub(colnamesubmat[i,1],colnamesubmat[i,2],colnames(Tab)[j])
      }
    }
  }
  
  
  for (i in 1:ncontrasts){
    nm=paste('Tab',i,sep='')
    decorate_annotation(nm, {
      grid.text(gsub('\\.','\n',substring(colnames(Tab)[i+1],7)), y = unit(1, "npc") + unit(2, "mm"), just = "bottom",hjust =0,rot=0,
                gp = gpar(box_col = "white", box_lwd =1, fontface=2,fontsize=8,fontfamily='mono'))
    })
  }
  
  
  dev.off()
  
}
