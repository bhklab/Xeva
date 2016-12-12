
if(1==2){

library(ComplexHeatmap)

mat = read.table(textConnection(
    ",s1,s2,s3
g1,snv;indel,snv,indel
g2,,snv;indel,snv
g3,snv,,indel;snv"), row.names = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
mat = as.matrix(mat)

mat
mat[2,1] = "indel;snv;snv;snv"
col = c(snv = "red", indel = "blue", mut="green")

#alter_fun = list(
#  snv = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["snv"], col = NA)),
#  indel = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["indel"], col = NA))
#)

alter_fun <-function(x, y, w, h, v)
  {
  vt = v[v==TRUE]
  lv = length(vt)
  #if(length(v[v==TRUE])==2){ grid.rect(x, y, w, h, gp = gpar(fill = "blue", col = NA))}
  grid.rect(x, y, w, h, gp = gpar(fill = "grey", col = NA))
  if(lv==0){return(NULL)}
  if(lv>0 & lv<4)
  {
    grid.rect(x, (y - h*0.5 + 1:lv/lv*h), w*0.95, h*(1/lv)*0.95,
              gp = gpar(fill = col[names(vt)], col = NA), just = "top")

  }
  if(lv==4){ grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "yellow", col = NA)) }

  }


oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],row_order=NULL, column_order= NULL,
          alter_fun = alter_fun, col = col)

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col,
          row_order=NULL, column_order= NULL)




################=============================================
mat = read.table(paste0(system.file("extdata", package = "ComplexHeatmap"),
                        "/tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"),
                 header = TRUE,stringsAsFactors=FALSE, sep = "\t")
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]

mat[1,1:3] = c("HOMDEL;", "AMP;", "MUT;")
mat[2,1:3] = c("HOMDEL;AMP;", "HOMDEL;MUT;", "")
mat[3,1:3] = c("AMP;MUT", "", "AMP;HOMDEL;MUT;")
mat = mat[1:5,1:5]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)

col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col,
          row_order=NULL, column_order= NULL,
          column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"),
                                      labels = c("Amplification", "Deep deletion", "Mutation")))



oncoPrintX(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col,
          row_order=NULL, column_order= NULL,
          column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"),
                                      labels = c("Amplification", "Deep deletion", "Mutation")))


################=============================================







##------creat barplot matrix --------
computeColumnStat <- function(dfx, idMap, groupBy)
{
  unqID = unique(idMap[,groupBy])
  rtx = data.frame(matrix(NA, nrow = length(unqID), ncol = 4))
  colnames(rtx) = c("CR", "PR", "SD", "PD")
  rownames(rtx) = unqID

  for(I in unqID)
  {
    Ix = subset(idMap, idMap[ , groupBy]==I)
    modx = dfx[, Ix$model.id]
    modv = as.vector(as.matrix(modx))
    modv = modv[!is.na(modv)]
    rtx[I, ]= sapply(c("CR", "PR", "SD", "PD"), function(x) sum(modv==x) )
  }
  return(rtx)
}

plotmRECIST <- function(df)
{
  mRcol = c("CR" = "blue", "PR" = "green", "SD" = "yellow", "PD"="red")

  library(ComplexHeatmap)
  groupBy = "biobase.id"

  data(pdxe)
  #setmRECIST(pdxe)<- setmRECIST(pdxe)
  df = getmRECIST(pdxe)
  df = df[1:500,]
  #df$model.mR = paste(df$model.id, df$mRECIST, sep=" = ")

  dfx = .castDataFram(df, row.var="drug.join.name", col.var = "model.id", value="mRECIST")

  idMap = df[, c(groupBy, "model.id")]
  idMap = BBmisc::sortByCol(idMap, c(groupBy, "model.id") )
  dfx = dfx[, idMap$model.id]

  colBarpltMat = computeColumnStat(dfx, idMap, groupBy)


  column_ha = HeatmapAnnotation(colBarpltMat = anno_barplot(colBarpltMat, axis = TRUE,
                                                    gp = gpar(fill = mRcol),
                                                    width = unit(2, "cm")))


  #Heatmap(mat, top_annotation = column_ha,
  #        top_annotation_height = unit(2, "cm"), km = 2) #+ row_ha
  oncoPrint(dfx, get_type = function(x) strsplit(x, ";")[[1]],
            alter_fun = alter_fun, col = col,
            column_title = "OncoPrint",
            heatmap_legend_param = list(title = "Alternations", at = c("CR", "PR", "SD", "PD"),
                                        labels = c("CR", "PR", "SD", "PD")),
            top_annotation = column_ha)


  ##-------------------------------------
  ##------creat barplot matrix


  foo1 = matrix(abs(rnorm(30)), ncol = 3)
  #foo1[1, ] = -foo1[1, ]
  column_ha = HeatmapAnnotation(foo1 = anno_barplot(foo1, axis = TRUE,
                                                    gp = gpar(fill = c("red", "blue", "green")),
                                                    width = unit(2, "cm")))
  #foo2 = matrix(abs(rnorm(24)), ncol = 2)
  #row_ha = rowAnnotation(foo2 = row_anno_barplot(foo2, axis = TRUE, axis_side = "top",
  #                                               gp = gpar(fill = c("red", "blue"))), width = unit(2, "cm"))

  Heatmap(mat, top_annotation = column_ha, top_annotation_height = unit(2, "cm"), km = 2) #+ row_ha



  ##-----------------------------------

  library(circlize)

  RColorBrewer::brewer.pal(10)

  annoCol = list()
  annoCol[[groupBy]] = colorRamp2(c(0, 20), c("white", "red"))


  df = data.frame(type = c(rep("a", 5), rep("b", 5)))
  ha = HeatmapAnnotation(df = df, col = list(type = c("a" =  "red", "b" = "blue")))
  draw(ha, 1:100)

  ha = HeatmapAnnotation(df = modelMap, col = annoCol)
  #draw(ha, 1:100)



  drgNames = c("HDM201", "BYL719 + LJM716", "BYL719", "CLR457", "untreated", "BKM120", "binimetinib", "LEE011")
  dfx = dfx[drgNames, ]
  dfx = dfx[, colSums(is.na(dfx)) != nrow(dfx)]

  dfx[is.na(dfx)] = ""
  dfx = as.matrix(dfx)

  #dfx = dfx[1:10,1:500]

  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    CR = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
    },
    PR = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "green", col = NA))
    },
    SD = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "yellow", col = NA))
    },
    PD = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
    }
  )


  col = c("CR" = "blue", "PR" = "green", "SD" = "yellow", "PD"="red")

  oncoPrint(dfx, get_type = function(x) strsplit(x, ";")[[1]],
            alter_fun = alter_fun, col = col,
            column_title = "OncoPrint",
            heatmap_legend_param = list(title = "Alternations", at = c("CR", "PR", "SD", "PD"),
                                        labels = c("CR", "PR", "SD", "PD")),
            top_annotation = hpx)










  ##----------------------------------------
  head(dfx)

  mRColPalette = c("CR" = "blue", "PR" = "green", "SD" = "yellow", "PD"="red")

  getColorForCell <- function(mod.mr)
  {
    mod.mr = c(mod.mr)
    mR = sapply(mod.mr, function(x){ strsplit(mod.mr, " = ")[[1]][2] })
    if(length(mR)>0 & all(is.na(mR))==FALSE)
    { return(mRColPalette[mR]) }
    #mR = strsplit(mod.mr, " = ")[[1]][2]
    #if(!is.na(mR)){ return(mR.col[mR]) }
    return("black")
  }

  alter_fun <-function(x, y, w, h, v)
  {
    vt = v[v==TRUE]
    lv = length(vt)
    grid.rect(x, y, w, h, gp = gpar(fill = "grey", col = NA))
    if(lv==0){return(NULL)}
    if(lv>0)
    {
      grid.rect(x, (y - h*0.5 + 1:lv/lv*h), w*0.95, h*(1/lv)*0.95,
                gp = gpar(fill = getColorForCell(names(vt)), col = NA), just = "top")

    }

  }


  oncoPrint(dfx, get_type = function(x) strsplit(x, ";")[[1]],row_order=NULL, column_order= NULL,
            alter_fun = alter_fun, col = mRColPalette)



  getTextforCell <- function(V)
  {
    #if(V==0) { return("") }
    return(sprintf("%1.1f", V))
  }



  heatX = Heatmap(mat, name = "Drug & Models", col=c("blue", "green", "yellow", "red"),
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  #column_dend_height = unit(5, "cm"),
                  show_row_dend = FALSE,
                  show_row_names = TRUE, row_names_side = "left",
                  show_column_names = TRUE, column_names_side = "top",
                  rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                  show_heatmap_legend = FALSE,
                  #top_annotation = ha1
                  cell_fun = function(j, i, x, y, width, height, col)
                  {
                    #grid.text(getTextforCell(mat[i, j]), x, y, gp = gpar(fontsize = 12,col = "black"))
                    grid.rect(x, y, width*0.9, height*0.4, gp = gpar(fill = "black", col = NA))
                  }  )



  draw(heatX, padding = unit(c(10, 10, 10, 10), "mm"))



}




plotHeatmapX <- function()
{

  library(ComplexHeatmap)
  getCol<- function(mr)
  {
    mRColPalette = c("CR" = "blue", "PR" = "green", "SD" = "yellow", "PD"="red")
    strsplit(mr, ";")[[1]]
  }



  mat = matrix(NA, nrow = 5, ncol = 4)
  mat[1,] = c("CR", "PR", "SD", "PD")
  mat[2,] = c("CR;PR",NA, "SD;PD", "PD")
  mat[3,] = c("CR;PR;CR","SD;PD;PD;SD", "PD;PD", "PR")
  mat[4,] = c("PR","CR", "SD", "PR")
  mat[5,] = c("CR","SD;SD", "PR;PD", "SD")

  nameSpc = unique(as.vector(mat))
  backgroundCol = rep("gray", length(nameSpc))


  getSqrPoint <- function(x,y,l)
  {
    X = c(x-l, x-l, x+l, x+l)
    Y = c(y-l, y+l, y+l, y-l)
    return(list(X=X,Y=Y))
  }

  cell_fun_old = function(j, i, x, y, width, height, fill)
  {
    #grid.polygon(x = rep(x[c("left", "right", "right", "left")], 2),
    #           y =  y[rep(c("bottom", "top", "bottom", "mid"), each = 2)],
    #           id = rep(1:2, each = 4),
    #           gp = gpar(fill = c(NA, "blue")))


    grid.polygon(#x=c(0, 0, 1, 1), y=c(0, 1, 1, 0),
      x=c(0.05, 0.05, 0.95, 0.95), y=c(0.05, 0.95, 0.95, 0.05),
               #id = rep(1:2, each = 4),
               gp = gpar(fill = "blue"))

    #z=list(j, i, x, y, width, height, fill)
    #saveRDS(z, file="tmp.Rda")
    #grid.rect(x, y, width*0.9, height*0.4, gp = gpar(fill = "black", col = NA))
    #grid.circle(x, y, r, gp = gpar(fill = "red", col = NA))

    #grid.polygon(#x=c(0, 0, 1, 1), y=c(0, 1, 1, 0),
    #  x=c(0.05, 0.05, 0.95, 0.95), y=c(0.05, 0.95, 0.95, 0.05),
    #  #id = rep(1:2, each = 4),
    #  gp = gpar(fill = "blue"))

  }

  cell_fun_old2 = function(j,i, x, y, w, h, v) {
    z=list(j,i, x, y, w, h, v)
    saveRDS(z, file="tmp.Rda")

    #w = convertWidth(w, "cm")*0.9
    #h = convertHeight(h, "cm")*0.9
    l = min(unit.c(w, h))*0.5

    #grid.circle(x, y, l*0.5, gp = gpar(fill = "red", col = NA))
    #sc = getSqrPoint(x,y,l)
    #w = z[[5]]; h=z[[6]]
    #w = convertWidth(w, "cm")*0.9
    #h = convertHeight(h, "cm")*0.9
    #l = min(unit.c(w, h))

    #x=z[[3]]; y=z[[4]];
    #x=convertX(z[[3]], "cm", valueOnly = FALSE)
    #y=convertX(z[[4]], "cm", valueOnly = FALSE)
    a=convertX(x-l, "npc", valueOnly = FALSE)
    b=convertX(x+l, "npc", valueOnly = FALSE)
    c=convertX(y-l, "npc", valueOnly = FALSE)
    d=convertX(y+l, "npc", valueOnly = FALSE)
    #b=x+l
    #c=y-l ; d=y+l
    #grid.polygon(x = c(x-l, x-l, x+l, x+l),
    #             y = c(y-l, y+l, y+l, y-l),
    #             gp = gpar(fill = "blue"))

    grid.polygon(x = c(a, a, b, b),
                 y = c(c, d, d, c),
                 gp = gpar(fill = "blue"))

    #grid.circle(c, d, l*0.5, gp = gpar(fill = "red", col = NA))



  }


z=list() ; saveRDS(z, file="tmp.Rda")
  cell_funX = function(j,i, x, y, w, h, v, value) {
    factR = 0.95
    wr=w*0.5*factR; hr=h*0.5*factR
    a=x-wr; b=x+wr
    c=y-hr; d=y+hr
    a=convertX(a, "npc", valueOnly = FALSE)
    b=convertX(b, "npc", valueOnly = FALSE)
    c=convertX(c, "npc", valueOnly = FALSE)
    d=convertX(d, "npc", valueOnly = FALSE)

    z=readRDS("tmp.Rda")
    z[[length(z)+1]] = list(j,i, x, y, w, h, a,b,c,d, value)
    saveRDS(z, file="tmp.Rda")

    grid.polygon(x = c(a, a, b, b),
                 y = c(c, d, d, c),
                 gp = gpar(fill = "blue"))
  }

  Heatmap(mat, name = "Drug & Models", col=backgroundCol,
          cluster_rows = FALSE, cluster_columns = FALSE,
          #column_dend_height = unit(5, "cm"),
          show_row_dend = FALSE,
          show_row_names = TRUE, row_names_side = "left",
          show_column_names = TRUE, column_names_side = "top",
          rect_gp = gpar(col = "white", lty = 1, lwd = 1),
          show_heatmap_legend = FALSE,
          cell_fun =function(j, i, x, y, width, height, fill)
          {         cell_funX(j, i, x, y, width, height, fill, mat[i,j])}
  )


          #cell_fun = function(j, i, x, y, width, height, fill)
          #{
          #  #grid.text(getTextforCell(mat[i, j]), x, y, gp = gpar(fontsize = 12,col = "black"))
          #  grid.rect(x, y, width*0.9, height*0.4, gp = gpar(fill = "black", col = NA))
          #
          #}
          #)







ht = Heatmap(pheudo, col = col, rect_gp = gpar(type = "none"),
               cluster_rows = FALSE, cluster_columns = FALSE, row_order = row_order, column_order = column_order,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 z = arr[i, j, ]
                 names(z) = dimnames(arr)[[3]]
                 af(x, y, width, height, z)
               },
               heatmap_legend_param = heatmap_legend_param)




}


}
