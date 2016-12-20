
getTestMat = function()
{
  mat = matrix(NA, nrow = 5, ncol = 4)
  mat[1,] = c("CR", "PR", "SD", "PD")
  mat[2,] = c("CR;PR",NA, "SD;PD", "PD")
  mat[3,] = c("CR;PR;CR","SD;PD;PD;SD", "PD;PD", "PR")
  mat[4,] = c("PR","CR", "SD", "PR")
  mat[5,] = c("CR","SD;SD", "PR;PD", "SD")
  rownames(mat)=paste0("Drug", 1:dim(mat)[1])
  colnames(mat)=paste0("Sample", 1:dim(mat)[2])
  mat
}

###-------------------------------------------------
.splitValue <- function(mRx, splitBy=";", sort=TRUE)
{
  if(is.na(mRx)){return(NA)}
  mRy = strsplit(mRx, splitBy)[[1]]
  if(sort==TRUE)
  {mRy = sort(mRy)}
  return(mRy)
}

getCellBoxCordi <- function(x0,x1,y0,y1, N)
{
  XV = rep(c(x0,x1,x1,x0), N)
  sqD = (y1-y0)/N
  Nseq = seq(y0,y1, sqD)
  Yv = unlist(lapply(Nseq, function (i) rep(i,4)))
  Li = length(Yv)-2
  YV = Yv[3:Li]
  XV = c(x0, x1, x1, x0, XV)
  YV = c(y0, y0, y1, y1, YV)
  return(list(x=XV, y=YV))
}

.custom_cell_fun <- function(x, y, w, h, value, colPalette, backgroundCol, splitBy=";", sort=TRUE)
{
  factR = 0.95
  wr=w*0.5*factR;   hr=h*0.5*factR
  x0=x-wr; x1=x+wr; y0=y-hr; y1=y+hr
  x0=convertX(x0, "npc", valueOnly = TRUE)
  x1=convertX(x1, "npc", valueOnly = TRUE)
  y0=convertX(y0, "npc", valueOnly = TRUE)
  y1=convertX(y1, "npc", valueOnly = TRUE)
  vx = .splitValue(value, splitBy=";", sort = sort)
  filCol = unlist(colPalette[vx])
  N=length(vx)
  cordXY = getCellBoxCordi(x0,x1, y0, y1, N)
  cordXY$x = unit(cordXY$x,"npc"); cordXY$y = unit(cordXY$y,"npc")
  grid.polygon(x = cordXY$x, y = cordXY$y,
               id = rep(1:(N+1), each = 4),
               gp = gpar(fill = c(NA, filCol), col = NA))
}



.calculatRowColStat <- function(mat, splitBy=";", scaleRow=TRUE, scaleCol=TRUE)
{
  cltab = list()
  for(I in 1:dim(mat)[1])
  {
    C = unlist(lapply(mat[I,], .splitValue, splitBy=splitBy))
    cltab[[I]] = as.vector(table(C), mode = "list")
  }

  rwtab = list()
  for(I in 1:dim(mat)[2])
  {
    R = unlist(lapply(mat[,I], .splitValue, splitBy=splitBy))
    rwtab[[I]] = as.vector(table(R), mode = "list")
  }

  creatDataFram <- function(inLst)
  {
  nColVal = unique(unlist(lapply(inLst, names)))
  rxt = data.frame(matrix(NA, nrow = length(inLst), ncol=length(nColVal)))
  colnames(rxt) = nColVal
  for(I in 1:length(inLst))
  {
    rx = sapply(nColVal, function(x){ w= inLst[[I]][[x]]
                                      if(is.null(w))
                                        {return(NA)} else{return(w)} })
    rxt[I,] = rx[nColVal]
  }
  return(rxt)
  }
  rwdf = creatDataFram(rwtab)
  cldf = creatDataFram(cltab)

  rwdf[is.na(rwdf)]=0
  cldf[is.na(cldf)]=0

  if(scaleRow==TRUE){ rwdf = t(apply(rwdf, 1, function(x)100*x/sum(x))) }
  if(scaleCol==TRUE){ cldf = t(apply(cldf, 1, function(x)100*x/sum(x))) }

  return(list(rowSt=cldf, colSt=rwdf))
}

creatSideBarPlot <- function(mat, colPalette, splitBy=";", scaleRow=TRUE, scaleCol=TRUE)
{
  rcDF = .calculatRowColStat(mat, splitBy, scaleRow=scaleRow, scaleCol=scaleCol)

  colorX = unlist(colPalette[colnames(rcDF$colSt)])
  colBar = anno_barplot(rcDF$colSt, which = "column", axis = TRUE, gp = gpar(fill = colorX))
  column_ha = HeatmapAnnotation(barplot = colBar)
  #foo1 = rcDF$colSt
  #column_ha = HeatmapAnnotation(foo1 = anno_barplot(foo1, axis = TRUE))


  colorX = unlist(colPalette[colnames(rcDF$rowSt)])
  rowbar = anno_barplot(rcDF$rowSt, which = "row", axis = TRUE, axis_side = "top", gp = gpar(fill = colorX))
  row_ha = rowAnnotation(row_anno_barplot=rowbar, width = unit(2, "cm"))
  return(list(colPlt= column_ha, rowPlt= row_ha))
}

.sortPlotMat <- function(mat, controlD =NA, control.col="green", drug.col="black")
{
  rwNM = rownames(mat); clNm = colnames(mat)
  rwNM = sort(rwNM); clNm = sort(clNm)

  ##---------for row ------------------------------------------
  controlD =c(controlD)
  if(length(controlD[!is.na(controlD)])>0)
  {
    nonCntr = rwNM[!(rwNM %in% controlD)]
    rwNMx = c(controlD, nonCntr)
    rwNameCol = c(rep(control.col, length(controlD)),
                  rep(drug.col, length(nonCntr) ))
  } else{
    rwNMx = rwNM
    rwNameCol = rep(drug.col, length(rwNM))
  }

  ##--------for column ------------------------------------------
  clNameCol = rep("black", length(clNm))
  rtx = list(mat= mat[rwNMx, clNm],
             row.name.col = rwNameCol,
             col.name.col = clNameCol)
  return(rtx)

}


##============================================================================
#' Plot mRECIST for models and drugs
#'
#' \code{plot.mRECIST} plots the mRECIST
#'
#' @param object The \code{Xeva} dataset
#' @param model.id The \code{model.id}
#'
#' @return plot
#'
#' @examples
#' data(pdxe)
#' #groupBy = "biobase.id"
#' control.name = c("untreated")
#' df = getmRECIST(pdxe)
#' df = df[1:500,]
#' plot.mRECIST(df,groupBy = "biobase.id", control.name = "untreated")
#'
#' @export
plot.mRECIST <- function(df,groupBy = "biobase.id", control.name = "untreated")
{

  control.name = c(control.name)
  mat = .castDataFram(df, row.var="drug.join.name", col.var = groupBy, value="mRECIST")
  matRC = .sortPlotMat(mat, controlD=control.name, control.col="green", drug.col="black")
  mat = as.matrix(matRC$mat)

  colPalette = list("CR" = "#4daf4a", "PR" = "#377eb8", "SD"= "#e41a1c", "PD"= "#984ea3")
  nameSpc = unique(as.vector(as.matrix(mat)))
  backgroundCol = "gray"
  bgCol = rep(backgroundCol, length(nameSpc))
  splitBy=";"
  sortCellValue = TRUE #FALSE

  sidePlt = creatSideBarPlot(mat, colPalette, splitBy=";", scaleRow=FALSE, scaleCol=FALSE)
  pltX = ComplexHeatmap::Heatmap(mat, name = "Drug & Models", col=bgCol,
                 top_annotation = sidePlt$colPlt, top_annotation_height = unit(2, "cm"),
                 cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE,
                 show_row_names = TRUE, row_names_side = "left",
                 row_names_gp = gpar(col = matRC$row.name.col),
                 show_column_names = TRUE, column_names_side = "top",
                 rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                 show_heatmap_legend = FALSE,
                 cell_fun =function(j, i, x, y, width, height, fill)
                 {.custom_cell_fun(x, y, width, height, mat[i,j], colPalette, fill, splitBy, sortCellValue)
                 }) + sidePlt$rowPlt


  colVec = unlist(colPalette)[names(colPalette)]
  HLeg = legendGrob(names(colPalette), pch=22,
                    gp=gpar(col = colVec, fill = colVec))
  #pdf(file="DATA-raw/mRECIST_plot.pdf", width=12, height=9)
  draw(pltX, heatmap_legend_list = list(HLeg), padding = unit(c(10, 10, 10, 10), "mm"))
  #dev.off()
}


