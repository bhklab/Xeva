
getTestMat = function()
{
  mat = matrix(NA, nrow = 5, ncol = 4)
  mat[1,] = c("CR", "PR", "SD", "PD")
  mat[2,] = c("CR;PR",NA, "SD;PD", "PD")
  mat[3,] = c("CR;PR;CR","SD;PD;PD;SD", "PD;PD", "PR")
  mat[4,] = c("PR","CR", "SD", "PR")
  mat[5,] = c("CR","SD;SD", "PR;PD", "SD")
  rownames(mat)=paste0("Drug", seq_len(dim(mat)[1]) )
  colnames(mat)=paste0("Sample", seq_len(dim(mat)[2]) )
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

.custom_cell_fun <- function(x, y, w, h, value, colPalette, backgroundCol,
                             splitBy=";", sort=TRUE)
{
  factR = 01.0
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
               id = rep(seq_len(N+1), each = 4),
               gp = gpar(fill = c(NA, filCol),
                         col = "#f0f0f0"))
}



.calculatRowColStat <- function(mat, splitBy=";", scaleRow=TRUE, scaleCol=TRUE)
{
  cltab = list()
  for(I in seq_len(dim(mat)[1]))
  {
    C = unlist(lapply(mat[I,], .splitValue, splitBy=splitBy))
    cltab[[I]] = as.vector(table(C), mode = "list")
  }

  rwtab = list()
  for(I in seq_len(dim(mat)[2]))
  {
    R = unlist(lapply(mat[,I], .splitValue, splitBy=splitBy))
    rwtab[[I]] = as.vector(table(R), mode = "list")
  }

  creatDataFram <- function(inLst)
  {
    nColVal = unique(unlist(lapply(inLst, names)))
    rxt = data.frame(matrix(NA, nrow = length(inLst), ncol=length(nColVal)))
    colnames(rxt) = nColVal
    for(I in seq_along(inLst))
    {
      rx = vapply(nColVal, function(x){w= inLst[[I]][[x]]
                                      if(is.null(w))
                                      {return(NA)} else{return(w)} },
                  FUN.VALUE = numeric(1))
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



.creatSideBarPlot <- function(mat, colPalette, splitBy=";", scaleRow=TRUE,
                              scaleCol=TRUE)
{
  rcDF = .calculatRowColStat(mat, splitBy, scaleRow=scaleRow, scaleCol=scaleCol)

  colorX = unlist(colPalette[colnames(rcDF$colSt)])

  colBar = ComplexHeatmap::anno_barplot(rcDF$colSt, which = "column",
                                        axis = TRUE, gp = grid::gpar(fill = colorX))

  column_ha = ComplexHeatmap::HeatmapAnnotation(barplot = colBar,
                                                height = grid::unit(2, "cm"),
                                show_annotation_name = FALSE)

  colorX = unlist(colPalette[colnames(rcDF$rowSt)])

  #rowbar = ComplexHeatmap::anno_barplot(rcDF$rowSt, which = "row", axis = TRUE,
  #                                      axis_side = "top",
  #                                      #axis_param = list(side = "top"),
  #                                      gp = gpar(fill = colorX))


  rowbar = tryCatch({
                     ComplexHeatmap::anno_barplot(rcDF$rowSt, which = "row", axis = TRUE,
                                        axis_side = "top", gp = grid::gpar(fill = colorX))
    }, error = function(e) {
    ComplexHeatmap::anno_barplot(rcDF$rowSt, which = "row", axis = TRUE,
                                 axis_param = list(side = "top"), gp = grid::gpar(fill = colorX))

      })


  row_ha = ComplexHeatmap::rowAnnotation(row_anno_barplot=rowbar, width = grid::unit(2, "cm"),
                         show_annotation_name = FALSE)
  return(list(colPlt= column_ha, rowPlt= row_ha))
}

.sortPlotMat <- function(mat, controlD, control.col, drug.col)
{
  ##-------first sort by number of NA --------------------------
  rowNa <- apply(mat, 1, function(x)sum(is.na(x)))
  colNa <- apply(mat, 2, function(x)sum(is.na(x)))
  mat <- mat[names(sort(rowNa)), names(sort(colNa))]
  ##------------------------------------------------------------
  rwNM <- rownames(mat); clNm <- colnames(mat)
  ##---------for row ------------------------------------------
  controlD =c(controlD)
  if(length(controlD[!is.na(controlD)])>0)
  {
    nonCntr <- rwNM[!(rwNM %in% controlD)]
    rwNMx <- c(controlD, nonCntr)
    rwNameCol <- c(rep(control.col, length(controlD)),
                  rep(drug.col, length(nonCntr) ))
  } else{
    rwNMx <- rwNM
    rwNameCol <- rep(drug.col, length(rwNM))
  }
  ##--------for column ------------------------------------------

  if(length(controlD[!is.na(controlD)])>0)
  {
    contMat <- mat[controlD[1], clNm]
    clNm <- names(sort(contMat, na.last = TRUE))
  }
  clNameCol <- rep("black", length(clNm))

  rtx = list(mat= mat[rwNMx, clNm],
             row.name.col = rwNameCol,
             col.name.col = clNameCol)
  return(rtx)

}

.mRcolPalette <- function(mr)
{
  cp <- list("CR" = "#0033CC", "CR-->PD" = "#3182bd", "CR-->-->PD" = "#bf8ef2",
             "PR" = "#1a9850", "PR-->PD" = "#91cf60", "PR-->-->PD" = "#bfb35a",
             "SD" = "#fed976", "SD-->PD" = "#ffeda0", "SD-->-->PD" = "#fed976",
             "PD"= "#e41a1c")
  colPal <- cp[mr]
  colPal <- colPal[names(cp)]
  colPal <- colPal[!is.na(names(colPal))]
  return(colPal)
}

##============================================================================
#' To plot mRECIST values
#'
#' \code{plotmRECIST} plots the mRECIST matrix obtained from \code{summarizeResponse}.
#'
#' @param mat The mRECIST matrix where rows are drugs and columns are patients.
#' @param control.name Name of the control.
#' @param control.col Color of the control.
#' @param drug.col Color of the drug names.
#' @param colPalette Color palette for mRECIST values.
#' @param name Title of the plot.
#' @param sort If matrix should be sorted before plotting.
#' @param row_fontsize Size of the row name font.
#' @param col_fontsize Size of the column name font.
#' @param legend_title Title for the legend.
#' @param draw_plot Default \code{TRUE} will plot the figure. If \code{FALSE}, return an object.
#' @return mRECIST plot.
#' @examples
#' data(brca)
#' brca.mr <- summarizeResponse(brca, response.measure = "mRECIST", group.by="patient.id")
#' plotmRECIST(brca.mr, control.name = "untreated")
#' @export
#' @import ComplexHeatmap
#' @import grid
#' @importFrom grid gpar
plotmRECIST <- function(mat, control.name = NA, control.col="#238b45", drug.col="black",
                        colPalette = NULL, name = "Drug & Models", sort=TRUE,
                        row_fontsize=12, col_fontsize=12, legend_title="Response",
                        draw_plot=TRUE)
{
  if(!is(mat, "matrix"))
  {mat <- as.matrix(mat)}

  unqMat <- as.character(unique(as.vector(mat)))
  unqMat <- unique(unlist(strsplit(unqMat, ";")))
  unqMat <- unqMat[!is.na(unqMat)]
  
  if(is.null(colPalette))
  {
    colPalette <- .mRcolPalette(unqMat)
  } else
  {
    colPre <- vapply(unqMat, function(x) is.null(colPalette[[x]]), FUN.VALUE = logical(1))
    if(any(colPre)==TRUE)
    {
      colAbName <- names(colPre[colPre==TRUE])
      colAb <- paste(colAbName,collapse = "\n")
      msg1 = sprintf("color for these values are not present in colPalette\n%s", colAb)
      stop(msg1)
    }
  }

  control.name <- c(control.name)
  if(sort==TRUE)
  {
    matRC <-.sortPlotMat(mat, controlD=control.name, control.col=control.col,
                        drug.col=drug.col)
    mat<-as.matrix(matRC$mat)
  } else
  { mat <- as.matrix(mat) }

  rowColors <- rep(drug.col, nrow(mat))
  names(rowColors) <- rownames(mat)
  if(!is.na(control.name)){ rowColors[control.name] <- control.col}

  nameSpc = unique(as.vector(as.matrix(mat)))
  backgroundCol = "gray"
  bgCol = rep(backgroundCol, length(nameSpc))
  splitBy <- ";"
  sortCellValue = TRUE
  sidePlt <- .creatSideBarPlot(mat, colPalette, splitBy=splitBy, scaleRow=FALSE, scaleCol=FALSE)

  maxRWN <- rownames(mat)[nchar(rownames(mat))==max(nchar(rownames(mat)))][1]
  rwSide <- grobWidth(textGrob(maxRWN)) + unit(0, "mm")
  pltX <- ComplexHeatmap::Heatmap(mat, name = name, col=bgCol,
                 top_annotation = sidePlt$colPlt,
                 cluster_rows = FALSE, cluster_columns = FALSE,
                 show_row_dend = FALSE, show_row_names = TRUE,
                 row_names_side = "left", row_names_max_width = rwSide,
                 row_names_gp = gpar(col = rowColors, fontsize = row_fontsize),
                 show_column_names = TRUE, column_names_side = "top",
                 column_names_gp = gpar(fontsize = col_fontsize),
                 column_title = name,
                 rect_gp = gpar(col = "white", lty = 1, lwd = 1),
                 show_heatmap_legend = FALSE,
                 cell_fun =function(j, i, x, y, width, height, fill)
                 {.custom_cell_fun(x, y, width, height, mat[i,j], colPalette, fill, splitBy,
                                   sortCellValue)
                 } #) + sidePlt$rowPlt
                 , right_annotation = sidePlt$rowPlt )

  colVec <- unlist(colPalette)[names(colPalette)]
  HLeg <- Legend(at = names(colPalette), title = legend_title,
                 legend_gp = gpar(col = colVec, fill = colVec))

  if(draw_plot==TRUE)
  {
    padding = unit(c(2,2,2,2), "mm")
    draw(pltX, heatmap_legend_list = list(HLeg), padding = padding)

  } else{ return(list(plot=pltX, legend = HLeg))}

}








