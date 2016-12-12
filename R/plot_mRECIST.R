library(ComplexHeatmap)
if(1==2){
####==============================================

library(ComplexHeatmap)
getTestMat = function()
{
  mat = matrix(NA, nrow = 5, ncol = 4)
  mat[1,] = c("CR", "PR", "SD", "PD")
  mat[2,] = c("CR;PR",NA, "SD;PD", "PD")
  mat[3,] = c("CR;PR;CR","SD;PD;PD;SD", "PD;PD", "PR")
  mat[4,] = c("PR","CR", "SD", "PR")
  mat[5,] = c("CR","SD;SD", "PR;PD", "SD")
  mat
}



###-------------------------------------------------
.splitmR <- function(mRx, sort=TRUE)
{
  if(is.na(mRx)){return(NA)}
  mRy = strsplit(mRx, ";")[[1]]
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

.cell_funXeva <- function(x, y, w, h, value, colPalette, backgroundCol)
{
  factR = 0.99
  wr=w*0.5*factR;   hr=h*0.5*factR
  x0=x-wr; x1=x+wr; y0=y-hr; y1=y+hr
  x0=convertX(x0, "npc", valueOnly = TRUE)
  x1=convertX(x1, "npc", valueOnly = TRUE)
  y0=convertX(y0, "npc", valueOnly = TRUE)
  y1=convertX(y1, "npc", valueOnly = TRUE)
  vx = .splitmR(value, sort = FALSE)
  filCol = unlist(colPalette[vx])
  N=length(vx)
  cordXY = getCellBoxCordi(x0,x1, y0, y1, N)
  cordXY$x = unit(cordXY$x,"npc"); cordXY$y = unit(cordXY$y,"npc")
  grid.polygon(x = cordXY$x, y = cordXY$y,
               id = rep(1:(N+1), each = 4),
               gp = gpar(fill = c(NA, filCol), col = NA))
}

plotHeatmapX <- function()
  {
    mat = getTestMat()
    colPalette = list("CR" = "blue", "PR" = "green", "SD" = "yellow", "PD"="red")
    nameSpc = unique(as.vector(mat))
    backgroundCol = "gray"
    bgCol = rep(backgroundCol, length(nameSpc))
    #z=list() ; saveRDS(z, file="cellFunValue.Rda")
    Heatmap(mat, name = "Drug & Models", col=bgCol,
            cluster_rows = FALSE, cluster_columns = FALSE,
            #column_dend_height = unit(5, "cm"),
            show_row_dend = FALSE,
            show_row_names = TRUE, row_names_side = "left",
            show_column_names = TRUE, column_names_side = "top",
            rect_gp = gpar(col = "white", lty = 1, lwd = 1),
            show_heatmap_legend = FALSE,
            cell_fun =function(j, i, x, y, width, height, fill)
            {.cell_funXeva(x, y, width, height, mat[i,j], colPalette, fill)
              } )


}


Heatmap(mat, name = "go", rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.rect(x, y, w, h, gp = gpar(fill = "#dcb35c", col = NA))
          if(i == 1) {
            grid.segments(x, y-h*0.5, x, y)
          } else if(i == nrow(mat)) {
            grid.segments(x, y, x, y+h*0.5)
          } else { grid.segments(x, y-h*0.5, x, y+h*0.5) }
          if(j == 1) {
            grid.segments(x, y, x+w*0.5, y)
          } else if(j == ncol(mat)) {
            grid.segments(x-w*0.5, y, x, y)
          } else { grid.segments(x-w*0.5, y, x+w*0.5, y)}

          if(i %in% c(4, 10, 16) & j %in% c(4, 10, 16)) {
            grid.points(x, y, pch = 16, size = unit(2, "mm"))
          }

          r = min(unit.c(w, h))*0.45
          if(is.na(mat[i, j])) {
          } else if(mat[i, j] == "W") {
            grid.circle(x, y, r, gp = gpar(fill = "red", col = "white"))
          } else if(mat[i, j] == "B") {
            grid.circle(x, y, r, gp = gpar(fill = "green", col = "black"))
          }
        },
        col = c("B" = "black", "W" = "white"),
        show_row_names = FALSE, show_column_names = FALSE,
        column_title = "One famous GO game",
        heatmap_legend_param = list(title = "Player", at = c("B", "W"),
                                    labels = c("player1", "player2"), grid_border = "black")
)




}
