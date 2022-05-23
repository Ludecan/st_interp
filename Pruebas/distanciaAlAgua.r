############################################################################
# Spatial Interpolation Lib by Pablo Alfaro (pabloalfaropineyro@gmail.com) #
############################################################################
# A complete library for performing spatial interpolation on space time    #
# series data                                                              #
# Allows incorporation of gridded auxiliary fields to assist on            #
# interpolation modelling via techniques from the Regression Kriging family#
# as well as provides functionalities for plotting and storing the results.#
# This library was developed as the graduation project for the author's    #
# Masters Degree on Mathematical Engineering                               #
############################################################################
# Copyright (C) 2012 Pablo Alfaro                                          #
# This program is free software: you can redistribute it and/or modify it  #
# under the terms of the GNU General Public License as published by the    #
# Free Software Foundation, either version 3 of the License, or (at your   #
# option) any later version.                                               #
# This program is distributed in the hope that it will be useful, but      #
# WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General #
# Public License for more details.                                         #
# You should have received a copy of the GNU General Public License along  #
# with this program. If not, see http://www.gnu.org/licenses/.             #
############################################################################

source("D:/Workspace/MCH2/MCH/MCHLib/R/Scripts/instalarPaquetes/instant_pkgs.r")
source("D:/Workspace/MCH2/MCH/MCHLib/R/Scripts/interpolar/interpolarEx.r")
source("D:/Workspace/MCH2/MCH/MCHLib/R/Scripts/interpolar/mapearEx.r")
instant_pkgs(c('Rcpp', 'raster'))


pathRasterLC <- 'C:/mch/ArchivosProcesosLocales/MascaraDeTierra/LCType_PuntaDelEsteCorregido.tif'
lc <- raster(pathRasterLC)
#shpBase <- cargarSHP(pathSHP = 'C:/mch/ArchivosProcesosLocales/CartografiaBase/uruguay_departamentos.shp', proj4strSHP = '+proj=longlat +datum=WGS84 +no_defs')

bbAreaAInterpolar <- bbox(coordsAInterpolar)
bbObservaciones <- bbox(coordsObservaciones)

# Obtengo los 4 vértices del área que engloba las coordenadas a interpolar y las observaciones
bbGeneral <- bbAreaAInterpolar
bbGeneral[1, 1] <- min(bbGeneral[1, 1], bbObservaciones[1, 1])
bbGeneral[2, 1] <- min(bbGeneral[2, 1], bbObservaciones[2, 1])
bbGeneral[1, 2] <- max(bbGeneral[1, 2], bbObservaciones[1, 2])
bbGeneral[2, 2] <- max(bbGeneral[2, 2], bbObservaciones[2, 2])

# Extiendo esa área 20% de su largo y ancho hacia cada borde
bbGeneralExt <- bbGeneral
pctExtender <- 0.2
extX <- diff(bbGeneralExt[1,]) * pctExtender
extY <- diff(bbGeneralExt[2,]) * pctExtender
bbGeneralExt[1,] <- bbGeneralExt[1,] + c(-extX, extX)
bbGeneralExt[2,] <- bbGeneralExt[2,] + c(-extY, extY)

if (!file.exists('C:/mch/ArchivosProcesosLocales/MascaraDeTierra/cuerposDeAguaUy_25MKm2.shp')) {
  bbToSpatialPolygon <- function(bb, crsOBJ) {
    coordsBB <- matrix(data=c(bb[1, 1], bb[2, 1], # abajo, izquierda
                              bb[1, 1], bb[2, 2], # arriba, izquierda
                              bb[1, 2], bb[2, 2], # arriba, derecha
                              bb[1, 2], bb[2, 1], # abajo, derecha
                              bb[1, 1], bb[2, 1]), # abajo, izquierda
                       ncol=2, byrow=T)
    p <- Polygon(coordsBB)
    ps <- Polygons(list(p), 1)
    sps <- SpatialPolygons(list(ps), proj4string=crsOBJ)
    return(sps)
  }
  
  bbPoligono <- bbToSpatialPolygon(bb = bbGeneralExt, crsOBJ = coordsAInterpolar@proj4string)
  bbPoligono <- spTransform(bbPoligono, crs(lc))
  
  lcRecortado <- crop(lc, bbPoligono)
  cuerposDeAgua <- sp::disaggregate(rasterToPolygons(lcRecortado, fun =  function(x) { x == 0 }, dissolve = T))
  cuerposDeAgua$area_MKm2 <- area(cuerposDeAgua) / 1E6
  cuerposDeAgua2 <- cuerposDeAgua[cuerposDeAgua$area_MKm2 >= 25,]
  
  spplot(cuerposDeAgua, zcol = 2)
  spplot(cuerposDeAgua2, zcol = 2)
  
  colnames(cuerposDeAgua@data)[1] <- 'LCType'
  colnames(cuerposDeAgua2@data)[1] <- 'LCType'
  
  guardarSHP(cuerposDeAgua, 'C:/mch/ArchivosProcesosLocales/MascaraDeTierra/cuerposDeAguaUy.shp')  
  guardarSHP(cuerposDeAgua2, 'C:/mch/ArchivosProcesosLocales/MascaraDeTierra/cuerposDeAguaUy_25MKm2.shp')
} else {
  cuerposDeAgua2 <- cargarSHP('C:/mch/ArchivosProcesosLocales/MascaraDeTierra/cuerposDeAguaUy_25MKm2.shp')
}


# Creo una grilla para el área extendida que contenga a coordsAInterpolar, con su misma resolución y que
# los bordes de coordsAInterpolar caigan exactamente en un borde da la grilla extendida
# Estos dos valores deben ser enteros
#(bbAreaAInterpolar[1, 1] - xmin) / cellSize[1]
#(bbAreaAInterpolar[2, 1] - ymin) / cellSize[2]
gr <- coordsAInterpolar@grid
cellSize <- gr@cellsize

xmin <- bbAreaAInterpolar[1, 1] - floor((bbAreaAInterpolar[1, 1] - bbGeneralExt[1, 1]) / cellSize[1]) * cellSize[1]
ymin <- bbAreaAInterpolar[2, 1] - floor((bbAreaAInterpolar[2, 1] - bbGeneralExt[2, 1]) / cellSize[2]) * cellSize[2]
cellsDim <- c(ceiling((bbGeneralExt[1, 2] - xmin) / cellSize[1]), ceiling((bbGeneralExt[2, 2] - ymin) / cellSize[2]))
grExtendida <- SpatialGridDataFrame(
  grid = GridTopology(cellcentre.offset = c(xmin, ymin), cellsize = cellSize, cells.dim = cellsDim),
  data = data.frame(distanciaAlAgua=rep(NA, cellsDim[1] * cellsDim[2])), 
  proj4string = coordsAInterpolar@proj4string)

cuerposDeAgua2 <- spTransform(cuerposDeAgua2, coordsAInterpolar@proj4string)
grExtendida$distanciaAlAgua <- t(gDistance(spgeom1 = as(grExtendida, 'SpatialPoints'), spgeom2 = gUnaryUnion(cuerposDeAgua2), byid = T))[,1]

escala <- crearEscalaEquiespaciada(grExtendida$distanciaAlAgua, nIntervalos = 8, brewerPal = 'PuBuGn', continuo = T)
mapearGrillaGGPlot(grilla=grExtendida, shpBase = spTransform(cuerposDeAgua, grExtendida@proj4string), escala = escala)

writeGDAL(dataset = grExtendida, fname = 'C:/mch/ArchivosProcesosLocales/MascaraDeTierra/distanciaAlAguaUy_25MKm2.tif', options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))

combinarRasters <- function(pathsRasters, archSalida='test.tif') {
  require(raster)
  pathsRasters <- dir('C:/mea/', pattern = '*.tif$', full.names = T)
  archSalida='C:/mea/GMTED2010_Mean_Uy.tif'
  rasters <- vector(mode = "list", length = length(pathsRasters))
  for (i in 1:length(pathsRasters)) {
    rasters[[i]] <- raster(pathsRasters[i])
  }
  lapply(rasters, range)
  
  rasters$filename <- archSalida
  rasters$overwrite <- TRUE
  m <- do.call(merge, rasters)
  plot(m)
  require(rgdal)
  n <- 10
}
