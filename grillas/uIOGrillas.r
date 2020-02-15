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

iFrame <- sys.nframe()
if (iFrame >= 3) { script.dir.uIOGrillas <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.uIOGrillas <- NULL }
while ((is.null(script.dir.uIOGrillas) || is.na(regexpr('uIOGrillas.r', script.dir.uIOGrillas, fixed=T)[1])) && iFrame >= 0) {
  script.dir.uIOGrillas <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.uIOGrillas)) { script.dir.uIOGrillas <- ''
} else { script.dir.uIOGrillas <- paste(dirname(script.dir.uIOGrillas), '/', sep='') }

source(paste0(script.dir.uIOGrillas, '../instalarPaquetes/instant_pkgs.r'), encoding = 'WINDOWS-1252')
instant_pkgs(c('sp', 'rgdal'))

# para todas las funciones grilla tiene que ser un objeto spatialPixelsDataFrame
guardarDefinicionGrilla <- function(archivoDefinicion, grilla) {
  if (!is.na(proj4string(grilla))) {
    p4str <- proj4string(grilla)
  } else {
    p4str <- ''
  }
  xs <- unique(grilla@coords[,1])
  ys <- unique(grilla@coords[,2])
  
  lineas <- array(data=c(
    paste('proj4string=', p4str, ';', sep=''),
    paste('ejeX=', paste(as.character(xs),collapse=','), ';', sep=''),
    paste('ejeY=', paste(as.character(xs),collapse=','), ';', sep='')))
  writeLines(lineas, con=archivoDefinicion)
}

guardarDatosGrilla <- function(archivoDatos, grilla) {
  writeBin(con=archivoDatos, object=grilla$value)
}

guardarGrillaArchDefYArchDatoBin <- function(archivoDefinicion, grilla) {
  guardarDefinicionGrilla(archivoDefinicion, grilla)
  # cambio la extensión por .bin
  archivoDatos <- paste(c(substr(archivoDefinicion, 1, tail(unlist(gregexpr("\\.", archivoDefinicion)), 1) - 1), '.bin'), collapse='')
  guardarDatosGrilla(archivoDatos, grilla)
}

guardarGrillaGDAL <- function(nombreArchivo, grilla) {
  # el tipo de archivo queda determinado por nombreArchivo, p.ej: si es .tif guarda un geotiff
  writeGDAL(grilla, nombreArchivo, options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
}

leerGrillaGDAL <- function(nombreArchivo) {
  # el tipo de archivo queda determinado por nombreArchivo, p.ej: si es .tif lee un geotiff
  return (readGDAL(nombreArchivo, silent = T))
}

leerDefinicionGrilla <- function(archivoDefinicion) {
  con <- file(archivoDefinicion, "r", blocking = FALSE)
  lineas <- readLines(con, n=3)
  close(con)
  
  tokens <- array(data=NA, dim=c(length(lineas), 2))
  for (i in 1:length(lineas)) {
    posIgual <- regexpr('=', lineas[i], fixed=T)[1]
    tokens[i, 1] <- substr(x=lineas[i], start=1, stop=posIgual - 1)
    tokens[i, 2] <- substr(x=lineas[i], start=posIgual+1, stop=nchar(lineas[i])-1)
  }
  
  ejeX <- as.double(unlist(strsplit(tokens[2,2], ",")))
  ejeY <- as.double(unlist(strsplit(tokens[3,2], ",")))
  
  grilla <- expand.grid(x=ejeX, y=ejeY)
  coordinates(grilla) <- c('x', 'y')
  if (tokens[1,2] != "")
    proj4string(grilla) <- tokens[1,2]
  gridded(grilla) <- T
  return (grilla)
}

leerDatosGrilla <- function(archivoDatos, grilla) {
  con <- file(archivoDefinicion, "rb")
  grilla$value <- readBin(con, double(), n=nrow(grilla))
  close(con)
  return (grilla)
}

leerGrilla <- function(archivoDefinicion) {
  grilla <- leerDefinicionGrilla(archivoDefinicion)
  # cambio la extensión por .bin
  archivoDatos <- paste(c(substr(archivoDefinicion, 1, tail(unlist(gregexpr("\\.", archivoDefinicion)), 1) - 1), '.bin'), collapse='')
  leerDatosGrilla(archivoDatos, grilla)
  return (grilla)
}

leerGrillaBin <- function(archiBin, naValue=NULL) {
  binFile <- file(archiBin, open="rb")
  nRows <- readBin(con=binFile, what='integer', n=1)
  if (nRows > 0) {
    nCols <- readBin(con=binFile, what='integer', n=1)
    values <- matrix(readBin(con=binFile, what='numeric', n=nRows*nCols, size=8), nrow=nRows, ncol=nCols, byrow=T)
    if (!is.null(naValue)) { values[values==naValue] <- NA }
  } else {
    nCols <- 0
    values <- NULL
  }
  xMin <- readBin(con=binFile, what='numeric', n=1)
  xMax <- readBin(con=binFile, what='numeric', n=1)
  yMin <- readBin(con=binFile, what='numeric', n=1)
  yMax <- readBin(con=binFile, what='numeric', n=1)
  
  #nChars <- readBin(con=binFile, what='integer', n=1)
  proj4stringGrilla <- readBin(con=binFile, what='character')
  close(binFile)
  x <- seq(from=xMin, to=xMax, length.out=nCols)
  y <- seq(from=yMin, to=yMax, length.out=nRows)  
  
  gr <- expand.grid(x=x, y=y)
  result <- data.frame(x=gr$x, y=gr$y, values=values)
  coordinates(result) <- c('x', 'y')
  proj4string(result) <- proj4stringGrilla
  gridded(result) <- T
  
  return (result)
}

leerDatosGrillaBin <- function(archiBin, naValue=NULL) {
  binFile <- file(archiBin, open="rb")
  nRows <- readBin(con=binFile, what='integer', n=1)
  if (nRows > 0) {
    nCols <- readBin(con=binFile, what='integer', n=1)
    values <- matrix(readBin(con=binFile, what='numeric', n=nRows*nCols, size=8), nrow=nRows, ncol=nCols, byrow=T)
    if (!is.null(naValue)) { values[values==naValue] <- NA }
  } else {
    values <- NULL
  }
  close(binFile)
  return (values)
}

guardarGrillaBin <- function(archiBin, grilla, naValue=NULL) {
  binFile <- file(archiBin, open="wb")
  
  coords <- coordinates(grilla)
  x <- unique(coords[,1])
  y <- unique(coords[,2])

  values <- grilla@data[,1]
  if (!is.null(naValue)) { values[is.na(values)] <- naValue }
  
  writeBin(length(y), con=binFile)
  writeBin(length(x), con=binFile)
  writeBin(values, con=binFile)
  writeBin(min(x), con=binFile)
  writeBin(max(x), con=binFile)
  writeBin(min(y), con=binFile)
  writeBin(max(y), con=binFile)
  
  proj4stringGrilla <- proj4string(grilla)
  writeBin(nchar(proj4stringGrilla), con=binFile)
  writeBin(proj4stringGrilla, con=binFile)
  close(binFile)
}

guardarRasterBin <- function(archiBin, grillaRaster, naValue=NULL) {
  binFile <- file(archiBin, open="wb")
  ext <- extent(grillaRaster)
  
  # R stores matrices in column-major order, so we need to transpose them to get row-major order
  # values <- c(t(matrix(data=grillaRaster@data@values, nrow=grillaRaster@nrows, ncol=grillaRaster@ncols)))
  values <- grillaRaster@data@values
  if (!is.null(naValue)) { values[is.na(values)] <- naValue }
  
  writeBin(grillaRaster@nrows, con=binFile)
  writeBin(grillaRaster@ncols, con=binFile)
  writeBin(values, con=binFile)
  writeBin(ext@xmin, con=binFile)
  writeBin(ext@xmax, con=binFile)
  writeBin(ext@ymin, con=binFile)
  writeBin(ext@ymax, con=binFile)
  proj4stringGrilla <- proj4string(grillaRaster)
  #writeBin(nchar(proj4stringGrilla), con=binFile)
  writeBin(proj4stringGrilla, con=binFile)
  close(binFile)
}

guardarSPobj <- function(archivoSalida, objSp, naValue=NULL, formatoSalida=c('kml', 'netCDF', 'GeoTiff')) {
  formatoSalida <- formatoSalida[1]
  if (formatoSalida == 'kml') {
    stop('uIOGrillas.guardarSPobj: formatoSalida ".kml" no implementado')
    #writeGDAL(objSp, changeFileExt(archivoI, '.kml'))
  } else  if (formatoSalida == 'netCDF') {
    stop('uIOGrillas.guardarSPobj: formatoSalida "netCDF" no implementado')
  } else  if (formatoSalida == 'GeoTiff') {
    writeGDAL(objSp, changeFileExt(archivoI, '.tif'), options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
  } else
    stop(paste('uIOGrillas.guardarSPobj: formatoSalida desconocido ', formatoSalida, sep=''))
}