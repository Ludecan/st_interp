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

# Busca en el stack el path relativo a getwd() del script para poder hacer los source correctamente
# Intenta con el frame actual - 3 primero que es donde estaba siempre cuando se hizo el programa
iFrame <- sys.nframe()
if (iFrame >= 3) { script.dir.parsearParamsMosaicAndClip <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.parsearParamsMosaicAndClip <- NULL }
while ((is.null(script.dir.parsearParamsMosaicAndClip) || is.na(regexpr('parsearParamsMosaicAndClip.r', script.dir.parsearParamsMosaicAndClip, fixed=T)[1])) && iFrame >= 0) {
  script.dir.parsearParamsMosaicAndClip <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.parsearParamsMosaicAndClip)) { script.dir.parsearParamsMosaicAndClip <- ''
} else { script.dir.parsearParamsMosaicAndClip <- paste(dirname(script.dir.parsearParamsMosaicAndClip), '/', sep='') }

source(paste(script.dir.parsearParamsMosaicAndClip, '../parsearParams/parsearParamsUtils.r', sep=''))

createParamsMosaic <- function(pathEjecucion='./',
                               pathProceso='',
                               escala=1,
                               offset=0) {
  return(list(pathEjecucion=pathEjecucion,
              pathProceso=pathProceso,
              escala=escala, offset=offset))
}

parsearParamsMosaic <- function(params) {
  return(getParamValuesFromConstructorParams(params, funcCrearParams=createParamsMosaic))
}


createParamsMosaicAndClip <- function(pathEjecucion='./',
                                      pathProceso='',
                                      proj4stringResultados,
                                      xMin, xMax, yMin, yMax,
                                      escala=1,
                                      offset=0) {
  return(list(pathEjecucion=pathEjecucion,
              pathProceso=pathProceso,
              proj4stringResultados=proj4stringResultados,
              xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax,
              escala=escala, offset=offset))
}

parsearParamsMosaicAndClip <- function(params) {
  return(getParamValuesFromConstructorParams(params, funcCrearParams=createParamsMosaicAndClip))
}

armarMosaicoI <- function(i, pathsTiles, escala=1, offset=0) {
  nombresArchivos <- character(ncol(pathsTiles) - 1)
  for (iArchivo in 1:length(nombresArchivos)) {
    nombresArchivos[iArchivo] <- nombreArchSinPathNiExtension(pathsTiles[i, iArchivo + 1])
  }
  
  extension <- getFileExt(pathsTiles[i, 2])
  
  n <- min(nchar(nombresArchivos[1]), nchar(nombresArchivos[2], keepNA = T), na.rm=T)
  iPrefijo <- 1
  prefijosHastaIIguales <- TRUE
  while (iPrefijo <= n && prefijosHastaIIguales) {
    j <- 2
    while (j <= length(nombresArchivos) & prefijosHastaIIguales) {
      prefijosHastaIIguales <- substr(nombresArchivos[1], 1, iPrefijo) == substr(nombresArchivos[j], 1, iPrefijo)
      j <- j + 1
    }
    
    if (prefijosHastaIIguales) iPrefijo <- iPrefijo + 1
  } 
  if (iPrefijo > n)
    iPrefijo <- n
  prefijo <- substr(nombresArchivos[1], 1, iPrefijo - 1)
  nombreArchivoSalida <- paste(dirname(pathsTiles[i, 2]), '/', prefijo, 
                               paste(gsub(pattern = prefijo, replacement = '', x = nombresArchivos, fixed = T), collapse='_'), 
                               '.', extension, sep='')
  
  if (!file.exists(nombreArchivoSalida)) {
    mosaico <- raster(pathsTiles[i, 2])
    for (j in 3:ncol(pathsTiles)) {
      rasterJ <- raster(pathsTiles[i, j])
      mosaico <- mosaic(mosaico, rasterJ, fun=mean)
    }
    
    if (escala != 1 | offset != 0) {
      valores <- getValues(mosaico)
      valores <- valores * escala + offset
      mosaico <- setValues(mosaico, valores)
    }
    writeRaster(mosaico, filename = nombreArchivoSalida, overwrite=TRUE, options = c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9'))
  }
  return (nombreArchivoSalida)
}
