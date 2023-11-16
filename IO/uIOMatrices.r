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

getDefaultSizeOf <- function(what='double') {
  if (what == 'numeric' | what == 'double' | what == 'complex') {
    return(8L)
  } else if (what == 'integer' | what == 'int') {
    return(4L)
  } else if (what == 'character') {
    return(2L)
  } else if (what == 'logical' | what == 'raw') {
    return(1L)
  } else 
    stop(paste('uIOMatrices.getDefaultSizeOf: parámetro what desconocido ', what, sep = ''))
}

salvarVectorABinarioSinDimensiones <- function(pathArchivo, vector, NAValue=-.Machine$double.xmax, recordSize=getDefaultSizeOf(typeof(vector))) {
  vector[is.na(vector)] <- NAValue
  writeBin(vector, pathArchivo, size=recordSize)
}

salvarMatrizABinarioSinDimensiones <- function(pathArchivo, matriz, NAValue=-.Machine$double.xmax, recordSize=getDefaultSizeOf(typeof(matriz)), byRow=T) {
  matriz[is.na(matriz)] <- NAValue
  if (byRow) { writeBin(as.vector(t(matriz)), pathArchivo, size=recordSize)
  } else { writeBin(as.vector(matriz), pathArchivo, size=recordSize) }
}

leerMatrizDeBinario <- function(pathArchivo, NAValue=-.Machine$double.xmax, byRow=T, tipoCon='file', what='double', size=getDefaultSizeOf(what), row_names=NULL, col_names=NULL) {
  if (tipoCon == 'file') { arch <- file(pathArchivo, "rb")
  } else if (tipoCon == 'gzfile') { arch <- gzfile(pathArchivo, "rb")
  } else { stop(paste('Tipo de conexion desconocida', tipoCon)) }
  
  tryCatch(expr = {
    nFilas <- readBin(arch, what='integer', n = 1L)
    nColumnas <- readBin(arch, what='integer', n = 1L)
    datos <- readBin(arch, what = what, n = nFilas * nColumnas, size = size)
  }, finally = close(arch))
  
  datos <- matrix(data = datos, nrow = nFilas, ncol = nColumnas, byrow = byRow)  
  datos[datos == NAValue] <- NA
  
  if (!is.null(row_names)) { rownames(datos) <- row_names }
  if (!is.null(col_names)) { colnames(datos) <- col_names }
  return(datos)
}

salvarMatrizABinario <- function(pathArchivo, matriz, NAValue=-.Machine$double.xmax, recordSize=getDefaultSizeOf(typeof(matriz)), byRow=T) {
  iNA <- is.na(matriz)
  if (any(iNA)) matriz[is.na(matriz)] <- NAValue
  f <- file(description=pathArchivo, open="wb", blocking=T)
  m <- nrow(matriz)
  writeBin(m, f)
  if (m > 0) {
    n <- ncol(matriz)
    writeBin(n, f)
    if (n > 0) {
      if (byRow) { writeBin(as.vector(t(matriz)), f, size=recordSize)
      } else { writeBin(as.vector(matriz), f, size=recordSize) }
    }
  }
  close(f)
}