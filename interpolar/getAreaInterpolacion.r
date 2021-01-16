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

getAreaInterpolacion <- function(params, proyeccionEntrada, proyeccionSalida, shpPais) {
  if (is.na(params$xMin) || is.na(params$xMax) || is.na(params$yMin) || is.na(params$yMax)) {
    if (is.null(shpPais)) {
      xRange <- shpPais@bbox[1,]
      yRange <- shpPais@bbox[2,]
    } else {
      xRange <- range(sp::coordinates(observaciones)[,1])
      yRange <- range(sp::coordinates(observaciones)[,2])
    }
    
    xMin <- xRange[1]
    xMax <- xRange[2]
    yMin <- yRange[1]
    yMax <- yRange[2]

    transformarRangos <- F
    aplicarHolgura <- T  
  } else {
    xMin <- params$xMin
    xMax <- params$xMax
    yMin <- params$yMin
    yMax <- params$yMax
    
    aplicarHolgura <- F
    transformarRangos <- T
  }

  if (transformarRangos) {
    esquinas <- c(xMin, yMax,  #topLeft
                  xMax, yMax, #topRight
                  xMin, yMin, #bottomLeft
                  xMax, yMin #bottomRight
                  )
    esquinas <- matrix(data=c(esquinas), ncol=2, byrow=T)
    esquinas <- data.frame(esquinas=esquinas)
    
    sp::coordinates(esquinas) <- c('X1', 'X2')
    proj4string(esquinas) <- proyeccionEntrada
    
    esquinas <- spTransform(esquinas, CRSobj=proyeccionSalida)
  
    xMin <- min(esquinas$X1)
    xMax <- max(esquinas$X1)
    yMin <- min(esquinas$X2)
    yMax <- max(esquinas$X2)    
  }
  if (aplicarHolgura) {  
    # dejo un borde de 1% para cada lado
    deltaX <- (xMax - xMin) * 0.01
    xMin <- trunc(xMin - deltaX)
    xMax <- ceiling(xMax + deltaX)
    deltaY <- (yMax - yMin) * 0.01
    yMin <- trunc(yMin - deltaY)
    yMax <- ceiling(yMax + deltaY)
  }
  
  return(list(xlim=c(xMin, xMax), ylim=c(yMin, yMax)))
}