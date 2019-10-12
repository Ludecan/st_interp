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

function(params, proyeccionEntrada, proyeccionSalida, shpPais) {
  if (is.na(params$xMin) || is.na(params$xMax) || is.na(params$yMin) || is.na(params$yMax)) {
    if (is.null(shpPais)) {
      xRange <- shpPais@bbox[1,]
      yRange <- shpPais@bbox[2,]
    } else {
      xRange <- range(coordinates(observaciones)[,1])
      yRange <- range(coordinates(observaciones)[,2])
    }
    
    params$xMin <- xRange[1]
    params$xMax <- xRange[2]
    params$yMin <- yRange[1]
    params$yMax <- yRange[2]

    transformarRangos <- F
    aplicarHolgura <- T  
  } else {
    aplicarHolgura <- F
    transformarRangos <- T
  }

  if (transformarRangos) {
    params$xMin <- -59
    params$xMax <- -53
    params$yMin <- -35
    params$yMax <- -30
    
    esquinas <- c(params$xMin, params$yMax,  #topLeft
                  params$xMax, params$yMax, #topRight
                  params$xMin, params$yMin, #bottomLeft
                  params$xMax, params$yMin #bottomRight
                  )
    esquinas <- matrix(data=c(esquinas), ncol=2, byrow=T)
    esquinas <- data.frame(esquinas)
    
    coordinates(esquinas) <- c('X1', 'X2')
    proj4string(esquinas) <- proyeccionEntrada
    
    esquinas <- spTransform(esquinas, CRSobj=proyeccionSalida)
  
    params$xMin <- min(esquinas$X1)
    params$xMax <- max(esquinas$X1)
    params$yMin <- min(esquinas$X2)
    params$yMax <- max(esquinas$X2)    
  }
  if (aplicarHolgura) {  
    # dejo un borde de 1% para cada lado
    deltaX <- (params$xMax - params$xMin) * 0.01
    params$xMin <- trunc(params$xMin - deltaX)
    params$xMax <- ceiling(params$xMax + deltaX)
    deltaY <- (params$yMax - params$yMin) * 0.01
    params$yMin <- trunc(params$yMin - deltaY)
    params$yMax <- ceiling(params$yMax + deltaY)
  }
}