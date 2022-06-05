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

seleccionarMetodoInterpolacion <- function(valoresObservados) {
  # adaptado de la descripción de interpolate de la documentación de intamap
  dataObs <- valoresObservados
  minObs <- min(dataObs)
  if (minObs <= 0)
    dataObs <- dataObs + minObs + sd(dataObs)
  
  g = boxcox(dataObs ~ 1,lambda=seq(-2.5,2.5,len=101),plotit=FALSE)$y
  tkfn <- fivenum(dataObs)
  rangoIntercuartil <- IQR(dataObs)
  test <- length(boxplot.stats(dataObs)$out)/length(dataObs) > 0.1 ||
    tkfn[3] - tkfn[2] < rangoIntercuartil/3 ||
    tkfn[4] - tkfn[3] < rangoIntercuartil/3 ||
    g[71] < sort(g)[91]
  
  if (test)
    return("copula")
  else
    return("automap")
}