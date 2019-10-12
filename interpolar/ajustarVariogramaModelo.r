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

ajustarVariogramaModelo <- function(variogramaEmpirico) {
  direcciones <- unique(variogramaEmpirico$dir.hor)
  
  if (length(direcciones) == 1) {
    # Busqueda del mejor variograma modelo
    # Estimaciones iniciales para el NL-OLS
    # estimamos el psill inicial como el valor en el medio entre la media y el máximo de la serie
    psillIni = (median(variogramaEmpirico$gamma) + max(variogramaEmpirico$gamma)) * 0.5
    # estimamos el rango inicial como la primer distancia a la cual se supera el psill
    rIni = variogramaEmpirico$dist[which(variogramaEmpirico$gamma > psillIni)[1]]
    # estimamos p inicial (para el modelo potencia) como 0.2
    pIni = 0.2
    
    # para el nugget hacemos la recta que pasa mejor (mínimos cuadrados) por los primeros 3 puntos, la estimación del
    # nugget será donde la recta corte el eje y (x=0), pero solo se usará el nugget si es mas del 50% del psill inicial
    nuggetIni <- as.numeric(predict(lm(gamma ~ dist, data=variogramaEmpirico[1:3]), newdata=data.frame(dist=c(0))))
    usarNugget <- nuggetIni > psillIni * 0.5
    
    if (usarNugget) {
      vmExp <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Exp', psill=psillIni, nugget=nuggetIni))
      vmSph <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Sph', psill=psillIni, nugget=nuggetIni))
      vmGau <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Gau', psill=psillIni, nugget=nuggetIni))
      vmWhi <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Mat', psill=psillIni, nugget=nuggetIni, kappa=1))
      vmMat <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Mat', psill=psillIni, nugget=nuggetIni, kappa=1.9))
      vmPow <- fit.variogram(variogramaEmpirico, vgm(range=pIni, 'Pow', psill=psillIni, nugget=nuggetIni))
      vmLin <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Lin', psill=psillIni, nugget=nuggetIni))
      # vmLog <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Log', psill=psillIni, nugget=nuggetIni))
    } else {
      vmExp <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Exp', psill=psillIni))
      vmSph <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Sph', psill=psillIni))
      vmGau <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Gau', psill=psillIni))
      vmWhi <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Mat', psill=psillIni, kappa=1))
      vmMat <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Mat', psill=psillIni, kappa=1.9))
      vmPow <- fit.variogram(variogramaEmpirico, vgm(range=pIni, 'Pow', psill=psillIni))
      vmLin <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Lin', psill=psillIni))
      # vmLog <- fit.variogram(variogramaEmpirico, vgm(range=rIni, 'Log', psill=psillIni))
    }
    
    # elección del mejor según mínimo error cuadrático medio
    variogramasModelo <- list(vmExp, vmSph, vmGau, vmWhi, vmMat, vmPow, vmLin)
    SErrs <- sapply(variogramasModelo, attr, which='SSErr')
    mejorVariogramaModelo <- variogramasModelo[[which.min(SErrs)]]
    
    # Grafico el variograma empírico y los distintos modelos ajustados
    # maxDist <- max(variogramaEmpirico$dist)
    # plot(variogramaEmpirico$dist, variogramaEmpirico$gamma, ylim=c(0, max(variogramaEmpirico$gamma)), xlab='Distancia', ylab='Semivarianza')
    # leyenda <- character(length(variogramasModelo))
    # for (i in 1:length(variogramasModelo)) {
    #  lines(variogramLine(variogramasModelo[[i]], maxDist), lty=i, col=i)
    #  leyenda[i] <- paste(variogramasModelo[[i]]$model, collapse=" ")
    #}
    #legend('topleft', col=1:length(variogramasModelo), lty=1:length(variogramasModelo), legend=leyenda)
    #rm(leyenda, vmExp, vmSph, vmGau, vmWhi, vmMat, vmPow, vmLin, variogramasModelo, SErri, SErrMejorVariogramaModelo)
    
    rm(variogramasModelo, SErrs, vmExp, vmSph, vmGau, vmWhi, vmMat, vmPow, vmLin, 
       psillIni, rIni, pIni, nuggetIni, usarNugget)
    
    return(mejorVariogramaModelo)
  } else {
    variogramasDirecciones <- array(data=NA, dim=length(direcciones))
    anisotropiaGeometrica <- F
    for (i in direcciones) {
      variogramaDireccionI <- variogramaEmpirico[variogramaEmpirico$dir.hor == i,]
      variogramasDirecciones[i] <- ajustarVariogramaModelo(variogramaDireccionI)
      if (i > 1) {
        anisotropiaGeometrica <- anisotropiaGeometrica &&
                                 variogramasDirecciones[i - 1]$model == variogramasDirecciones[i]$model &&
                                 variogramasDirecciones[i - 1]$psill * 0.9 <= variogramasDirecciones[i]$psill &&
                                 variogramasDirecciones[i]$psill <= variogramasDirecciones[i - 1]$psill * 1.1 &&
                                 variogramasDirecciones[i - 1]$range * 0.9 <= variogramasDirecciones[i]$range &&
                                 variogramasDirecciones[i]$range <= variogramasDirecciones[i - 1]$range * 1.1      
      }
    }
    
    if (anisotropiaGeometrica) {
      
    }
  }
}