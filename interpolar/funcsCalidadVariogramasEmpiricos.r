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

# A maximizar
gradoDeMonotoneidadCreciente <- function(variograma) {
  grado <- 0
  for (i in 1:(nrow(variograma)-1)) {
    if (variograma$gamma[i] <= variograma$gamma[i+1]) {
      grado <- grado * 2 + 1
    } else {
      grado <- grado * 2
    }
  }
  return(grado / (2^(nrow(variograma)-1) - 1))
}

rangoEstimado <- function(variograma) {
  psillEstimado = (median(variograma$gamma) + max(variograma$gamma)) * 0.5
  i <- max(which(variograma$gamma>=psillEstimado)[1]-1,1)
  
  if (i < nrow(variograma)-1) {
    return (variograma$dist[i] + (psillEstimado - variograma$gamma[i]) / (variograma$gamma[i+1] - variograma$gamma[i]) * (variograma$dist[i+1] - variograma$dist[i]))
  } else {
    return (variograma$dist[i])
  }
}

# A minimizar
pendienteAlInicio <- function(variograma, maxIntervalo=3) {
  # fuerzo la recta que pasa por el 0 y lo mejor posible por los 3 primeros puntos del variograma
  modelo <- lm(gamma~dist-1,data=variograma[1:maxIntervalo])
  return(modelo$coefficients[1])
}

pendienteGeneral <- function(variograma) {
  # fuerzo la recta que pasa por el 0 y lo mejor posible por todos los puntos del variograma
  modelo <- lm(gamma~dist-1,data=variograma)
  return(modelo$coefficients[1])
}

varianzaEstimadaAcumulada <- function(variograma) {
  return(sum(variograma$np/variograma$dist^2))
}

varianzaAlrededorDePendienteGeneral <- function(variograma) {
  # fuerza a pasar por el 0 y toma la recta que mejor ajuste a los puntos
  # la varianza es el error cuadrático medio del ajuste, minimizarla implica tomar el variograma empírico 
  # que mejor se ajuste a la función usada como formula en lm. 
  modelo <- lm(gamma~dist-1,data=variograma,weights=variograma$np/variograma$dist^2)
  tendencia <- predict.lm(modelo,newdata=data.frame(dist=variograma$dist))
  return(var(variograma$gamma-tendencia))
}

varianzaAlrededorDeModeloExponencial <- function(variograma) {
  # Por mas que fuerze a un modelo exponencial, tiene más la forma que queremos, es más parecido que una recta
  # tiene problemas pq no siempre converge. Hay que revisar eso
  psillEstimado = (median(variograma$gamma) + max(variograma$gamma)) * 0.5
  i <- max(which(variograma$gamma>=psillEstimado)[1]-1,1)
  if (i < nrow(variograma)-1) {
    rangoEstimado <- (variograma$dist[i] + (psillEstimado - variograma$gamma[i]) / (variograma$gamma[i+1] - variograma$gamma[i]) * (variograma$dist[i+1] - variograma$dist[i]))
  } else {
    rangoEstimado <- (variograma$dist[i])
  }
  modelo <- nls(gamma~C*(1-I(exp(-3*dist/a))), data=variograma, start=list(a=rangoEstimado, C=psillEstimado), 
                control=list(maxiter = 100, tol = 1e-04, minFactor = 1/1024, printEval = FALSE, warnOnly = T))
  tendencia <- predict(modelo,newdata=data.frame(dist=variograma$dist))
  if (modelo$convInfo$isConv) {
    return(var(variograma$gamma-tendencia))
  } else {
    return(.Machine$double.xmax)
  }
}

ratioNuggetPsill <- function(variograma) {
  # estimate psill and nugget
  nuggetEst <- as.numeric(predict(lm(gamma~dist,data=variograma[1:3]),newdata=data.frame(dist=c(0))))
  psillEst <- (median(variograma$gamma) + max(variograma$gamma)) * 0.5  
  return (nuggetEst/psillEst) 
}
