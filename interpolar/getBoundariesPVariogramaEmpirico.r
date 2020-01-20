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
if (iFrame >= 3) { script.dir.getBoundariesPVariogramaEmpirico <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.getBoundariesPVariogramaEmpirico <- NULL }
while ((is.null(script.dir.getBoundariesPVariogramaEmpirico) || is.na(regexpr('getBoundariesPVariogramaEmpirico.r', script.dir.getBoundariesPVariogramaEmpirico, fixed=T)[1])) && iFrame >= 0) {
  script.dir.getBoundariesPVariogramaEmpirico <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.getBoundariesPVariogramaEmpirico)) { script.dir.getBoundariesPVariogramaEmpirico <- ''
} else { script.dir.getBoundariesPVariogramaEmpirico <- paste0(dirname(script.dir.getBoundariesPVariogramaEmpirico), '/') }

getProporcionDeLaDiagonalValida <- function(nObservaciones) {
  if (nObservaciones <= 8) { return(1)
  } else if (nObservaciones <= 15) { return(0.9)
  } else if (nObservaciones <= 50) { return(0.8)
  } else if (nObservaciones <= 80) { return(0.7)
  } else if (nObservaciones <= 150) { return(0.5)
  } else { return(0.35) }
}

getBoundariesPVariogramaEmpirico <- function(fml, observaciones, minDivDist=5, 
                                             maxDivDist=14, proporcionDeLaDiagonalValida=0) {
  require(gstat)
  source(paste0(script.dir.getBoundariesPVariogramaEmpirico, 'funcsCalidadVariogramasEmpiricos.r'))
  
  if (proporcionDeLaDiagonalValida <= 0) { proporcionDeLaDiagonalValida <- getProporcionDeLaDiagonalValida(nrow(observaciones)) }
  
  longlat <- !is.projected(observaciones)
  if(is.na(longlat)) longlat <- FALSE
  diagonal <- spDists(t(bbox(observaciones)), longlat=longlat)[1,2] * proporcionDeLaDiagonalValida # times the length of the central axis through the area
  
  # para funcs de maximización
  #valorOptimo <- -.Machine$double.xmax
  #for (i in maxDivDist:minDivDist) {
  #  ancho <- diagonal / i
  #  boundaries <- seq(from=ancho, to=diagonal, by=ancho)
  #  v <- variogram(fml, observaciones, boundaries=boundaries)
  #  valor <- gradoDeMonotoneidadCreciente(v)
  #  #valor <- rangoEstimado(v)
  #  if (valor > valorOptimo) {
  #    res <- boundaries
  #    valorOptimo <- valor
  #  }
  #}

  # try different amounts of intervals
  # para funcs de minimización
  valorOptimo <- .Machine$double.xmax
  i<-maxDivDist
  i<- i-1
  for (i in maxDivDist:minDivDist) {
    ancho <- diagonal / i
    boundaries <- seq(from=ancho, to=diagonal, by=ancho)
    v <- variogram(fml, observaciones, boundaries=boundaries)
    #plot(v)
    # valor <- pendienteGeneral(v) 
    # valor <- varianzaEstimadaAcumulada(v)
    
    # the best empirical variogram is the one with the least squared error around the best fit line that goes through 0,0
    # this way we get a variogram which honors the actual observation as much as it can, is monotonically growing and is given by the data
    # ask Tufa about this approach
    valor <- gradoDeMonotoneidadCreciente(variograma=v)
    # varianzaAlrededorDePendienteGeneral(variograma=v)
    # valor <- varianzaAlrededorDeModeloExponencial(v)
    # valor <- ratioNuggetPsill(variograma=v)
    if (valor < valorOptimo) {
      res <- boundaries
      valorOptimo <- valor
    }
  }
  return(res)
}

getBoundariesPVariogramaEmpiricoV2 <- function(fml,observaciones,minDivDist=5, 
                                               maxDivDist=14, proporcionDeLaDiagonalValida=0) {
  # the best empirical variogram is the one with the least squared error around the trend going through 0,0
  # this way the "ideal" variogram honors the actual observation as much as it can, is monotonically growing and is given by the data
  require(gstat)
  if (proporcionDeLaDiagonalValida <= 0) { proporcionDeLaDiagonalValida <- getProporcionDeLaDiagonalValida(nrow(observaciones)) }
  
  longlat <- !is.projected(observaciones)
  if(is.na(longlat)) longlat <- FALSE
  diagonal <- spDists(t(bbox(observaciones)), longlat=longlat)[1,2] * proporcionDeLaDiagonalValida # times the length of the central axis through the area

  vCloud <- variogram(fml, observaciones, cloud=T, cutoff=diagonal)
  # logarithmic trendline through 0 at 0.001
  # f <- function(x,a, b) {a * (log(x, base=b) - a * log(0.01, base=b)}
  # f <- function(x, a) {a * log(x/0.001)}
  # modeloTendencia <- nls(gamma ~ f(dist, a), start=c(a=1), data=vCloud)
  
  # linear trendline through 0
  # the best empirical variogram is the one with the least squared error around the best fit line that goes through 0,0
  # this way we get a variogram which honors the actual observation as much as it can, is monotonically growing and is given by the data
  modeloTendencia <- lm(gamma~dist-1, data=vCloud)
  # plot(vCloud$dist, vCloud$gamma)
  # lines(sort(vCloud$dist), predict(modeloTendencia, newdata=data.frame(x=sort(sort(vCloud$dist)))), col="blue")

  # try different amounts of intervals
  valorOptimo <- .Machine$double.xmax
  #i<-maxDivDist
  #i <- maxDivDist
  #i <- 5
  for (i in maxDivDist:minDivDist) {
    ancho <- diagonal / i
    boundaries <- seq(from=ancho, to=diagonal, by=ancho)
    v <- variogram(fml, observaciones, boundaries=boundaries)
    tendencia <- predict(modeloTendencia, newdata=data.frame(dist=v$dist))
    
    # print(plot(v))
    # plot(x=v$dist, y=tendencia, ylim=c(0,max(tendencia)))
    # lines(x=v$dist, y=tendencia)
    
    valor <- sum((v$gamma-tendencia)^2)
    if (valor < valorOptimo) {
      res <- boundaries
      valorOptimo <- valor
    }
  }
  
  return(res)
}

getBoundariesPVariogramaEmpiricoV3 <- function(fml, observaciones, proporcionDeLaDiagonalValida=0) {
  require(gstat)
  if (proporcionDeLaDiagonalValida <= 0) { proporcionDeLaDiagonalValida <- getProporcionDeLaDiagonalValida(nrow(observaciones)) }
  
  longlat <- !is.projected(observaciones)
  if(is.na(longlat)) longlat <- FALSE
  diagonal <- spDists(t(bbox(observaciones)), longlat=longlat)[1,2] * proporcionDeLaDiagonalValida # times the length of the central axis through the area
  
  # return all distances in the variogram cloud, equivalent to fitting directly to the cloud
  vCloud <- variogram(fml, observaciones, cloud=T, cutoff=diagonal)
  return(sort(unique(vCloud$dist)))
}

fitLinearSplitVariogram <- function(x, y) {
  # la funcion a ajustar es la siguiente
  # f(x) <- { ax  si x < Cx
  #         { aCx   si x >= Cx
  # pasa por el origen con pendiente a hasta Cx y luego de Cx sigue una constante horizontal igual a aCx
  
  f <- function (Cx) {
    a <- coefficients(lm(y ~ ifelse(x < Cx, x, Cx) - 1))[1]
    return(ifelse(x < Cx, a * x, a * Cx))
  }
  r2 <- function(Cx) { sum((y - f(Cx))^2) }
  
  rango <- range(x)
  # Como mínimo Cx vale la distancia a la que hay al menos 30 elementos o la mínima distancia mas 25% del recorrido
  minimo <- max(rango[1] + (rango[2] - rango[1]) * 0.25, sort(x)[30], na.rm=T)
  # Como máximo puede valer el máximo de x menos 10% del recorrido
  maximo <- rango[2] - (rango[2] - rango[1]) * 0.1

  # sapply(seq(minimo, maximo, length.out = 100), r2)
  
  res <- optimize(r2, interval=c(minimo, maximo))
  best_Cx <- res$minimum

  a <- coefficients(lm(y ~ ifelse(x < best_Cx, x, best_Cx) - 1))[1]
  res <- function(x) {
    return(ifelse(x < best_Cx, a * x, a * best_Cx))
  }
  # plot(x, y)
  # plot(x, res(x))
  return(res)
}

getBoundariesPVariogramaEmpiricoV4 <- function(fml, observaciones, minDivDist=5, maxDivDist=14, 
                                               proporcionDeLaDiagonalValida=0, cutoff=NA) {
  # the best empirical variogram is the one with the least squared error around a linear variogram fit to the variogram cloud and with no nugget
  # this way the "ideal" variogram honors the actual observation as much as it can, is monotonically growing and is given by the data
  require(gstat)
  if (proporcionDeLaDiagonalValida <= 0) { proporcionDeLaDiagonalValida <- getProporcionDeLaDiagonalValida(nrow(unique(coordinates(observaciones)))) }
  
  if (is.na(cutoff)) {
    longlat <- !is.projected(observaciones)
    if(is.na(longlat)) longlat <- FALSE
    cutoff <- max(spDists(observaciones, longlat=longlat)) * proporcionDeLaDiagonalValida
  } 
  
  vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  modeloTendencia <- fitLinearSplitVariogram(x = vCloud$dist, y = vCloud$gamma)
  
  
  
  # colores <- brewer.pal(4, 'Set1')  
  # plot(vCloud$dist, vCloud$gamma, col=colores[1], pch='*', ylim=c(0, max(predict(modeloTendencia, newdata=data.frame(x=vCloud$dist)), vCloud$gamma)))
  # lines(sort(vCloud$dist), predict(modeloTendencia, newdata=data.frame(x=sort(sort(vCloud$dist)))), col='black')

  # i <- maxDivDist
  # i <- round((minDivDist + maxDivDist)/2)
  # i <- minDivDist
  # i <- 5
  
  # try different amounts of intervals
  valorOptimo <- .Machine$double.xmax
  for (i in maxDivDist:minDivDist) {
    ancho <- cutoff / i
    boundaries <- seq(from=ancho, to=cutoff, by=ancho)
    
    v <- variogram(fml, observaciones, boundaries=boundaries)
    tendencia <- modeloTendencia(v$dist)
    
    # if (i == maxDivDist) { print(points(x = v$dist, y=v$gamma, col=colores[2], pch=15))
    # } else if (i == round((minDivDist + maxDivDist)/2)) { print(points(x = v$dist, y=v$gamma, col=colores[3], pch=16))
    # } else if (i == minDivDist) { print(points(x = v$dist, y=v$gamma, col=colores[4], pch=17)) }
    
    valor <- sum((v$gamma-tendencia)^2)
    if (valor < valorOptimo) {
      res <- boundaries
      valorOptimo <- valor
    }
  }
  
  return(res)
}

getBoundariesPVariogramaEmpiricoV5 <- function(fml, observaciones, minNIntervalos=5, maxNIntervalos=14, 
                                               proporcionDeLaDiagonalValida=0, cutoff=NA, 
                                               minProporcionPuntosPorBin=1/maxNIntervalos) {
  # the best empirical variogram is the one with the least squared error around a linear variogram fit to the variogram cloud and with no nugget
  # but also keep only the divisions with at least minProporcionPuntosPorBin in each bin
  # this way the "ideal" variogram honors the actual observation as much as it can, is monotonically growing and is given by the data
  require(gstat)
  if (is.na(cutoff)) {
    if (proporcionDeLaDiagonalValida <= 0) { proporcionDeLaDiagonalValida <- getProporcionDeLaDiagonalValida(nrow(unique(coordinates(observaciones)))) }
    longlat <- !is.projected(observaciones)
    if(is.na(longlat)) longlat <- FALSE
    cutoff <- max(spDists(observaciones, longlat=longlat)) * proporcionDeLaDiagonalValida # times the length of the central axis through the area
  } 
  
  vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  # Si el cutoff indicado no deja al menos 30 pares de puntos en la nube, agrando el cutoff
  if (autoAjustarCutoff && (nrow(vCloud) < 30)) {
    longlat <- !is.projected(observaciones)
    if (is.na(longlat)) longlat <- FALSE
    cutoff <- max(spDists(observaciones, longlat=longlat))
    vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  }  
  modeloTendencia <- fitLinearSplitVariogram(x = vCloud$dist, y = vCloud$gamma)
  minNParesPuntos <- minProporcionPuntosPorBin * length(vCloud$dist)
  
  # colores <- brewer.pal(4, 'Set1')  
  # plot(vCloud$dist, vCloud$gamma, col=colores[1], pch='*', ylim=c(0, max(predict(modeloTendencia, newdata=data.frame(x=vCloud$dist)), vCloud$gamma)))
  # lines(sort(vCloud$dist), predict(modeloTendencia, newdata=data.frame(x=sort(sort(vCloud$dist)))), col='black')
  
  # try different amounts of intervals
  valorOptimo <- .Machine$double.xmax
  sortedDist <- sort(vCloud$dist)
  
  i <- maxNIntervalos
  for (i in maxNIntervalos:minNIntervalos) {
    nPuntosPorIntervalo <- length(sortedDist) / i
    
    if (i == minNIntervalos || nPuntosPorIntervalo >= minNParesPuntos) {
      iDists <- round(seq.int(from=nPuntosPorIntervalo, to=length(sortedDist), length.out = i))
      boundaries <- sortedDist[iDists]
      
      v <- variogram(fml, observaciones, boundaries=boundaries)
      tendencia <- modeloTendencia(v$dist)
      
      # if (i == minNIntervalos + 2) { print(points(x = v$dist, y=v$gamma, col=colores[2], pch=15))
      # } else if (i == minNIntervalos + 1) { print(points(x = v$dist, y=v$gamma, col=colores[3], pch=16))
      # } else if (i == minNIntervalos) { print(points(x = v$dist, y=v$gamma, col=colores[4], pch=17)) }
      
      valor <- sum((v$gamma-tendencia)^2)
      if (valor < valorOptimo) {
        res <- boundaries
        valorOptimo <- valor
      }
    }
  }

  return(res)
}

getBoundariesPVariogramaEmpiricoV6 <- function(fml, observaciones, minNIntervalos=5, maxNIntervalos=14, 
                                               proporcionDeLaDiagonalValida=0, cutoff=NA, autoAjustarCutoff=T,
                                               minProporcionPuntos1erBin=1/maxNIntervalos) {
  # the best empirical variogram is the one with the least squared error around a linear variogram fit to the variogram cloud and with no nugget
  # but also keep at least minProporcionPuntos1erBin % in the 1st bin
  # this way the "ideal" variogram honors the actual observation as much as it can, is monotonically growing and is given by the data  
  require(gstat)
  if (is.na(cutoff)) {
    if (proporcionDeLaDiagonalValida <= 0) { proporcionDeLaDiagonalValida <- getProporcionDeLaDiagonalValida(nrow(unique(coordinates(observaciones)))) }
    longlat <- !is.projected(observaciones)
    if (is.na(longlat)) longlat <- FALSE
    cutoff <- max(spDists(observaciones, longlat=longlat)) * proporcionDeLaDiagonalValida
  }
  # proporcionDeLaDiagonalValida <- 1
  vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  # plot(vCloud)
  # Si el cutoff indicado no deja al menos 30 pares de puntos en la nube, agrando el cutoff
  if (autoAjustarCutoff && (nrow(vCloud) < 30)) {
    longlat <- !is.projected(observaciones)
    if (is.na(longlat)) longlat <- FALSE
    cutoff <- max(spDists(observaciones, longlat=longlat))
    vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  }
  modeloTendencia <- fitLinearSplitVariogram(x = vCloud$dist, y = vCloud$gamma)
  
  sDist <- sort(vCloud$dist)
  minDist <- sDist[round(minProporcionPuntos1erBin * length(vCloud$dist))]
  
  # try different amounts of intervals
  valorOptimo <- .Machine$double.xmax
  # i <- maxNIntervalos
  for (i in maxNIntervalos:minNIntervalos) {
    boundaries <- seq(from=minDist, to=cutoff, length.out = i)
    
    v <- variogram(fml, observaciones, boundaries=boundaries)
    tendencia <- modeloTendencia(v$dist)
    
    pesos <- v$np / v$dist^2 / sum(v$np / v$dist^2)
    #pesos <- v$np / v$dist / sum(v$np / v$dist)
    valor <- sum(pesos * (v$gamma-tendencia)^2)
    
    # print(c(i, valor))
    #print(ggplot(data = data.frame(dist = c(0, v$dist), semivariance = c(0, v$gamma), tendencia = c(0, tendencia)), aes(dist)) + geom_point(aes(y=semivariance), colour="black") + geom_line(aes(y=tendencia), colour="red") + scale_x_continuous(breaks=seq(from = 0, to = cutoff, length.out = 11)) + scale_y_continuous(breaks=seq(from = 0, to = max(v$gamma), length.out = 11)) + labs(title = paste('i = ', i, ', valor = ', valor)))
    
    if (valor < valorOptimo) {
      res <- boundaries
      valorOptimo <- valor
    }
  }
  
  return(res)
}

getBoundariesPVariogramaEmpiricoV7 <- function(fml, observaciones, minNIntervalos=5, maxNIntervalos=14, 
                                               minProporcionPuntos1erBin=1/maxNIntervalos, autoAjustarCutoff=T) {
  # the best empirical variogram is the one with the least squared error around a linear variogram fit to the variogram cloud and with no nugget
  # but also keep at least minProporcionPuntos1erBin % in the 1st bin and weight the rmse to give priority to bins with more observations and closer distances
  # this way the "ideal" variogram honors the actual observation as much as it can, is monotonically growing and is given by the data  
  require(gstat)
  longlat <- !is.projected(observaciones)
  if (is.na(longlat)) longlat <- FALSE
  dists <- spDists(observaciones, longlat = longlat)
  dists <- dists[upper.tri(dists)]
  x <- sort(dists)
  y <- sapply(x, FUN = function(x) { return(sum(dists < x)) })
  # plot(x, y / x)
  cutoff <- x[which.max(y / x)]
  vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  # plot(vCloud)
  # Si el cutoff indicado no deja al menos 30 pares de puntos en la nube, agrando el cutoff
  if (autoAjustarCutoff && (nrow(vCloud) < 30)) {
    cutoff <- max(dists)
    vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  }
  modeloTendencia <- fitLinearSplitVariogram(x = vCloud$dist, y = vCloud$gamma)
  
  sDist <- sort(vCloud$dist)
  minDist <- sDist[round(minProporcionPuntos1erBin * length(vCloud$dist))]
  
  # try different amounts of intervals
  valorOptimo <- .Machine$double.xmax
  # i <- 7
  for (i in maxNIntervalos:minNIntervalos) {
    boundaries <- seq(from=minDist, to=cutoff, length.out = i)
    
    v <- variogram(fml, observaciones, boundaries=boundaries)
    tendencia <- modeloTendencia(x=v$dist)
    # print(ggplot(data = data.frame(dist = v$dist, semivariance = v$gamma, tendencia = tendencia), aes(dist)) + geom_point(aes(y=semivariance), colour="black") + geom_line(aes(y=tendencia), colour="red"))
    
    pesos <- v$np / v$dist^2 / sum(v$np / v$dist^2)
    valor <- sum(pesos * (v$gamma-tendencia)^2)
    if (valor < valorOptimo) {
      res <- boundaries
      valorOptimo <- valor
    }
  }
  
  return(res)
}

getBoundariesPVariogramaEmpiricoV8 <- function(fml, observaciones, minNIntervalos=5, maxNIntervalos=50, 
                                               proporcionDeLaDiagonalValida=0, cutoff=NA, autoAjustarCutoff=T,
                                               minNPPorBin=ceil(length(observaciones)^(2/3))) {
  # the best empirical variogram is the one with the least squared error around a linear variogram fit to the variogram cloud and with no nugget
  # but also keep at least minProporcionPuntos1erBin % in the 1st bin
  # this way the "ideal" variogram honors the actual observation as much as it can, is monotonically growing and is given by the data  
  require(gstat)
  require(zoo)
  if (is.na(cutoff) || is.infinite(cutoff)) {
    if (proporcionDeLaDiagonalValida <= 0) { proporcionDeLaDiagonalValida <- getProporcionDeLaDiagonalValida(nrow(unique(coordinates(observaciones)))) }
    longlat <- !is.projected(observaciones)
    if (is.na(longlat)) longlat <- FALSE
    cutoff <- max(spDists(observaciones, longlat=longlat)) * proporcionDeLaDiagonalValida
  }
  
  vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  # plot(vCloud)
  # Si el cutoff indicado no deja al menos 30 pares de puntos en la nube, agrando el cutoff
  if (autoAjustarCutoff && (nrow(vCloud) < 30)) {
    longlat <- !is.projected(observaciones)
    if (is.na(longlat)) longlat <- FALSE
    cutoff <- max(spDists(observaciones, longlat=longlat))
    vCloud <- variogram(fml, observaciones, cloud=T, cutoff=cutoff)
  }
  
  sDist <- sort(vCloud$dist)
  logSDist <- log(sDist)
  sqrtSDist <- sqrt(sDist)
  
  invProporcionObservaciones <- 2
  semik = round(length(vCloud$gamma) / invProporcionObservaciones) %/% 2
  modeloTendencia <- rollmean(vCloud$gamma, k = 2 * semik + 1, fill = 'extend')
  modeloTendencia <- modeloTendencia * logSDist / logSDist[length(logSDist) - 1]
  # modeloTendencia <- modeloTendencia * sqrtSDist / sqrtSDist[length(sqrtSDist) - 1]
  
  # modeloTendencia <- rollmean(c(rep(0, semik), vCloud$gamma), k = 2 * semik + 1)
  # modeloTendencia <- c(modeloTendencia, rep(modeloTendencia[length(modeloTendencia)], semik))
  # modeloTendencia <- modeloTendencia * log(sDist) / max(log(sDist))
  # modeloTendencia <- modeloTendencia * mean(vCloud$gamma) / mean(modeloTendencia)
  # plot(vCloud$gamma)
  # lines(modeloTendencia, col='red')
  # lines( modeloTendencia * logSDist / logSDist[length(logSDist) - 1], col='blue')
  # lines( modeloTendencia * sqrtSDist / sqrtSDist[length(sqrtSDist) - 1], col='green')
  
  gammas <- numeric(length(modeloTendencia))
  
  # try different amounts of intervals
  valorOptimo <- .Machine$double.xmax
  # i <- maxNIntervalos
  # maxNIntervalos <- 100
  i <- 24
  for (i in maxNIntervalos:minNIntervalos) {
    boundaries <- seq(from=cutoff/i, to=cutoff, length.out = i)
    #boundaries <- sDist[round(seq(from=1/i, to=1, length.out = i) * length(vCloud$dist))]
    
    v <- variogram(fml, observaciones, boundaries=boundaries)
    
    if (min(v$np) >= minNPPorBin | i == minNIntervalos) {
      j <- 1
      k <- 1
      while (j <= length(sDist)) {
        while (j <= length(sDist) & (sDist[j] <= v$dist[k] | k == length(v$dist))) {
          gammas[j] = v$gamma[k]
          j <- j + 1
        }
        k <- k + 1
      }
      
      valor <- sum((gammas - modeloTendencia)^2 / sqrtSDist)
      
      # print(ggplot(data = data.frame(dist = c(0, sDist), semivariance = c(0, gammas), tendencia = c(0, modeloTendencia)), aes(dist)) + geom_point(aes(y=semivariance), colour="black") + geom_line(aes(y=tendencia), colour="red") + scale_x_continuous(breaks=seq(from = 0, to = cutoff, length.out = 11)) + scale_y_continuous(breaks=seq(from = 0, to = max(v$gamma), length.out = 11)) + labs(title = paste('i = ', i, ', valor = ', valor)))
      # print(c(i, valor, ifelse(valor < valorOptimo, 'New best', '')))
      
      if (valor < valorOptimo) {
        res <- boundaries
        valorOptimo <- valor
      }     
    }
  }
  
  return(res)
}
