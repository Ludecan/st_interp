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
if (iFrame >= 3) { script.dir.qcTests <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.qcTests <- NULL }
while ((is.null(script.dir.qcTests) || is.na(regexpr('qcTests.r', script.dir.qcTests, fixed=T)[1])) && iFrame >= 0) {
  script.dir.qcTests <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.qcTests)) { script.dir.qcTests <- ''
} else { script.dir.qcTests <- paste(dirname(script.dir.qcTests), '/', sep='') }

source(paste(script.dir.qcTests, '../instalarPaquetes/instant_pkgs.r', sep=''))
source(paste(script.dir.qcTests, '../interpolar/interpolarEx.r', sep=''))
instant_pkgs(c('sp', 'robustbase'))

# Códigos para los distintos tipos de outliers detectables por los métodos
TTO_SinProblemasDetectados = 0L
TTO_OutlierPorLoBajo = 1L
TTO_OutlierPorLoAlto = 2L
TTO_PrecipitacionAislada=3L
TTO_SequedadAislada=4L
TTO_ValorNulo = -1L 
TTO_SinDatosSuficientesEnLaEstacion = -2L
TTO_SinDatosSuficientesEnVecinos = -3L
TTO_ValorConfirmado = -4L

stringsTTO <- c('Valor Confirmado', 'Sin Datos Suficientes En Vecinos', 'Sin Datos Suficientes En La Estación', 
                'Valor Nulo', 'Sin Problemas Detectados', 'Outlier Por Lo Bajo', 'Outlier Por Lo Alto', 
                'Precipitación Aislada', 'Sequedad Aislada')

TTipoOutlierToString <- function(tipoOutlier) {
  return(stringsTTO[tipoOutlier - TTO_ValorConfirmado + 1])
}

tiposOutliersValoresSospechosos <- c(TTO_OutlierPorLoBajo, TTO_OutlierPorLoAlto, TTO_PrecipitacionAislada, TTO_SequedadAislada)

between <- function(x, minInclusivo, maxExclusivo) { return(minInclusivo <= x & x < maxExclusivo) }

simpleInvDistanceWeightingI <- function(dataI, distsI, distPower=2, RnRMaskThreshold=0.3) {
  pesos <- 1 / (distsI ^ distPower)
  pesos <- pesos / sum(pesos)
  if (RnRMaskThreshold > 0) { rnr <- as.integer(as.integer(dataI > 0) %*% pesos > RnRMaskThreshold)
  } else { rnr <- 1L }
    
  return(as.numeric(as.numeric(dataI) %*% pesos * rnr))
}

simpleInvDistanceWeighting <- function(data, dists, idpRange=seq(1, 15, by = 0.5), RnRMaskThreshold=0.3) {
  #data=valoresVecinos
  #dists=distsVecinos
  # iFecha <- 24
  # iFecha <- nrow(data)-1
  # iFecha <- which(row.names(data) == '2017-10-01 00:00')

  res <- numeric(nrow(data))
  rmses <- numeric(length = length(idpRange))
  iNoNA <- !is.na(data)
  iFecha <- 29
  for (iFecha in 1:nrow(data)) {
    #print(iFecha)
    dataI <- as.numeric(data[iFecha, iNoNA[iFecha,]])
    distsI <- dists[iNoNA[iFecha,]]
    estimsCV <- numeric(length(dataI))
    
    if (length(dataI) > 1 & var(dataI) > 1E-3) {
      for (idp in seq_along(idpRange)) {
        #print(paste(iFecha, idp))
        for (i in seq_along(dataI)) {
          #print(paste(iFecha, idp, i))
          estimsCV[i] <- simpleInvDistanceWeightingI(dataI = dataI[-i], distsI = distsI[-i], 
                                                     distPower = idpRange[idp], RnRMaskThreshold = RnRMaskThreshold)
        }
        rmses[idp] <- mean((estimsCV - dataI)^2)
      }
    } else { rmses[] <- 0 }
    
    res[iFecha] <- simpleInvDistanceWeightingI(dataI = dataI, distsI = distsI, distPower = idpRange[which.min(rmses)], RnRMaskThreshold = RnRMaskThreshold)
  }
  
  return(res)
}

deteccionGradienteEnPuntos <- function(coordsObservaciones, iPuntoATestear, maxDist, 
                                       iesPuntosVecinos=which(gWithinDistance(coordsObservaciones[iPuntoATestear,], coordsObservaciones, dist = maxDist, byid = T)),
                                       umbralGradiente=5, mapear=FALSE, shpBase=NULL) {
  # iPuntoATestear <- i
  # coordsObservaciones$value <- - coordsObservaciones$value
  # row.names(coordsObservaciones)[iPuntoATestear]
  iesPuntosVecinos <- iesPuntosVecinos[!is.na(coordsObservaciones$value[iesPuntosVecinos])]
  coordsObservacionesVecinos <- coordsObservaciones[iesPuntosVecinos, ]
  # mapearPuntosGGPlot(coordsObservacionesVecinos, continuo = T, dibujar = F, zcol='value')
  # coordsObservaciones$Nombre <- row.names(coordsObservaciones)
  # mapearPuntosConEtiquetasGGPlot(puntos = coordsObservacionesVecinos, shpBase = shpBase, dibujar = F, zcol='Nombre')
  
  if (length(coordsObservacionesVecinos) >= 4 && var(coordsObservacionesVecinos$value) > 1E-6) {
    coords <- coordinates(coordsObservacionesVecinos)
    
    if (F) {
      rX <- range(coords[,1])
      dX <- rX[2] - rX[1]
      rY <- range(coords[,2])
      dY <- rY[2] - rY[1]
      maxD <- max(dX, dY)
      
      coords[,1] <- (coords[,1] - rX[1]) / maxD
      coords[,2] <- (coords[,2] - rY[1]) / maxD
    } else { 
      rX <- c(0, 0)
      rY <- c(0, 0)
      maxD <- 1 
    }
    
    dfAux <- data.frame(x = coords[,1], y = coords[,2], value = coordsObservacionesVecinos$value)
    if (!tryCatch(expr = { modelo <- rlm(formula = as.formula('value~x+y+1'), data=dfAux, maxit=50) 
    TRUE },
    warning = function(x) { FALSE }, error = function(x) { FALSE })) {
      modelo <- lm(formula = as.formula('value~x+y+1'), data=dfAux)
    }
    
    coefX <- coefficients(modelo)['x']
    coefY <- coefficients(modelo)['y']
    m <- as.numeric(coefY / coefX)
    creceConX <- (sign(coefX) > 0 & sign(coefY) > 0) | (sign(coefX) > 0 & sign(coefY) < 0)
    
    mPerpendicular <- -1/m
    invM <- -mPerpendicular
    
    # grilla <- grillaSobreBoundingBox(shpBase)
    # grilla <- SpatialGridDataFrame(grilla, data=data.frame(value=predict(modelo, newdata=data.frame(x=(coordinates(grilla)[,1] - rX[1]) / maxD, y=(coordinates(grilla)[,2] - rY[1]) / maxD))))
    # mapearGrillaGGPlot(grilla, shpBase, continuo = T, dibujar = F)
    
    # Sea r la recta definida por y = m * x + n y P el punto (x1, y1)
    # La perpendicular a r, r2, por P cumple que
    # y = -1/m * x + n2
    # y1 = -1/m * x1 + n2
    # Entonces
    # n2 = y1 + 1/m * x1
    # El punto donde intersectan r y r2 cumple
    # m * x + n = -1/m * x + n2
    # Entonces, el x de la proyección ortogonal de P sobre r es:
    # x = (n2 - n) / (m + 1/m)
    # n podría ser cualquiera ya que solo nos importa la dirección pero para que todo quede dentro del área del gráfico
    # definimos n para que r pase por el punto medio del área del gráfico Q => 
    # n = yq + -m * xq
    Q <- (rowMeans(bbox(coordsObservaciones)) - c(rX[1], rY[1])) / maxD
    n <- Q[2] - m * Q[1]
    
    n2s <- coords[,2] + invM * coords[,1]
    xs <- as.numeric(sort((n2s - n) / (m + invM)))
    
    iAux <- 1:(length(xs)-1)
    # Tomo los medios de los intervalos entre los xs originales proyectados sobre la dirección perpendicular al gradiente
    # y calculo sus ys
    nuevosXs <- (xs[iAux] + xs[iAux+1]) * 0.5
    nuevosYs <- nuevosXs * m + n
    
    # cbind(nuevosXs, nuevosYs)
    # cbind(row.names(coordsObservacionesVecinos), coordsObservacionesVecinos$value)
    
    #diferenciaEntreMediasDeCadaLado <- numeric(length(nuevosXs))
    medioValores <- mean(range(coordsObservacionesVecinos$value))
    nDatosConSignoCorrectoACadaLado <- matrix(NA_integer_, nrow=length(nuevosXs), ncol=2)
    # mediasACadaLado <- matrix(NA_real_, nrow=length(nuevosXs), ncol=2)
    # sdsACadaLado <- matrix(NA_real_, nrow=length(nuevosXs), ncol=2)
    # j <- 7
    for (j in seq_along(nuevosXs)) {
      #xs[j] * maxD + rX[1]
      x0y0 <- c(nuevosXs[j], nuevosYs[j])
      
      if ((!creceConX & mPerpendicular >= 0) | (creceConX & mPerpendicular < 0)) {
        iesLadoCalido <- coords[, 2] >= mPerpendicular * (coords[, 1] - x0y0[1]) + x0y0[2]
        iesLadoFrio <- !iesLadoCalido
      } else {
        iesLadoCalido <- coords[, 2] < mPerpendicular * (coords[, 1] - x0y0[1]) + x0y0[2]
        iesLadoFrio <- !iesLadoCalido
      }
      
      #print(mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpBase, zcol = 'value', continuo = T, dibujarTexto = T,
      #                         dibujar = F, titulo = j) + geom_abline(slope = mPerpendicular, intercept = x0y0[1]/m + x0y0[2]))
      
      nDatosConSignoCorrectoACadaLado[j, 1] <- sum(coordsObservacionesVecinos$value[iesLadoCalido] >= medioValores)
      nDatosConSignoCorrectoACadaLado[j, 2] <- sum(coordsObservacionesVecinos$value[iesLadoFrio] <= medioValores)
      # mediasACadaLado[j, 1] <- mean(coordsObservacionesVecinos$value[iesLadoCalido])
      # mediasACadaLado[j, 2] <- mean(coordsObservacionesVecinos$value[iesLadoFrio])
      # sdsACadaLado[j, 1] <- sd(coordsObservacionesVecinos$value[iesLadoCalido])
      # sdsACadaLado[j, 2] <- sd(coordsObservacionesVecinos$value[iesLadoFrio])
      
      #if (sum(iesLadoCalido) > 1 & sum(iesLadoFrio) > 1) {
      #  diferenciaEntreMediasDeCadaLado[j] <- mean(coordsObservacionesVecinos$value[iesLadoCalido]) - mean(coordsObservacionesVecinos$value[iesLadoFrio])
      #} else { diferenciaEntreMediasDeCadaLado[j] <- -Inf }
    }
    #j <- which(diferenciaEntreMediasDeCadaLado == max(diferenciaEntreMediasDeCadaLado))
    nDatosConSignoCorrectoACadaLado <- rowSums(nDatosConSignoCorrectoACadaLado)
    j <- which(nDatosConSignoCorrectoACadaLado == max(nDatosConSignoCorrectoACadaLado))[1]
    x0y0 <- c(nuevosXs[j], nuevosYs[j])
    
    intercept <- -mPerpendicular * (x0y0[1] * maxD + rX[1]) + (x0y0[2] * maxD + rY[1])
    # intercept <- x0y0[1]/m + x0y0[2]
    #intercept <- (x0y0[1] * (rX[2] - rX[1]) + rX[1])/m + (x0y0[2] * (rY[2] - rY[1]) + rY[1])
    
    if (FALSE) {
      shpBase <- cargarSHP('C:/mch/ArchivosProcesosLocales/CartografiaBase/uruguay_departamentos.shp')
      shpBase <- spTransform(shpBase, CRS(proj4string(coordsObservaciones)))
      xyLims <- getXYLims(c(coordsObservaciones, shpBase), ejesXYLatLong = F)
      grilla <- grillaSobreBoundingBox(shpBase)
      grilla <- SpatialGridDataFrame(grilla, data=data.frame(value=predict(modelo, newdata=data.frame(x=coordinates(grilla)[,1], y=coordinates(grilla)[,2]))))
      
      caja <- bbox(grilla)
      medioCoords <- rowMeans(caja)
      x1 <- medioCoords[1] - (caja[,2] - caja[,1]) * 0.05
      x2 <- medioCoords[1] + (caja[,2] - caja[,1]) * 0.05
      intercept = -m * medioCoords[1] + medioCoords[2]
      # intercept = -m * (medioCoords[1] * maxD + rX[1]) + medioCoords[2] * maxD + rY[1]
      
      # Grilla con el gradiente
      p <- mapearGrillaGGPlot(grilla, shpBase, continuo = T, dibujar = F, xyLims = xyLims,
                              titulo=format(x = fechasObservaciones[iFecha], format='%Y-%m-%d %H:%M'), 
                              coordsObservaciones = coordsObservaciones, dibujarPuntosObservaciones = F)
      p
      
      # Dirección del gradiente, sin ubicación
      df <- data.frame(x1 = x1, x2 = x2, y1 = m * x1 + intercept, y2 = m * x2 + intercept, value=NA)
      p <- mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpBase, xyLims = xyLims, zcol = 'value',
                              continuo = T, dibujarTexto = T, dibujar = F) + 
        geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df, arrow = arrow(length=unit(0.30,"cm")), size=1) +
        geom_abline(slope = mPerpendicular, intercept = -mPerpendicular * medioCoords[1] + medioCoords[2])
      p
      ggsave(p, file='D:/testsMCH/SRT/2-Gradiente2.png', dpi=90, width = 630 / 90, height = 630 / 90, units = 'in', type='cairo')
      
      ggsave(p, file='D:/testsMCH/SRT/3-DireccionGradiente.png', dpi=90, width = 630 / 90, height = 630 / 90, units = 'in', type='cairo')
      
      # Dirección perpendicular del gradiente, sin ubicación
      p <- mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpBase, xyLims = xyLims, zcol = 'value',
                              continuo = T, dibujarTexto = T, dibujar = F) + 
        geom_abline(slope = m, intercept = -m * medioCoords[1] + medioCoords[2])
      p
      ggsave(p, file='D:/testsMCH/SRT/4-PosicionGradiente1.png', dpi=90, width = 630 / 90, height = 630 / 90, units = 'in', type='cairo')    
      
      # x = (x1/m + y1 - c) / (m + 1/m)
      # Ubicaciones potenciales, proyección sobre perpendicular
      coordsOrig <- coordinates(coordsObservacionesVecinos)
      xsAux <- as.numeric((coordsOrig[, 1] * invM + coordsOrig[, 2] - intercept) / (m + invM))
      ysAux <- xsAux * m + intercept
      
      df2 <- data.frame(x1 = coordsOrig[, 1], x2 = xsAux, y1 = coordsOrig[, 2], y2 = ysAux, value=NA)
      
      p2 <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df2, arrow = arrow(length=unit(0.30,"cm")), size=0.5)
      p2
      ggsave(p2, file='D:/testsMCH/SRT/4-PosicionGradiente2.png', dpi=90, width = 630 / 90, height = 630 / 90, units = 'in', type='cairo')
      
      # Ubicaciones potenciales, puntos
      xsAux <- as.numeric(sort((coordsOrig[, 1] * invM + coordsOrig[, 2] - intercept) / (m + invM)))
      iAux <- 1:(length(xsAux)-1)
      # Tomo los medios de los intervalos entre los xs originales proyectados y calculo sus ys
      nuevosXs <- (xsAux[iAux] + xsAux[iAux+1]) * 0.5
      nuevosYs <- nuevosXs * m + intercept
      
      df3 <- data.frame(x=nuevosXs, y=nuevosYs, value=NA)
      p3 <- p + geom_point(data = df3)
      p3
      ggsave(p3, file='D:/testsMCH/SRT/4-PosicionGradiente3.png', dpi=90, width = 630 / 90, height = 630 / 90, units = 'in', type='cairo')
      
      # Gradiente final, ubicación y dirección
      i2 <- -mPerpendicular * (x0y0[1] * maxD + rX[1]) + (x0y0[2] * maxD + rY[1])
      # i2 <- -mPerpendicular * (x0y0[1]) + (x0y0[2])
      p4 <- mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpBase, xyLims=xyLims, zcol = 'value', continuo = T, dibujarTexto = T,
                               dibujar = F) + geom_abline(slope = mPerpendicular, intercept = i2)
      p4
      ggsave(p4, file='D:/testsMCH/SRT/4-PosicionGradiente4.png', dpi=90, width = 630 / 90, height = 630 / 90, units = 'in', type='cairo')
    }  
    # row.names(coordsObservaciones)[iPuntoATestear]
    
    coordsI <- as.numeric(coordinates(coordsObservaciones)[iPuntoATestear, ])
    coordsI[1] <- (coordsI[1] - rX[1]) / maxD
    coordsI[2] <- (coordsI[2] - rY[1]) / maxD
    
    if ((!creceConX & mPerpendicular >= 0) | (creceConX & mPerpendicular < 0)) {
      iesLadoCalido <- coords[, 2] >= mPerpendicular * (coords[, 1] - x0y0[1]) + x0y0[2]
      iesLadoFrio <- !iesLadoCalido
      iPuntoEnLadoCalido <- coordsI[2] >= mPerpendicular * (coordsI[1] - x0y0[1]) + x0y0[2]
    } else {
      iesLadoCalido <- coords[, 2] < mPerpendicular * (coords[, 1] - x0y0[1]) + x0y0[2]
      iesLadoFrio <- !iesLadoCalido
      iPuntoEnLadoCalido <- coordsI[2] < mPerpendicular * (coordsI[1] - x0y0[1]) + x0y0[2]
    }
    
    hayGradiente <- sum(iesLadoCalido) > 1 & sum(iesLadoFrio) > 1 & 
      mean(coordsObservacionesVecinos$value[iesLadoCalido]) - mean(coordsObservacionesVecinos$value[iesLadoFrio]) > umbralGradiente
    
    iIPuntoATestear <- which(row.names(coordsObservaciones)[iPuntoATestear] == names(iesLadoCalido))
    if (length(iIPuntoATestear) > 0) {
      iesLadoCalido <- iesLadoCalido[-iIPuntoATestear]
      iesLadoFrio <- iesLadoFrio[-iIPuntoATestear]
    }
    
    if (mapear) {
      if (hayGradiente) {
        print(mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpBase, zcol = 'value', continuo = T, dibujarTexto = T,
                                 dibujar = F) + geom_abline(slope = mPerpendicular, intercept = intercept))
      } else {
        print(mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpBase, zcol = 'value', continuo = T, dibujarTexto = T,
                                 dibujar = F))
      }
    }
  } else {
    hayGradiente <- FALSE
    mPerpendicular <- NA
    intercept <- NA
    iesLadoCalido <- integer(0)
    iesLadoFrio <- integer(0)
    iPuntoEnLadoCalido <- FALSE
  }
  
  return(list(hayGradiente=hayGradiente, pendiente=mPerpendicular, intercept=intercept, 
              iesLadoCalido=iesLadoCalido, iesLadoFrio=iesLadoFrio, 
              iPuntoEnLadoCalido=iPuntoEnLadoCalido))
}

deteccionGradienteEnMatrizDeObservaciones <- function(coordsObservaciones, fechasObservaciones, valoresObservaciones,
                                                      iColumnasATestear=seq.int(from = 1, to = ncol(valoresObservaciones), by = 1),
                                                      itsATestear=seq.int(from = 1, to = nrow(valoresObservaciones), by = 1),
                                                      maxDistKm=250, umbralAceptacionGradiente=5, archivoSalida=NULL, 
                                                      verbose=interactive() && nCoresAUsar == 1, nCoresAUsar=0) {
  deteccionGradienteEnMatrizDeObservacionesI <- function(i, coordsObservaciones, fechasObservaciones, valoresObservaciones, 
                                                         iesPuntosVecinos, umbralAceptacionGradiente, mapear, shpBase) {
    # i <- which(colnames(valoresObservaciones) == 'EMA.CANELA')
    # print(i)
    ij <- arrayInd(ind = i, .dim = dim(valoresObservaciones))
    # valoresObservaciones[1,i,drop=F]
    iesPuntosVecinosJ <- iesPuntosVecinos[[ij[1,2]]]
    
    coordsObservaciones$value <- valoresObservaciones[ij[1,1],]
    iNoNA <- !is.na(coordsObservaciones$value)
    aux <- coordsObservaciones[iNoNA,]
    
    if (length(aux) >= 4 && !is.na(valoresObservaciones[ij[1,1], ij[1,2]])) {
      gradiente <- deteccionGradienteEnPuntos(coordsObservaciones = coordsObservaciones, iPuntoATestear = ij[1,2],
                                              maxDist = maxDist, iesPuntosVecinos = iesPuntosVecinosJ,
                                              umbralGradiente = umbralAceptacionGradiente)
      return(gradiente$hayGradiente)
    } else { return(FALSE) }
  }
  
  maxDist <- distKmToP4Str(proj4string(coordsObservaciones), maxDistKm)
  iesPuntosVecinos <- getIVecinosAMenosDeMaxDist(coordsObservaciones = coordsObservaciones, maxDist = maxDist)
  # factorRMSE <- qnorm(1 - (1 - conf) / 2)
  
  if (nCoresAUsar <= 0) nCoresAUsar <- min(detectCores(T, T), ncol(valoresObservaciones))
  mapear=FALSE
  shpBase=NULL
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.qcTests'))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = source(paste(script.dir.qcTests, 'qcTests.r', sep='')))
    
    # structure(vapply(m, y, numeric(1)), dim=dim(m))
    hayGradiente <- parSapply(cl, X = 1:length(valoresObservaciones), FUN = deteccionGradienteEnMatrizDeObservacionesI, 
                              coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, 
                              valoresObservaciones=valoresObservaciones, iesPuntosVecinos=iesPuntosVecinos, 
                              umbralAceptacionGradiente=umbralAceptacionGradiente, mapear=mapear, shpBase=shpBase)
    
    stopCluster(cl)
  } else {
    hayGradiente <- vapply(X = 1:length(valoresObservaciones), FUN = deteccionGradienteEnMatrizDeObservacionesI, 
                           coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, 
                           valoresObservaciones=valoresObservaciones, iesPuntosVecinos=iesPuntosVecinos, 
                           umbralAceptacionGradiente=umbralAceptacionGradiente, mapear=mapear, shpBase=shpBase, 
                           FUN.VALUE = logical(1))
  }
  
  return(structure(hayGradiente, dim=dim(valoresObservaciones), dimnames=dimnames(valoresObservaciones)))
}

mapearResultadosDeteccionOutliers <- function(test, carpetaSalida=NULL, coordsObservaciones, valoresObservaciones, 
                                              shpBase, tamaniosPuntos=5, tiposOutliersDeInteres=tiposOutliersValoresSospechosos, 
                                              nCoresAUsar=0) {
  mapearResultadosDeteccionOutliersI <- function(i, test, carpetaSalida, coordsObservaciones, valoresObservaciones, 
                                                 shpBase, xyLims, tamaniosPuntos=5) {
    # i <- 1
    # i <- which(test$estacion=='Laguna del Sauce' & test$fecha=='2017-05-10 13:00')
    fecha <- test$fecha[i]
    iFecha <- which(row.names(valoresObservaciones) == fecha)
    if (!is.null(carpetaSalida)) {
      strFecha <- gsub(pattern = ' ', replacement = '_', x = gsub(pattern = '-|/|:', replacement = '', x = test$fecha[i]))
      archivoSalida <- paste(carpetaSalida, test$estacion[i], '_', strFecha, '.png', sep = '')
    } else { archivoSalida <- NULL }
    coordsObservaciones$value <- valoresObservaciones[iFecha, ]
    
    iOutliers <- which(test$fecha == fecha)
    
    if (length(iOutliers) > 1) {
      # Si hay más de un outlier los ordeno por latitud, longitud
      iEstaciones <- match(test$estacion[iOutliers], row.names(coordsObservaciones))
      coordsIEstaciones <- coordinates(coordsObservaciones)[iEstaciones,, drop=F]
      
      orden <- order(coordsIEstaciones[,2], coordsIEstaciones[,1], decreasing = c(TRUE, FALSE))
      iOutliers <- iOutliers[orden]
    }
    
    titulo <- paste(fecha, '\n', 
                    paste(test$estacion[iOutliers], ': Valor = ', test$valor[iOutliers], 
                          '. Estimado = ', round(test$estimado[iOutliers], digits = 1),
                          '. StdDif = ', round(test$stdDif[iOutliers], digits = 1), 
                          '. ', TTipoOutlierToString(test$tipoOutlier[iOutliers]), sep='', collapse = '\n'))
    
    titulo <- gsub(pattern = 'Estimado = NA. ', replacement = '', titulo, fixed = T)
    titulo <- gsub(pattern = 'StdDif = NA. ', replacement = '', titulo, fixed = T)
    
    iPuntosOutliers <- which(row.names(coordsObservaciones) %in% test$estacion[iOutliers]) 
    puntosAResaltar <- coordsObservaciones[iPuntosOutliers,]
    
    # nomArchResultados <- 'C:/testsMCH/EjemploSRT/1-Vecindad.png'
    # titulo <- paste(fecha, '\n', paste(test$estacion[iOutliers], ': Valor = ', test$valor[iOutliers], '.', sep=''))
    mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpBase, nomArchResultados = archivoSalida, dibujarTexto = T, 
                       zcol='value', titulo = titulo, xyLims = xyLims, dibujar = is.null(archivoSalida), tamaniosPuntos=tamaniosPuntos, 
                       contornearPuntos = nrow(valoresObservaciones) <= 50, puntosAResaltar=puntosAResaltar, widthPx = 1024, 
                       heightPx = 1024, continuo=T)
    return(NULL)
  }
  
  test <- test[test$tipoOutlier %in% tiposOutliersDeInteres, ]
  #shpBase <- shpMask$shp
  #carpetaSalida <- 'Resultados/QC/SRT/'
  test <- test[order(test$estacion, test$fecha), ]
  xyLims <- getXYLims(c(coordsObservaciones, shpBase), ejesXYLatLong = T)
  
  if (nCoresAUsar <= 0) nCoresAUsar <- min(detectCores(T, T), ncol(valoresObservaciones))
  
  if (nrow(test) > 0) {
    if (nCoresAUsar > 1 && !is.null(carpetaSalida)) {
      cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
      clusterExport(cl, varlist = c('script.dir.qcTests'))
      if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
      clusterEvalQ(cl = cl, expr = { 
        source(paste(script.dir.qcTests, 'qcTests.r', sep='')) 
        source(paste(script.dir.qcTests, '../interpolar/mapearEx.r', sep='')) 
      })
      parLapplyLB(cl = cl, seq.int(from = 1, to = nrow(test), by = 1), fun = mapearResultadosDeteccionOutliersI, test=test, carpetaSalida=carpetaSalida, 
                  coordsObservaciones=coordsObservaciones, valoresObservaciones=valoresObservaciones, shpBase=shpBase, xyLims=xyLims, tamaniosPuntos=tamaniosPuntos)
      stopCluster(cl)
    } else {
      lapply(seq.int(from = 1, to = nrow(test), by = 1), FUN = mapearResultadosDeteccionOutliersI, test=test, carpetaSalida=carpetaSalida, 
             coordsObservaciones=coordsObservaciones, valoresObservaciones=valoresObservaciones, shpBase=shpBase, xyLims=xyLims, tamaniosPuntos=tamaniosPuntos)
    }
  }
}

mapearResultadosDeteccionOutliersV2 <- function(test, carpetaSalida=NULL, coordsObservaciones, valoresObservaciones, 
                                                shpBase, tamaniosPuntos=5, tiposOutliersDeInteres=tiposOutliersValoresSospechosos, 
                                                nCoresAUsar=0) {
  mapearResultadosDeteccionOutliersV2I <- function(i, fechas, test, carpetaSalida, coordsObservaciones, valoresObservaciones, 
                                                   shpBase, xyLims, tamaniosPuntos, tiposOutliersDeInteres) {
    # i <- 1
    # i <- which(fechas=='2017-12-09 00:00')
    
    fecha <- fechas[i]
    if (!is.null(carpetaSalida)) {
      strFecha <- gsub(pattern = ' ', replacement = '_', x = gsub(pattern = '-|/|:', replacement = '', x = fecha))
      archivoSalida <- paste(carpetaSalida, strFecha, '.png', sep = '')
    } else { archivoSalida <- NULL }
    iFecha <- which(row.names(valoresObservaciones) == fecha)
    coordsObservaciones$value <- valoresObservaciones[iFecha, ]
    
    iOutliers <- which(test$fecha == fecha & test$tipoOutlier %in% tiposOutliersDeInteres)
    
    if (length(iOutliers) > 0) {
      if (length(iOutliers) > 1) {
        # Si hay más de un outlier los ordeno por latitud, longitud
        iEstaciones <- match(test$estacion[iOutliers], row.names(coordsObservaciones))
        coordsIEstaciones <- coordinates(coordsObservaciones)[iEstaciones,, drop=F]
        
        orden <- order(coordsIEstaciones[,2], coordsIEstaciones[,1], decreasing = c(TRUE, FALSE))
        iOutliers <- iOutliers[orden]
      } 
      
      titulo <- paste(fecha, '\n', 
                      paste(test$estacion[iOutliers], ': Valor = ', test$valor[iOutliers], 
                            '. Estimado = ', round(test$estimado[iOutliers], digits = 1),
                            '. StdDif = ', round(test$stdDif[iOutliers], digits = 1), 
                            '. ', TTipoOutlierToString(test$tipoOutlier[iOutliers]), sep='', collapse = '\n'))
      titulo <- gsub(pattern = 'Estimado = NA. ', replacement = '', titulo, fixed = T)
      titulo <- gsub(pattern = 'StdDif = NA. ', replacement = '', titulo, fixed = T)

      iPuntosOutliers <- which(row.names(coordsObservaciones) %in% test$estacion[iOutliers]) 
      puntosAResaltar <- coordsObservaciones[iPuntosOutliers,]
    } else { 
      titulo <- paste(fecha, '\nNo se encontraron valores sospechosos.')
      puntosAResaltar <- NULL
    }
    
    # nomArchResultados <- 'C:/testsMCH/EjemploSRT/1-Vecindad.png'
    # titulo <- paste(fecha, '\n', paste(test$estacion[iOutliers], ': Valor = ', test$valor[iOutliers], '.', sep=''))
    mapearPuntosGGPlot(puntos = coordsObservaciones, shpBase = shpBase, nomArchResultados = archivoSalida, dibujarTexto = T, 
                       zcol='value', titulo = titulo, xyLims = xyLims, dibujar = is.null(archivoSalida), tamaniosPuntos=tamaniosPuntos, 
                       contornearPuntos = nrow(valoresObservaciones) <= 50, puntosAResaltar=puntosAResaltar, widthPx = 1024, 
                       heightPx = 1024, continuo=T)
    return(NULL)
  }
  
  #test <- testSRT
  #shpBase <- shpMask$shp
  #carpetaSalida <- 'Resultados/QC/SRT/'
  test <- test[order(test$fecha, test$estacion), ]
  xyLims <- getXYLims(c(coordsObservaciones, shpBase), ejesXYLatLong = T)
  
  fechas <- unique(test$fecha)
  
  if (nCoresAUsar <= 0) nCoresAUsar <- min(detectCores(T, T), ncol(valoresObservaciones))
  
  if (nrow(test) > 0) {
    if (nCoresAUsar > 1 && !is.null(carpetaSalida)) {
      cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
      clusterExport(cl, varlist = c('script.dir.qcTests'))
      if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
      clusterEvalQ(cl = cl, expr = { 
        source(paste(script.dir.qcTests, 'qcTests.r', sep='')) 
        source(paste(script.dir.qcTests, '../interpolar/mapearEx.r', sep='')) 
      })
      parLapplyLB(cl = cl, seq_along(fechas), fun = mapearResultadosDeteccionOutliersV2I, fechas=fechas, test=test, carpetaSalida=carpetaSalida, 
                  coordsObservaciones=coordsObservaciones, valoresObservaciones=valoresObservaciones, shpBase=shpBase, xyLims=xyLims, 
                  tamaniosPuntos=tamaniosPuntos, tiposOutliersDeInteres=tiposOutliersDeInteres)
      stopCluster(cl)
    } else {
      lapply(seq_along(fechas), FUN = mapearResultadosDeteccionOutliersV2I, fechas=fechas, test=test, carpetaSalida=carpetaSalida, 
             coordsObservaciones=coordsObservaciones, valoresObservaciones=valoresObservaciones, shpBase=shpBase, xyLims=xyLims, 
             tamaniosPuntos=tamaniosPuntos, tiposOutliersDeInteres=tiposOutliersDeInteres)
    }
  }
}

ejecutarReemplazosSRT <- function(test, valoresObservaciones) {
  # Se efectua el reemplazo segun test$reemplazar
  # Los que tienen 0 se dejan como están
  # Los que tienen 1 se reemplazan por NA
  # Los que tienen 2 se reemplazan por el estimado de SRT
  # Los que tienen 3 se reemplazan por test$valorReemplazo
  # Los que tienen 4 se dejan como están pero además no volveran a ser chequeados en futuras ejecuciones de SRT
  iNAs <- which(test$reemplazar == 1)
  i <- iNAs[1]
  for (i in iNAs) valoresObservaciones[test$fecha[i], test$estacion[i]] <- NA
  
  iEstimados <- which(test$reemplazar == 2 && !is.na(test$estimado))
  for (i in iEstimados) valoresObservaciones[test$fecha[i], test$estacion[i]] <- test$estimado[i]
  
  iValorReemplazo <- which(test$reemplazar == 3 && !is.na(test$valorReemplazo))
  for (i in iValorReemplazo) valoresObservaciones[test$fecha[i], test$estacion[i]] <- test$valorReemplazo[i]
  
  return(valoresObservaciones)
}

obtenerValoresProtegidos <- function(test, valoresObservaciones) {
  iProtegidos <- which(test$reemplazar == 4)
  res <- matrix(FALSE, nrow = nrow(valoresObservaciones), ncol = ncol(valoresObservaciones))
  i <- 1
  for (i in seq_along(iProtegidos)) {
    iProtegido <- iProtegidos[i]
    
    iEstacion <- which(colnames(valoresObservaciones) == test$estacion[iProtegido])
    iFecha <- which(rownames(valoresObservaciones) == test$fecha[iProtegido])
    res[iFecha, iEstacion] <- TRUE
  }
  return(res)
}

getIVecinosAMenosDeMaxDist <- function(coordsObservaciones, maxDist=50, filtrarDistanciaCero=FALSE) {
  # Retorna para cada punto i en coordsObservaciones una lista con los índices de los demás puntos de coordsObservaciones
  # que estén a menos de maxDist de i. Las unidades de distancia son las especificadas en la 
  # proj4string de coordsObservaciones
  dists <- spDists(coordsObservaciones, longlat = FALSE)
  # x <- dists[277,]
  ies <- apply(dists, 1, FUN = function(x, maxStn, maxDist, elvDiff_m, filtrarDistanciaCero) {
    if (filtrarDistanciaCero) { ies <- which(x > 0 & x <= maxDist)
    } else { ies <- which(x <= maxDist) }
    
    ies <- ies[order(x[ies])]
    return(ies)
  }, maxStn=maxStn, maxDist=maxDist, elvDiff_m=elvDiff_m, filtrarDistanciaCero=filtrarDistanciaCero)
  if (!filtrarDistanciaCero) {
    ies <- lapply(seq_along(ies), FUN = function(x, ies) { return(setdiff(ies[[x]], x)) }, ies=ies)  
  }
  return(ies)
}

getICuadrantes <- function(iPunto, iesVecinos, coords) {
  derecha <- coords[iesVecinos[[iPunto]], 1] > coords[iPunto, 1]
  arriba <- coords[iesVecinos[[iPunto]], 2] > coords[iPunto, 2]
  return(derecha * 2 + arriba + 1)
}

getICuadrante <- function(iVecino, iPunto, coords) {
  #coords[c(iVecino, iPunto), ]
  if (coords[iVecino, 1] <= coords[iPunto, 1]) { # El vecino está a la izquierda
    if (coords[iVecino, 2] <= coords[iPunto, 2]) { return(1L) # El vecino está a la izquierda, abajo
    } else { return(2L) } # El vecino está a la izquierda, arriba
  } else {
    # El vecino está a la derecha
    if (coords[iVecino, 2] <= coords[iPunto, 2]) { return(3L) # El vecino está a la derecha, abajo
    } else { return(4L) } # El vecino está a la derecha, arriba
  }
}

spatialRegressionTest <- function(coordsObservaciones, fechasObservaciones, valoresObservaciones,
                                  iColumnasATestear=seq.int(from = 1, to = ncol(valoresObservaciones), by = 1),
                                  itsATestear=seq.int(from = 1, to = nrow(valoresObservaciones), by = 1),
                                  ventana=min(30, nrow(valoresObservaciones) %/% 2), minStn=2, maxStn=5, 
                                  maxDistKm=250, elvDiff_m=200, zcolElv='Altitud', minValAbs=-40, maxValAbs=60, 
                                  minAdjR2=0.5, factorRMSE=2.5, usarDeteccionDeGradientes=FALSE, archivoSalida=NULL, 
                                  verbose=interactive() && nCoresAUsar == 1, nCoresAUsar=0, 
                                  datosProtegidos = matrix(data = FALSE, nrow = nrow(valoresObservaciones), ncol = ncol(valoresObservaciones))) {
  spatialRegressionTestI <- function(i, ies, itsATestear, coordsObservaciones, fechasObservaciones, 
                                     dfValoresObservaciones, ventana=30, minStn=2, maxStn=5, maxDist=250, 
                                     elvDiff_m=200, zcolElv='Altitud', minValAbs=-40, maxValAbs=60, 
                                     minAdjR2=0.5, factorRMSE=2.5, usarDeteccionDeGradientes=FALSE, 
                                     verbose=FALSE, iNoNAs, datosProtegidos) {
    # i <- iColumnasATestear[1]
    # i <- which(colnames(valoresObservaciones) == 'Aeropuerto.Carrasco')
    # i <- which(colnames(valoresObservaciones) == 'Aeropuerto.Melilla')
    # i <- which(colnames(valoresObservaciones) == 'Colonia')
    # i <- which(colnames(valoresObservaciones) == 'Durazno')
    # i <- which(colnames(valoresObservaciones) == 'Laguna.del.Sauce')
    # i <- which(colnames(valoresObservaciones) == 'Melo')
    # i <- which(colnames(valoresObservaciones) == 'Paso.de.los.Toros')
    # i <- which(colnames(valoresObservaciones) == 'Prado')
    # i <- which(colnames(valoresObservaciones) == 'Treinta.y.Tres')
    # iFecha <- itsATestear[1]
    # iFecha <- iFecha + 1
    # iFecha <- which(row.names(valoresObservaciones) == '2017-05-01 20:00')
    # iFecha <- which(row.names(valoresObservaciones) == '2014-01-21 15:00')
    # iFecha <- which(row.names(valoresObservaciones) == '2014-01-02 14:00')
    # iFecha <- which(row.names(valoresObservaciones) == '2017-05-05 15:00')
    # iFecha <- which(row.names(valoresObservaciones) == '2017-09-20 20:00')
    # iFecha <- which(row.names(valoresObservaciones) == '2017-09-02 07:00')
    # iFecha <- which(!is.na(valoresObservaciones[,i]))[1]
    # iFecha <- which(tiposOutliers %in% c(1,2))[1]
    # iFecha <- 1
    
    n <- length(itsATestear)
    estimados <- numeric(n)
    tiposOutliers <- integer(n)
    stdDifs <- numeric(n)
    # ies[[i]] <- base::setdiff(1:length(coordsObservaciones), i)
    n <- 1
    for (iFecha in itsATestear) {
      val <- dfValoresObservaciones[iFecha, i]
      if (verbose) print(paste(i, '.', colnames(dfValoresObservaciones)[i], ', ', iFecha, '. ', row.names(dfValoresObservaciones)[iFecha], ', ', val, sep=''))
      
      tsVentana <- getVentana(ti = iFecha, nT = length(fechasObservaciones), tamanioSemiVentana = ventana)
      # Saco las fechas que son NA en la estación que se está chequeando
      tsAUsarVentana <- tsVentana$tsVentana[iNoNAs[tsVentana$tsVentana, i]]
      rm(tsVentana)
      
      if (length(tsAUsarVentana) >= ventana) {
        # La fecha en cuestión la dejo siempre, aunque sea NA, para poder hacer el estimado
        tsAUsarVentana <- base::union(tsAUsarVentana, iFecha)
        iTiEnTsVentana <- which(iFecha == tsAUsarVentana)
        # Idea para considerar el ciclo diario, no dio buenos resultados, empeoró los RMSE y Adj R^2
        # tsAUsarVentana <- tsAUsarVentana[tsAUsarVentana %% 24 == tsAUsarVentana[iTiEnTsVentana] %% 24]
        # iTiEnTsVentana <- which(iFecha == tsAUsarVentana)
        
        # Me fijo cuantas fechas tienen disponibles los vecinos en la ventana, 
        # exijo al menos (semi)ventana o 100 fechas disponibles y disponer del dato de iFecha
        bVecinoConDatosSuficientes <- iNoNAs[tsAUsarVentana[iTiEnTsVentana], ies[[i]]]
        iesConDatosSuficientes <- ies[[i]][bVecinoConDatosSuficientes]
        # Hasta acá iesConDatosSuficientes tiene los que no son NA en iFecha
        
        if (length(iesConDatosSuficientes) >= 4 & usarDeteccionDeGradientes) {
          # Solo se considerará gradiente si al menos hay 2 puntos a cada lado del gradiente, así que al menos
          # debe haber 4 vecinos para calcular el gradiente
          coordsObservaciones@data[,'value'] <- as.numeric(dfValoresObservaciones[iFecha, ])
          
          gradiente <- deteccionGradienteEnPuntos(coordsObservaciones, iPuntoATestear = i, maxDist = maxDist, 
                                                  iesPuntosVecinos = c(iesConDatosSuficientes, i),
                                                  mapear=F, shpBase=shpBase, umbralGradiente = 4)
          
          if (gradiente$hayGradiente) {
            if (gradiente$iPuntoEnLadoCalido) { iesConDatosSuficientes <- iesConDatosSuficientes[gradiente$iesLadoCalido]
            } else { iesConDatosSuficientes <- iesConDatosSuficientes[gradiente$iesLadoFrio] }
          }
        }
        
        bVecinoConDatosSuficientes <- colSums(iNoNAs[tsAUsarVentana, iesConDatosSuficientes, drop=F]) >= max(30, ventana / 3)
        iesConDatosSuficientes <- iesConDatosSuficientes[bVecinoConDatosSuficientes]
        rm(bVecinoConDatosSuficientes)
        
        # plot(iNoNAs[tsAUsarVentana, iesConDatosSuficientes[[1]]])
        
        if (length(iesConDatosSuficientes) >= minStn) {
          # x <- iesConDatosSuficientes[1]
          lms <- lapply(iesConDatosSuficientes, FUN = function(x) {
            return(lm(formula = as.formula(paste(paste(colnames(dfValoresObservaciones)[c(i, x)], collapse = '~'), '+1', sep='')), 
                      data = dfValoresObservaciones, subset = tsAUsarVentana[iNoNAs[tsAUsarVentana, x]]))
          })
          
          # res <- round(do.call(rbind, lapply(lms, FUN = function(x) { return(data.frame(a=coefficients(x)[2], b=coefficients(x)[1], adjR2=summary(x)$adj.r.squared, rmse=mean(x$residuals^2))) })), 2)
          # res[order(res$rmse),]
          # dim(res)
          
          adjR2s <- sapply(lms, FUN = function(x) { return(summary(x)$adj.r.squared) })
          # cbind(colnames(dfValoresObservaciones)[iesConDatosSuficientes], adjR2s)
          iModelos <- adjR2s >= minAdjR2 & adjR2s < 1 & !is.nan(adjR2s)
        } else {
          iModelos <- FALSE
        }
        
        if (sum(iModelos) >= minStn) {
          lms <- lms[iModelos]
          iesConDatosSuficientes <- iesConDatosSuficientes[iModelos]
          mses <- sapply(lms, FUN = function(x) { return(mean(x$residuals^2)) })
          # cbind(colnames(dfValoresObservaciones)[iesConDatosSuficientes], mses)[order(mses),]
          
          if (length(mses) > maxStn) {
            iModelos <- order(mses, decreasing = FALSE)[seq.int(from = 1, to = maxStn, by = 1)]
            # colnames(dfValoresObservaciones[,iesConDatosSuficientes])[iModelos]
            # sort(mses)
            # which(mses <= mean(sort(mses, decreasing = FALSE)[1:3]) + 3 * sd(sort(mses, decreasing = FALSE)[1:3]))
            # which(mses <= median(sort(mses, decreasing = FALSE)[1:3]) + 3 * mad(sort(mses, decreasing = FALSE)[1:3]))
            
            lms <- lms[iModelos]
            mses <- mses[iModelos]
            iesConDatosSuficientes <- iesConDatosSuficientes[iModelos]
          }
          
          #if (length(mses) > 3) {
          #  mejoresMSEs <- sort(mses, decreasing = FALSE)[1:3]
          #  iModelos <- which(mses <= mean(mejoresMSEs) + 2.5 * sd(mejoresMSEs))
          #  #which(mses <= median(sort(mses, decreasing = FALSE)[1:3]) + 3 * mad(sort(mses, decreasing = FALSE)[1:3]))
          #  
          #  lms <- lms[iModelos]
          #  mses <- mses[iModelos]
          #  iesConDatosSuficientes <- iesConDatosSuficientes[iModelos]
          #}
          
          #dfValoresObservaciones[tsAUsarVentana[iTiEnTsVentana], iesConDatosSuficientes, drop=F]
          fitValues <- sapply(X = seq_along(lms), FUN = function(x) { 
            predict(object = lms[[x]], 
                    newdata=dfValoresObservaciones[tsAUsarVentana[iTiEnTsVentana], iesConDatosSuficientes[x], drop=F])
          })
          
          invMSEs <- 1 / mses
          pesos <- invMSEs / sum(invMSEs)
          
          # cbind(t(dfValoresObservaciones[tsAUsarVentana[iTiEnTsVentana], iesConDatosSuficientes]), fitValues, pesos)
          # round(cbind(t(dfValoresObservaciones[tsAUsarVentana[iTiEnTsVentana], iesConDatosSuficientes]), fitValues, pesos),2)
          
          # Versión Original paper 2005
          # Hago los cuadrados preservando el signo
          # estimado <- sum((fitValues^2) * sign(fitValues) * pesos)
          # estimado <- sqrt(abs(estimado)) * sign(estimado)
          
          # Versión del paper 2012
          estimado <- sum(fitValues * pesos)
          
          #res <- data.frame(estacion=character(n), fecha=character(n), valor=numeric(n), estimado=numeric(n), 
          #                  tipoOutlier=integer(n), stdDif=numeric(n), reemplazar=integer(n),
          #                  valorReemplazo=numeric(n), stringsAsFactors = F)
          if (is.na(val)) {
            # El dato en la estacion i, fecha iFecha es NA
            estimados[n] <- estimado
            tiposOutliers[n] <- TTO_ValorNulo
            stdDifs[n] <- NA_real_
          } else  {
            # Versión original de Hubbard 2005 y 2012
            meanInvMSEs <- mean(invMSEs)
            
            # Versión revisada, en vez de hacer el promedio simple hace la suma ponderada de los MSEs
            # meanInvMSEs <- sum(invMSEs * pesos)
            
            # Esta otra trata a los residuos como distribuciones normales y a su MSE como la varianza.
            # La varianza de la combinación lineal de normales es la sumatoria de las varianzas por los 
            # pesos de la combinación al cuadrado
            # s <- sum(sqrt(mses * pesos^2))
            s <- 1 / sqrt(meanInvMSEs)
            
            stdDif <- (val - estimado) / s
            if (estimado - factorRMSE * s > val) { tipoOutlier <- TTO_OutlierPorLoBajo
            } else if (estimado + factorRMSE * s < val) { tipoOutlier <- TTO_OutlierPorLoAlto
            } else { tipoOutlier <- TTO_SinProblemasDetectados }
            
            estimados[n] <- estimado
            tiposOutliers[n] <- tipoOutlier
            stdDifs[n] <- stdDif
          }
        } else {
          # No hay al menos minStn vecinos con AdjR^2 > minAdjR2 y con al menos ventana valores disponibles en el período de la ventana
          estimados[n] <- NA_real_
          tiposOutliers[n] <- TTO_SinDatosSuficientesEnVecinos
          stdDifs[n] <- NA_real_
        }
      } else {
        # La estación no tiene al menos ventana valores disponibles en el período de la ventana
        estimados[n] <- NA_real_
        tiposOutliers[n] <- TTO_SinDatosSuficientesEnLaEstacion
        stdDifs[n] <- NA_real_      
      }
      
      if (datosProtegidos[iFecha, i]) tiposOutliers[n] <- TTO_ValorConfirmado
      n <- n + 1
    }
    
    #res <- data.frame(estacion=character(n), fecha=character(n), valor=numeric(n), estimado=numeric(n), 
    #                  tipoOutlier=integer(n), stdDif=numeric(n), reemplazar=integer(n),
    #                  valorReemplazo=numeric(n), stringsAsFactors = F)
    return(data.frame(estacion=row.names(coordsObservaciones)[i], fecha=row.names(dfValoresObservaciones)[itsATestear], 
                      valor=dfValoresObservaciones[itsATestear, i], estimado=estimados, tipoOutlier=tiposOutliers, 
                      stdDif=stdDifs, reemplazar=0L, valorReemplazo=NA_real_, stringsAsFactors = F))
  }
    
  maxDist <- distKmToP4Str(proj4string(coordsObservaciones), maxDistKm)
  ies <- getIVecinosAMenosDeMaxDist(coordsObservaciones = coordsObservaciones, maxDist = maxDist)
  
  dfValoresObservaciones <- as.data.frame(valoresObservaciones)
  dfValoresObservaciones[dfValoresObservaciones < minValAbs | dfValoresObservaciones > maxValAbs] <- NA
  iNoNAs <- !is.na(dfValoresObservaciones)
  
  # factorRMSE <- qnorm(1 - (1 - conf) / 2)
  
  if (nCoresAUsar <= 0) nCoresAUsar <- min(detectCores(T, T), ncol(valoresObservaciones))
  
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.qcTests'))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = source(paste(script.dir.qcTests, 'qcTests.r', sep='')))
    test <- parLapplyLB(cl = cl, X=iColumnasATestear, fun = spatialRegressionTestI, 
                        ies=ies, itsATestear=itsATestear, coordsObservaciones=coordsObservaciones, 
                        fechasObservaciones=fechasObservaciones, dfValoresObservaciones=dfValoresObservaciones, 
                        ventana=ventana, minStn=minStn, maxStn=maxStn, maxDist=maxDist, elvDiff_m=elvDiff_m, 
                        zcolElv=zcolElv, minAdjR2=minAdjR2, factorRMSE=factorRMSE, 
                        usarDeteccionDeGradientes=usarDeteccionDeGradientes, verbose=verbose, iNoNAs=iNoNAs, 
                        datosProtegidos=datosProtegidos)
    stopCluster(cl)
  } else {
    # verbose <- T
    test <- lapply(X=iColumnasATestear, FUN = spatialRegressionTestI, 
                   ies=ies, itsATestear=itsATestear, coordsObservaciones=coordsObservaciones, 
                   fechasObservaciones=fechasObservaciones, dfValoresObservaciones=dfValoresObservaciones, 
                   ventana=ventana, minStn=minStn, maxStn=maxStn, maxDist=maxDist, elvDiff_m=elvDiff_m, 
                   zcolElv=zcolElv, minAdjR2=minAdjR2, factorRMSE=factorRMSE, 
                   usarDeteccionDeGradientes=usarDeteccionDeGradientes, verbose=verbose, iNoNAs=iNoNAs, 
                   datosProtegidos=datosProtegidos)
  }
  
  # test[[1]][215,]
  # Código para profiling
  #system.time({
  #Rprof()
  #lele <- spatialRegressionTestI(i = 5,
  #                               ies=ies, coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, 
  #                               dfValoresObservaciones=dfValoresObservaciones, ventana=ventana, minStn=minStn, maxStn=maxStn, 
  #                               maxDist=maxDist, elvDiff_m=elvDiff_m, zcolElv=zcolElv, minAdjR2=minAdjR2, f=f, iNoNAs = iNoNAs)
  #Rprof(NULL)
  #})
  #summaryRprof()
  
  test <- do.call(rbind, test)
  row.names(test) <- paste(make.names(test$estacion), test$fecha, sep='_')
  
  if (!is.na(minValAbs)) {
    iBajos <- which(valoresObservaciones[itsATestear, iColumnasATestear] < minValAbs)
    test$valor[iBajos] <- valoresObservaciones[itsATestear, iColumnasATestear][iBajos]
    test$tipoOutlier[iBajos] <- TTO_OutlierPorLoBajo
  }
  if (!is.na(maxValAbs)) {
    iAltos <- which(valoresObservaciones[itsATestear, iColumnasATestear] > maxValAbs)
    test$valor[iAltos] <- valoresObservaciones[itsATestear, iColumnasATestear][iAltos]
    test$tipoOutlier[iAltos] <- TTO_OutlierPorLoAlto
  }
  
  if (!is.null(archivoSalida)) {
    dir.create(dirname(archivoSalida), showWarnings = F, recursive = T)
    write.table(test, archivoSalida, sep = '\t', col.names = T, row.names = T)
  }
  # test <- read.table('Resultados/srt.tsv', header = T, sep = '\t', row.names = 1, stringsAsFactors = F)
  # test[test$tipoOutlier %in% c(1,2),]
  # test[test$tipoOutlier %in% c(-3),]
  # test[test$tipoOutlier %in% c(-1),]
  # sum(test$tipoOutlier %in% c(1,2))
  # sum(test$tipoOutlier %in% c(-3))
  # sum(test$tipoOutlier %in% c(-2))
  # sum(test$tipoOutlier %in% c(-1))
  
  # test <- test[test$tipoOutlier %in% tiposOutliersDeInteres,]
  return(test)
}

testEspacialPrecipitacionI <- function(i, iesVecinos, coordsObservaciones, fechasObservaciones, dfValoresObservaciones,
                                       iNoNAs,
                                       itsATestear, ispMax, ispObs, isdMin, isdObs, isdQ1, fInf, fSup,
                                       maxDist, minNCuadrantes, minNVecinosPorCuadrante, datosProtegidos, verbose) {
  # i <- 102
  # iFecha <- 56
  # i <- which(colnames(dfValoresObservaciones) == 'Moirones')
  # iFecha <- which(row.names(dfValoresObservaciones) == '2016-08-05 00:00')
  # mapearPuntosGGPlot(SpatialPointsDataFrame(geometry(coordsObservaciones), data = data.frame(value=as.numeric(dfValoresObservaciones[iFecha,]))), shpBase = shpBase, continuo = T, dibujar = F)  
  iesVecinosI <- iesVecinos[[i]]
  estimacionPobs <- apply(dfValoresObservaciones[, iesVecinosI], MARGIN = 1, FUN = median, na.rm=T)
    
  nCuantiles <- 5
  cuantilesPositivos <- c(0, quantile(estimacionPobs[estimacionPobs > 0], 
                                      probs=seq(from=0, to=1, length.out = nCuantiles), na.rm=T))
    
  rango <- range(estimacionPobs[estimacionPobs > 0], na.rm=T)
  nClasif <- 6
  clasif <- c(0, seq(from=rango[1] + diff(rango) * 0.05, to= rango[1] + diff(rango) * 0.5, length.out = nClasif-1))
  clasifAux <- c(clasif[2:length(clasif)], Inf)
  clasif <- clasif[1:(length(clasif)-1)]
    
  limites <- matrix(data = NA, nrow = nClasif, ncol = 3)
    
  limitesAdjTukey <- limites
  limitesTukey <- limites
  limitesMedianRule <- limites
  limitesAdjMedian <- limites
  
  valoresVecinos <- as.numeric(na.omit(unlist(dfValoresObservaciones[, iesVecinosI])))
  medCouple <- mc(valoresVecinos[valoresVecinos > 0])
  
  mStats <- matrix(nrow = length(clasif), ncol = 6)
  colnames(mStats) <- c('mediana', 'Q1', 'Q3', 'mc', 'IQR', 'n')
  
  iClasif <- 4
  for (iClasif in 1:length(clasif)) {
    iFechas <- between(x = estimacionPobs, minInclusivo = clasif[iClasif], maxExclusivo = clasifAux[iClasif])
    # sum(iFechas, na.rm = T)
    valoresClasif <- as.numeric(na.omit(unlist(dfValoresObservaciones[iFechas, iesVecinosI])))
    valoresClasif <- valoresClasif[valoresClasif > 0]
    # plot(density(valoresClasif))
    # hist(valoresClasif)
    # sort(unique(valoresClasif))
    mediana <- median(valoresClasif)
    cuartiles <- quantile(valoresClasif, probs=c(0.25, 0.75))
    # medCouple <- mc(valoresClasif)
    rangoIntercuartil <- cuartiles[2] - cuartiles[1]
    if (iClasif > 1) rangoIntercuartil <- max(rangoIntercuartil, mStats[iClasif - 1, 5])
    if (medCouple >= 0) { coefs <- c(-4, 3)
    } else { coefs <- c(-3, 4) }
    
    # hist(log(valoresVecinos))
    mStats[iClasif, ] <- c(mediana, cuartiles, medCouple, rangoIntercuartil, length(valoresClasif))
    
    limitesAdjTukey[iClasif, ] <- c(clasif[iClasif], cuartiles[1] - fInf * exp(coefs[1] * medCouple) * rangoIntercuartil, cuartiles[2] + fSup * exp(coefs[2] * medCouple) * rangoIntercuartil)
    limitesTukey[iClasif, ] <- c(clasif[iClasif], cuartiles[1] - fInf * rangoIntercuartil, cuartiles[2] + fSup * rangoIntercuartil)
    limitesMedianRule[iClasif, ] <- c(clasif[iClasif], mediana - fInf * rangoIntercuartil, mediana + fSup * rangoIntercuartil)
    limitesAdjMedian[iClasif, ] <- c(clasif[iClasif], mediana - fInf * exp(coefs[1] * medCouple) * rangoIntercuartil, mediana + fSup * exp(coefs[2] * medCouple) * rangoIntercuartil)
  }
  
  n <- length(itsATestear)
  tiposOutliers <- integer(n)
  # stdDifs <- numeric(n)
  # ies[[i]] <- base::setdiff(1:length(coordsObservaciones), i)
  n <- 1
  for (iFecha in itsATestear) {
    Pobs <- dfValoresObservaciones[iFecha, i]
    # verbose <- TRUE
    if (verbose) print(paste(i, '.', colnames(dfValoresObservaciones)[i], ', ', iFecha, '. ', row.names(dfValoresObservaciones)[iFecha], ', ', Pobs, sep=''))

    if (!datosProtegidos[iFecha, i]) {
      if (!is.na(Pobs)) {
        iNoNAsVecinos <- which(!is.na(dfValoresObservaciones[iFecha, iesVecinosI]))
        cuadrantesVecinosNoNA <- cuadrantesVecinos[iNoNAsVecinos]
        nVecinosPorCuadrante <- rep(0L, minNCuadrantes)
        for (iCuadrante in seq_along(cuadrantesVecinosNoNA)) nVecinosPorCuadrante[cuadrantesVecinosNoNA[iCuadrante]] <- nVecinosPorCuadrante[cuadrantesVecinosNoNA[iCuadrante]] + 1L

        # Si hay algún vecino no NA y la cantidad de vecinos por cuadrante es mayor a minNVecinosPorCuadrante en el cuadrante con menos
        # valores disponibles
        if (min(nVecinosPorCuadrante) >= minNVecinosPorCuadrante) {
          valoresVecinos <- as.numeric(dfValoresObservaciones[iFecha, iesVecinosI[iNoNAsVecinos]])
          
          if (Pobs >= ispObs && max(valoresVecinos) < ispMax) {
            tiposOutliers[iFecha] <- TTO_PrecipitacionAislada
            stdDifs[iFecha] <- NA_real_
          } else if (Pobs < isdObs && min(valoresVecinos) >= isdMin && quantile(valoresVecinos, probs=0.25) >= isdQ1) {
            tiposOutliers[iFecha] <- TTO_SequedadAislada
            stdDifs[iFecha] <- NA_real_
          } else {
            # simpleInvDistanceWeighting(iFecha, data = dfValoresObservaciones[, iesVecinosI], dists = distsVecinos)
            # valoresVecinos <- as.numeric(unlist(dfValoresObservaciones[iFechas, iesVecinosI[iNoNAsVecinos]]))
            iClasif <- between(median(valoresVecinos), clasif, clasifAux)
            # iClasif <- between(simpleInvDistanceWeighting(iFecha, data = dfValoresObservaciones[, iesVecinosI], dists = distsVecinos), clasif, clasifAux)
            limitesI <- limitesTukey[iClasif]
            
            #cuartiles <- quantile(valoresVecinos, probs=c(0.25, 0.75), na.rm = T)
            #medCouple <- mc(valoresVecinos, na.rm = T)
            
            #if (medCouple >= 0) { coefs <- c(-3.5, 4)
            #} else { coefs <- c(-4, 3.5) }
            #limites <- c(cuartiles[1] - f * exp(coefs[1] * medCouple) * rangoIntercuartil, 
            #             cuartiles[2] + f * exp(coefs[2] * medCouple) * rangoIntercuartil)
                        
            #if (F) {
            #  mediana <- median(valoresVecinos, na.rm = T)
            #  amplitud <- diff(range(valoresVecinos, na.rm = T))
            #  
            #  rangoIntercuartil <- IQR(valoresVecinos, na.rm = T)
            #  
            #  f <- 1.5
            #  c(mediana - (amplitud / 2 + 30), mediana + 2 * (amplitud / 2 + 30))
            #  c(mediana - amplitud * f, mediana + 2 * f * amplitud)
            #  c(cuartiles[1] - f * rangoIntercuartil, cuartiles[2] + f * rangoIntercuartil)
            #  c(mediana - f * rangoIntercuartil, mediana + f * rangoIntercuartil)
            #}
            
            #desv <- (Pobs - mediana) / rangoIntercuartil
            if (Pobs < limitesI[2]) { tiposOutliers[iFecha] <- TTO_OutlierPorLoBajo
            } else if (Pobs > limitesI[3]) { tiposOutliers[iFecha] <- TTO_OutlierPorLoAlto
            } else { tiposOutliers[iFecha] <- TTO_SinProblemasDetectados }
            stdDifs[iFecha] <- NA_real_
          }
        } else {
          tiposOutliers[iFecha] <- TTO_SinDatosSuficientesEnVecinos
          stdDifs[iFecha] <- NA_real_
        }
      } else {
        tiposOutliers[iFecha] <- TTO_ValorNulo
        stdDifs[iFecha] <- NA_real_
      }
    } else {
      tiposOutliers[iFecha] <- TTO_ValorConfirmado
      stdDifs[iFecha] <- NA_real_
    }
  }
  return(data.frame(estacion=row.names(coordsObservaciones)[i], fecha=row.names(dfValoresObservaciones)[itsATestear], 
                    valor=dfValoresObservaciones[itsATestear, i], estimado=NA_real_, tipoOutlier=tiposOutliers, 
                    stdDif=stdDifs, reemplazar=0L, valorReemplazo=NA_real_, stringsAsFactors = F))
}

testEspacialPrecipitacionIV2 <- function(
    i, iesVecinos, coordsObservaciones, fechasObservaciones, dfValoresObservaciones, iNoNAs, 
    itsATestear, ispMax, ispObs, isdMin, isdObs, isdEstMin, fInf, fSup, amplitudMin, maxDist, 
    minNCuadrantes, minNVecinosPorCuadrante, datosProtegidos, verbose) {
  # i <- 1
  # iFecha <- 31
  # i <- which(colnames(dfValoresObservaciones) == 'Villa.25.de.Mayo')
  # iFecha <- which(row.names(dfValoresObservaciones) == '2017-10-20 00:00')
  # mapearPuntosGGPlot(SpatialPointsDataFrame(geometry(coordsObservaciones), data = data.frame(value=as.numeric(dfValoresObservaciones[iFecha,]))), shpBase = shpBase, continuo = T, dibujar = F, dibujarTexto = T, tamaniosPuntos = 2)  
  # print(i)
  iesVecinosI <- iesVecinos[[i]]
  if (length(iesVecinosI) >= minNCuadrantes) {
    valoresVecinos <- dfValoresObservaciones[, iesVecinosI, drop=F]
    iNoNAsVecinos <- iNoNAs[, iesVecinosI, drop=F]
    
    cuadrantesVecinos <- getICuadrantes(iesVecinos = iesVecinos, iPunto = i, coords = coordinates(coordsObservaciones))
    nVecinosPorCuadrante <- matrix(data = 0L, nrow = nrow(dfValoresObservaciones), ncol = 4)
    for (iFechaAux in 1:nrow(dfValoresObservaciones)) {
      cuadrantesVecinosI <- cuadrantesVecinos[iNoNAsVecinos[iFechaAux,]]
      for (iCuadrante in seq_along(cuadrantesVecinosI))
        nVecinosPorCuadrante[iFechaAux, cuadrantesVecinosI[iCuadrante]] <- nVecinosPorCuadrante[iFechaAux, cuadrantesVecinosI[iCuadrante]] + 1L
    }
    iFechasConVecinosSuficientes <- rowSums(nVecinosPorCuadrante >= minNVecinosPorCuadrante) >= minNCuadrantes
    #iFechasConVecinosSuficientes <- rowMins(nVecinosPorCuadrante) >= minNVecinosPorCuadrante
    
    coords <- coordinates(coordsObservaciones)[c(i, iesVecinosI),]
    distsVecinos <- spDistsN1(pts = coords[2:nrow(coords), ,drop=F], pt = coords[1, ,drop=F])
    
    estimado <- simpleInvDistanceWeighting(data = valoresVecinos, dists=distsVecinos)
    #estimado <- sapply(X = 1:nrow(dfValoresObservaciones), FUN = simpleInvDistanceWeighting, data = valoresVecinos, dists=distsVecinos)
    #estimadoMedian <- apply(X = valoresVecinos, MARGIN = 1, FUN = median, na.rm=T)
    estimado[!iFechasConVecinosSuficientes] <- NA
    #estimadoMedian[!iFechasConVecinosSuficientes] <- NA
    #estimado <- numeric(length = length(estimadoIDW))
    
    n <- length(itsATestear)
    tiposOutliers <- integer(n)
    stdDifs <- numeric(n)
    for (iFecha in itsATestear) {
      Pobs <- dfValoresObservaciones[iFecha, i]
      # verbose <- TRUE
      if (verbose) print(paste(i, '.', colnames(dfValoresObservaciones)[i], ', ', iFecha, '. ', row.names(dfValoresObservaciones)[iFecha], ', ', Pobs, sep=''))
      
      if (!datosProtegidos[iFecha, i]) {
        if (!is.na(Pobs)) {
          if (iFechasConVecinosSuficientes[iFecha]) {
            valoresVecinosI <- as.numeric(valoresVecinos[iFecha, iNoNAsVecinos[iFecha,]])
            
            if (Pobs >= ispObs && max(valoresVecinosI) <= ispMax) {
              tiposOutliers[iFecha] <- TTO_PrecipitacionAislada
              stdDifs[iFecha] <- NA_real_
            } else if (Pobs <= isdObs && min(valoresVecinosI) >= isdMin && estimado[iFecha] >= isdEstMin) {
              tiposOutliers[iFecha] <- TTO_SequedadAislada
              stdDifs[iFecha] <- NA_real_
            } else {
              amplitud <- diff(range(valoresVecinosI))
              #amplitud <- iqr(valoresVecinosI)
              #amplitud <- diff(quantile(valoresVecinosI, probs=c(0.15, 0.85)))
              
              if (amplitud > amplitudMin) {
                #if (amplitud > 10) { estimado[iFecha] <- estimadoIDW[iFecha]
                #} else { estimado[iFecha] <- estimadoMedian[iFecha] }
                stdDifs[iFecha] <- (Pobs - estimado[iFecha]) / amplitud
                if (stdDifs[iFecha] < -fInf) { tiposOutliers[iFecha] <- TTO_OutlierPorLoBajo
                } else if (stdDifs[iFecha] > fSup) { tiposOutliers[iFecha] <- TTO_OutlierPorLoAlto
                } else { tiposOutliers[iFecha] <- TTO_SinProblemasDetectados }
              } else {
                # Si ninguno de los vecinos muestra precipitación y Pobs no fue precipitación aislada
                # asumimos que no hay problemas
                tiposOutliers[iFecha] <- TTO_SinProblemasDetectados 
              }
            }
          } else {
            tiposOutliers[iFecha] <- TTO_SinDatosSuficientesEnVecinos
            stdDifs[iFecha] <- NA_real_
          }
        } else { tiposOutliers[iFecha] <- TTO_ValorNulo }
      } else { tiposOutliers[iFecha] <- TTO_ValorConfirmado }
    }
    return(data.frame(estacion=row.names(coordsObservaciones)[i], fecha=row.names(dfValoresObservaciones)[itsATestear], 
                      valor=dfValoresObservaciones[itsATestear, i], estimado=estimado, tipoOutlier=tiposOutliers, 
                      stdDif=stdDifs, reemplazar=0L, valorReemplazo=NA_real_, stringsAsFactors = F))
  } else {
    return(data.frame(estacion=row.names(coordsObservaciones)[i], fecha=row.names(dfValoresObservaciones)[itsATestear], 
                      valor=dfValoresObservaciones[itsATestear, i], estimado=NA_real_, tipoOutlier=TTO_SinDatosSuficientesEnVecinos, 
                      stdDif=NA_real_, reemplazar=0L, valorReemplazo=NA_real_, stringsAsFactors = F))
  }
}

#' Applies the isolated precipitation/dryness and too large deviation tests
#' We define a neighbourhoud of radius maxDistKm around the target station. For any given date where
#' there are at least minNVecinosPorCuadrante in atl least minNCuadrantes the tests will be run.
#' The isolated precipitation test checks for the value being observed wether, if it's larger than a 
#' given threshold, all neighbours have close to zero precipitation.
#' Conversely the isolated dryness test checks wether for a close to zero value, all neighbours have
#' precipitation above a threshold.
#' The large deviations test uses an estimate obtained from neighbours, which currently is 
#' implemented using an inversely distance weighted interpolation, as well as the amplitude in
#' neighbouring observations, the max - min value in neighbours.
#' If the difference between the observed value and the estimate divided by the amplitude exceeds
#' or is smaller than a certain ratio, the value is considered a high or low valued outlier 
#' respectively.
#' @param coordsObservaciones SpatialPointsDataframe object containing station locations.
#' @param fechasObservaciones dates to which the rows in valoresObservaciones correspond.
#' @param valoresObservaciones matrix containing actual rainfall observations. Columns should be
#' the time series of station observations.
#' @param iColumnasATestear index of the columns in valoresObservaciones to apply the tests to.
#' @param itsATestear index of the rows in valoresObservaciones to apply the tests to.
#' @param ispMax maximum neighbour value to consider a case of isolated precipitation. If at least 
#' one neighbour has a larger value than this, the test will reject isolated precipitation.
#' @param ispObs minimum observed value to consider a case of isolated precipitation. If the value
#' being tested is smaller than this it won't be considered isolated precipitation.
#' @param isdMin minimum neighbour value to consider a case of isolated dryness. If the minimum 
#' neighbour value is smaller than this the test will reject isolated dryness.
#' @param isdObs maximum observed value to consider a case of isolated dryness. If the value being
#' tested is larger than this it won't be considered isolated dryness.
#' @param isdEstMin minimum neighbour derived estimate to consider a case of isolated dryness. If 
#' the neigbour based precipitation estimation is smaller than this it won't be considered isolated
#' dryness.
#' @param fInf ratio of difference between observation and neighbour derived estimate to observed 
#' neighbour amplitude to consider the case a low outlier. If the observation minus the estimate 
#' divided by the amplitude (max - min) observed in the neighbours is smaller than this ratio, the
#' value is considered a low valued outlier.
#' @param fSup ratio of difference between observation and neighbour derived estimate to observed 
#' neighbour amplitude to consider the case a high outlier. If the observation minus the estimate 
#' divided by the amplitude (max - min) observed in the neighbours is larger than this ratio, the
#' value is considered a high valued outlier.
#' @param amplitudMin minimum amplitude observed in neighbours to consider the too large deviations
#' test. If the maximum minus minimum value observed in neighbours is smaller than this the test
#' will not consider the value an outlier.
#' @param minValAbs minimum observable value. Values smaller than this will be considered low valued
#' outliers.
#' @param maxValAbs maximum observable value. Values larger than this will be considered high valued
#' outliers.
#' @param maxDistKm maximum distance in km to consider another station a neighbour.
#' @param minNCuadrantes minimum number of quadrants around the station with observations to apply
#' the tests.
#' @param minNVecinosPorCuadrante minimum number of neighbours in each quadrant with observations to
#' apply the tests.
#' @param datosProtegidos used to feed back the model with corrections. If datosProtegidos[i, j] is
#' true, valoresObservaciones[i, j] will be assumed to have been manually verified and will not be
#' tested for.
#' @param archivoSalida path to store the test results in.
#' @param nCoresAUsar number of CPU cores to run the tests with. Values less than or equal to 0 will
#' use all available CPU cores.
#' @param verbose if verbose and nCoresAUsar == 1, print progress information in console
#' @return a dataframe of test results containing for each observation: station name, date, value,
#' neighbourhood based estimate, test result (tipoOutlier), standardized difference (difference 
#' between observation and estimate divided by amplitude).
testEspacialPrecipitacion <- function(
  coordsObservaciones, fechasObservaciones, valoresObservaciones, 
  iColumnasATestear=seq.int(from = 1, to = ncol(valoresObservaciones), by = 1), 
  itsATestear=1:nrow(valoresObservaciones), ispMax=0.3, ispObs=8, isdMin=1, isdObs=0.3, isdEstMin=3,
  fInf=1, fSup=3, amplitudMin=1, minValAbs=0, maxValAbs=450, maxDistKm=50, minNCuadrantes=4, 
  minNVecinosPorCuadrante=1, 
  datosProtegidos=matrix(data=FALSE, nrow=nrow(valoresObservaciones), ncol=ncol(valoresObservaciones)),
  archivoSalida=NULL, nCoresAUsar=0, verbose=FALSE) {

  maxDist <- distKmToP4Str(proj4string(coordsObservaciones), maxDistKm)
  iesVecinos <- getIVecinosAMenosDeMaxDist(coordsObservaciones = coordsObservaciones, maxDist = maxDist)
  # cuadrantesIes <- lapply(seq_along(iesVecinos), FUN = iCuadrantes, iesVecinos=ies, coords=coords)
  
  if (F) {
    # TO-DO: Para permitir rotar los angulos de los cuadrantes
    diffs <- coords[ies[[i]], ] - coords[i, ]
    anguloCuadrantes = 4 * pi / minNCuadrantes
    
    anguloCuadrantes * 90 / pi
    angulos * 90 / pi
    
    # angulosConVecinos <- atan2(diffs[,1], diffs[,2])
  }
  
  dfValoresObservaciones <- as.data.frame(valoresObservaciones)
  dfValoresObservaciones[dfValoresObservaciones < minValAbs | dfValoresObservaciones > maxValAbs] <- NA
  iNoNAs <- !is.na(dfValoresObservaciones)
  
  if (nCoresAUsar <= 0) nCoresAUsar <- min(detectCores(T, T), ncol(valoresObservaciones))
  if (nCoresAUsar > 1) {
    cl <- makeCluster(getOption('cl.cores', nCoresAUsar))
    clusterExport(cl, varlist = c('script.dir.qcTests'))
    if (exists(x = 'setMKLthreads')) { clusterEvalQ(cl = cl, expr = setMKLthreads(1)) }
    clusterEvalQ(cl = cl, expr = source(paste(script.dir.qcTests, 'qcTests.r', sep='')))
    test <- parLapplyLB(
      cl = cl, X=iColumnasATestear, fun = testEspacialPrecipitacionIV2, iesVecinos=iesVecinos, 
      coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, 
      dfValoresObservaciones=dfValoresObservaciones, iNoNAs=iNoNAs, itsATestear=itsATestear, 
      ispMax=ispMax, ispObs=ispObs, isdMin=isdMin, isdObs=isdObs, isdEstMin=isdEstMin, fInf=fInf, 
      fSup=fSup, amplitudMin=amplitudMin, maxDist=maxDist, minNCuadrantes=minNCuadrantes, 
      minNVecinosPorCuadrante=minNVecinosPorCuadrante, datosProtegidos=datosProtegidos, 
      verbose=verbose)
    stopCluster(cl)
  } else {
    test <- lapply(
      X=iColumnasATestear, FUN = testEspacialPrecipitacionIV2, iesVecinos=iesVecinos, 
      coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, 
      dfValoresObservaciones=dfValoresObservaciones, iNoNAs=iNoNAs, itsATestear=itsATestear, 
      ispMax=ispMax, ispObs=ispObs, isdMin=isdMin, isdObs=isdObs, isdEstMin=isdEstMin, fInf=fInf, 
      fSup=fSup, amplitudMin=amplitudMin, maxDist=maxDist, minNCuadrantes=minNCuadrantes, 
      minNVecinosPorCuadrante=minNVecinosPorCuadrante, datosProtegidos=datosProtegidos, 
      verbose=verbose)
  }
  
  # test[[1]][215,]
  # Código para profiling
  #system.time({
  #Rprof()
  #lele <- spatialRegressionTestI(i = 5,
  #                               ies=ies, coordsObservaciones=coordsObservaciones, fechasObservaciones=fechasObservaciones, 
  #                               dfValoresObservaciones=dfValoresObservaciones, ventana=ventana, minStn=minStn, maxStn=maxStn, 
  #                               maxDist=maxDist, elvDiff_m=elvDiff_m, zcolElv=zcolElv, minAdjR2=minAdjR2, f=f, iNoNAs = iNoNAs)
  #Rprof(NULL)
  #})
  #summaryRprof()
  
  test <- do.call(rbind, test)
  row.names(test) <- paste(make.names(test$estacion), test$fecha, sep='_')
  
  if (!is.na(minValAbs)) test$tipoOutlier[which(valoresObservaciones[itsATestear, iColumnasATestear] < minValAbs)] <- 1
  if (!is.na(maxValAbs)) test$tipoOutlier[which(valoresObservaciones[itsATestear, iColumnasATestear] > maxValAbs)] <- 2
  
  if (!is.null(archivoSalida)) {
    dir.create(dirname(archivoSalida), showWarnings = F, recursive = T)
    write.table(test, archivoSalida, sep = '\t', col.names = T, row.names = T)
  }
  
  return(test)
}
