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
if (iFrame >= 3) { script.dir.verificacionPronosticos <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.verificacionPronosticos <- NULL }
while ((is.null(script.dir.verificacionPronosticos) || is.na(regexpr('verificacionPronosticos.r', script.dir.verificacionPronosticos, fixed=T)[1])) && iFrame >= 0) {
  script.dir.verificacionPronosticos <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.verificacionPronosticos)) { script.dir.verificacionPronosticos <- ''
} else { script.dir.verificacionPronosticos <- paste0(dirname(script.dir.verificacionPronosticos), '/') }

source(paste0(script.dir.verificacionPronosticos, '../Graficas/graficas.r'))

# "Constantes"
StatNames <- c(
  'ME', 'MAE', 'MAD', 'MSE', 'VarDif', 'RMSE', 'L2/3', 'L2/5', 'Corr', 'RankCorr', 'CorrAnom', 
  'RankCorrAnom', 'Cant. Datos')
dfInfoValidationStats <- data.frame(
  StatNames = c(
    'ME', 'MAE', 'MAD', 'MSE', 'VarDif', 'RMSE', 'L2/3', 'L2/5', 'Corr', 'RankCorr', 'CorrAnom', 
    'RankCorrAnom', 'Cant. Datos', "POD", "FAR", "FBS"),
  LongStatNames=c(
    'Mean Error', 'Mean Absolute Error', 'Mean Absolute Deviation', 'Mean Squared Error', 
    'Difference Variance', 'Root Mean Squared Error', '2/3 Norm', 'L2/5 Norm', 
    'Pearson\'s Coefficient of Correlation', 'Spearman\'s Coefficient of Correlation', 
    'Anomaly Correlation (Pearson)', 'Anomaly Correlation (Spearman)', 
    'Cantidad de Datos Utilizados', 'Probability of Detection', 'False Alarm Rate', 
    'Frequency Bias'),
  valoresPerfectos = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 'max', 1, 0, 1),
  peoresValores = c('max', 'max', 'max', 'max', 'max', 'max', 'max', 'max', 0, 0, 0, 0, 0, 0, 1, 'max'),
  minEscala = c('min', 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0),
  medioEscala = c(0, NA, NA, NA, NA, NA, NA, NA, 0, 0, 0, 0, NA, 0.5, 0.5, NA),
  maxEscala = c(
    'max', 'max', 'max', 'max', 'max', 'max', 'max', 'max', 1, 1, 1, 1, 'max', 1, 1, 'max'),
  invertirColores = c(F, T, T, T, T, T, T, T, F, F, F, F, F, F, F, F), stringsAsFactors = F)

#cbind(StatNames, valoresPerfectos, peoresValores, minEscala, medioEscala, maxEscala)

# Funciones
calcValidationStatisticsEx <- function(pronostico, observacion, climatologia) {
  # Saco los nulos para acelerar las operaciones
  i <- !is.na(pronostico) & !is.na(observacion)
  CantDatos <- sum(i) 
  
  if (CantDatos > 0) {
    p <- pronostico[i]
    o <- observacion[i]
    
    #write.table(cbind(p, o), file = 'pruebaValidationStats.tsv', sep = '\t', dec = '.')
    # Auxiliares
    dif <- p - o
    absDif <- abs(dif)
    
    # Calculo los estadísticos
    ME <- mean(dif)
    MAE <- mean(absDif)
    MAD <- mad(dif)
    MSE <- mean(dif^2)
    RMSE <- sqrt(MSE)
    
    L2_3 <- mean(absDif ^ (2/3)) ^ (3/2)
    L2_5 <- mean(absDif ^ (2/5)) ^ (5/2)
  } else {
    ME <- NA
    MAE <- NA
    MAD <- NA
    MSE <- NA
    RMSE <- NA
    L2_3 <- NA
    L2_5 <- NA
  }
  
  if (CantDatos > 1) { VarDif <- var(dif)
  } else { VarDif <- NA }

  if (CantDatos > 2) {
    if (var(p) > 0 && var(o) > 0) {
      Corr <- cor(p, o, method = 'pearson')
      RankCorr <- cor(p, o, method = 'spearman')
    } else {
      Corr <- NA
      RankCorr <- NA
    }
    
    if (!is.null(climatologia)) {
      # Anomalías
      clim <- climatologia[i]
      a_p <- p - clim
      a_o <- o - clim
      
      if (var(a_p) > 0 && var(a_o) > 0) {
        CorrAnom <- cor(a_p, a_o, method = 'pearson')
        RankCorrAnom <- cor(a_p, a_o, method = 'spearman')
      } else {
        CorrAnom <- NA
        RankCorrAnom <- NA
      }
    } else {
      CorrAnom <- NA
      RankCorrAnom <- NA
    }
  } else {
    Corr <- NA
    RankCorr <- NA
    CorrAnom <- NA
    RankCorrAnom <- NA
  }  
  
  #MAE_Clim <- mean(abs(a_o))
  #SS_MAE <- 1 - (MAE / MAE_Clim)
  #MSE_Clim <- mean(a_o^2)
  #SS_MSE <- 1 - (MSE / MSE_Clim)
  
  return(c(ME, MAE, MAD, MSE, VarDif, RMSE, L2_3, L2_5, Corr, RankCorr, CorrAnom, RankCorrAnom, CantDatos))
}

calcRainfallDetectionStatistics <- function(
    pronostico, observacion, threshold=NULL, positive='TRUE') {
  # threshold <- thresholds[1]
  idx <- !is.na(pronostico) & !is.na(observacion)

  if (!is.null(threshold)) {
    p <- pronostico[idx] > threshold
    o <- observacion[idx] > threshold
  } else {
    p <- pronostico[idx]
    o <- threshold[idx]
  }
  p <- as.factor(p)
  o <- as.factor(o)
  
  if (any(o == positive)) {
    FBS <- sum(p == positive) / sum(o == positive)
    POD <- sum(p == positive & o == positive) / sum(o == positive)
  } else {
    FBS <- NA_real_
    POD <- NA_real_
  }
  
  if (any(p == positive)) {
    FAR <- sum(p == positive & o != positive) / sum(p == positive)
  } else {
    FAR <- NA_real_
  }

  return(structure(c(POD, FAR, FBS), names=c('POD', 'FAR', 'FBS')))
}

calcRainfallDetectionMultiThresholds <- function(
    pronostico, observacion, thresholds, positive='TRUE') {
  rainfallDetectionStats <- sapply(thresholds, function(x) {
    calcRainfallDetectionStatistics(
      pronostico=pronostico, observacion=observacion, threshold=x, positive=positive)
  })
  
  resNames <- apply(
    expand.grid(rownames(rainfallDetectionStats), as.character(thresholds)), 
    MARGIN = 1, FUN = paste, sep='', collapse='_')
  return(structure(as.vector(rainfallDetectionStats), names=resNames))
}

calcRainfallDetectionStatisticsEspacial <- function(
    pronosticos, observaciones, thresholds, positive='TRUE') {
  res <- matrix(NA, nrow=ncol(observaciones), ncol = 3 * length(thresholds))
  for (i in 1:ncol(observaciones)) {
    aux <- calcRainfallDetectionMultiThresholds(
      pronostico=pronosticos[,i], observacion=observaciones[,i], thresholds=thresholds, 
      positive=positive)
    res[i, ] <- aux
  }
  colnames(res) <- names(aux)
  rownames(res) <- colnames(observaciones)
  return(as.data.frame(res))
}

calcValidationStatistics <- function(pronostico, observacion, climatologia) {
  return(as.data.frame(t(setNames(calcValidationStatisticsEx(pronostico, observacion, climatologia), 
                                  StatNames))))
}

calcValidationStatisticsTemporal <- function(pronosticos, observaciones, climatologias) {
  res <- matrix(NA, nrow=nrow(observaciones), ncol = length(StatNames))
  for (i in 1:nrow(observaciones)) {
    res[i, ] <- calcValidationStatisticsEx(
      pronostico=pronosticos[i, ], observacion=observaciones[i, ], climatologia=climatologias[i, ])
  }
  colnames(res) <- StatNames
  rownames(res) <- rownames(observaciones)
  return(as.data.frame(res))
}

calcValidationStatisticsEspacial <- function(pronosticos, observaciones, climatologias) {
  res <- matrix(NA, nrow=ncol(observaciones), ncol = length(StatNames))
  for (i in 1:ncol(observaciones)) {
    res[i, ] <- calcValidationStatisticsEx(
      pronostico=pronosticos[,i], observacion=observaciones[,i], climatologia=climatologias[, i])
  }
  colnames(res) <- StatNames
  rownames(res) <- colnames(observaciones)
  return(as.data.frame(res))
}

calcValidationStatisticsOverall <- function(nombreModelo, pronosticos, observaciones, climatologias) {
  res <- calcValidationStatistics(
    pronostico=as.vector(pronosticos), observacion=as.vector(observaciones), 
    climatologia=as.vector(climatologias))
  rownames(res) <- nombreModelo
  return(res)
}

calcAllValidationStatistics <- function(nombreModelo, pronosticos, observaciones, climatologias) {
  statsOverall <- calcValidationStatisticsOverall(nombreModelo, pronosticos, observaciones, climatologias)
  statsTemporales <- calcValidationStatisticsTemporal(pronosticos, observaciones, climatologias)
  statsEspaciales <- calcValidationStatisticsEspacial(pronosticos, observaciones, climatologias)
  
  return(list(nombreModelo=nombreModelo, statsOverall=statsOverall, statsTemporales=statsTemporales, statsEspaciales=statsEspaciales))
}

saveAndPlotValidationStatistics <- function(validationStats, fechas, carpetaSalida='Resultados/Validacion/', 
                                            coordsObservaciones, shpBase=NULL) {
  crearDirectoriosSiNoExisten(carpetaSalida)
  archivoStatsOverall <- paste0(carpetaSalida, 'ValidationStats_', validationStats$nombreModelo, '.txt')
  write.table(x = round(validationStats$statsOverall, 3), sep = '\t', file=archivoStatsOverall)  
  
  archivoStatsTemporales <- paste0(carpetaSalida, 'ValidationStatsTemporales_', validationStats$nombreModelo, '.png')
  
  png(archivoStatsTemporales, width=1920, height=1080)
  par(mfrow=c(4,3))
  for (i in 1:ncol(validationStats$statsTemporales))
    plot(x= fechas, validationStats$statsTemporales[, i], type='l', col='blue', lwd=0.5, main=colnames(validationStats$statsTemporales)[i])
  par(mfrow=c(1,1))
  dev.off()
  
  archivoStatsEspaciales <- paste0(carpetaSalida, 'ValidationStatsEspaciales_', validationStats$nombreModelo, '.png')
  png(archivoStatsEspaciales, width=3940, height=2160)
  
  png(width=3940, height=2160)
  plots <- list()
  length(plots) <- ncol(validationStats$statsEspaciales)
  #i <- 1
  for (i in 1:ncol(validationStats$statsEspaciales)) {
    coordsObservaciones$values <- validationStats$statsEspaciales[, i]
    
    plots[[i]] <- mapearPuntosConEtiquetasGGPlot(puntos = coordsObservaciones, shpBase = shpBase, 
                                                 zcol='values', titulo = colnames(validationStats$statsEspaciales)[i])
  }
  do.call(grid.arrange, c(plots, nrow=3, ncol=4))
  dev.off()
}

calcAndPlotAllValidationStatistics <- function(nombreModelo, fechas, pronosticos, observaciones, climatologias, 
                                               carpetaSalida='Resultados/Validacion/', coordsObservaciones, shpBase=NULL) {
  
  validationStats <- calcAllValidationStatistics(nombreModelo, pronosticos, observaciones, climatologias)
  saveAndPlotValidationStatistics(validationStats = validationStats, fechas = fechas, carpetaSalida = carpetaSalida, 
                                  coordsObservaciones = coordsObservaciones, shpBase = shpBase)
  return(validationStats)
}

calcAndPlotAllValidationStatisticsV2 <- function(
    fechas, pronosticos, observaciones, climatologias, carpetaSalida='Resultados/Validacion2/', 
    coordsObservaciones, shpBase=NULL, xyLims=NULL, 
    nColsPlots=min(3, length(ordenModelosPorColumnas)), ordenModelosPorColumnas=names(pronosticos), 
    tamaniosPuntos = 4, tamanioFuentePuntos = 3, tamanioFuenteEjes = 15, tamanioFuenteTitulo=22) {
  dir.create(carpetaSalida, showWarnings = FALSE, recursive = TRUE)
  if (is.null(ordenModelosPorColumnas)) ordenModelosPorColumnas <- names(pronosticos)
  
  primerYUltimoNoNA <- function(x) {
    iNoNAs <- which(apply(x, 1, FUN = function(x) { any(!is.na(x)) } ))
    return(c(iNoNAs[1], iNoNAs[length(iNoNAs)]))
  }
  primerosYUltimosNoNA <- sapply(pronosticos, FUN = primerYUltimoNoNA)
  primerNoNA <- min(primerosYUltimosNoNA[1,])
  ultimoNoNA <- max(primerosYUltimosNoNA[2,])
  
  iNoNAs <- primerNoNA:ultimoNoNA
  fechas <- fechas[iNoNAs]
  observaciones <- observaciones[iNoNAs,]
  climatologias <- climatologias[iNoNAs,]
  
  validationStats <- calcValidationStatisticsMultipleModels(
    valoresObservaciones=observaciones, cvs=pronosticos, 
    climatologias=climatologias, pathResultados=carpetaSalida)
  
  write.table(
    x=validationStats$validationStatsOverall, 
    file=paste0(carpetaSalida, 'validationStatsOverall.tsv'), sep='\t', dec='.', row.names=TRUE, 
    col.names=TRUE)
  
  ####### Espaciales ########
  plotValidationStatsEspaciales(
    coordsObservaciones=coordsObservaciones, 
    statsEspaciales=validationStats$validationStatsEspaciales[ordenModelosPorColumnas], 
    carpetaSalida=carpetaSalida, shpBase=shpBase, xyLims=xyLims, nColsPlots=nColsPlots, 
    ordenModelosPorColumnas=ordenModelosPorColumnas, tamaniosPuntos=tamaniosPuntos,
    tamanioFuentePuntos=tamanioFuentePuntos, tamanioFuenteEjes=tamanioFuenteEjes, 
    tamanioFuenteTitulo=tamanioFuenteTitulo
  )

  ####### Temporales ########
  plotValidationStatsTemporales(
    fechas=fechas, statsTemporales=validationStats$validationStatsTemporales[ordenModelosPorColumnas], 
    carpetaSalida=carpetaSalida, nColsPlots=nColsPlots, 
    ordenModelosPorColumnas=ordenModelosPorColumnas
  )
  
  return(validationStats)
}

plotCVStationByStation <- function(iEstacion, fechas, valoresObservaciones, pronosticos, clasesFechas, 
                                   formatoFechas="%Y-%m-%d %H:%M") {
  # iEstacion <- which(estaciones$Nombre == 'Artigas')
  valoresEstacionI <- matrix(nrow = length(fechas), ncol=length(pronosticos) + 1)
  valoresEstacionI[,1] <- valoresObservaciones[iNoNAs, iEstacion]
  for (i in seq_along(pronosticos)) {
    valoresEstacionI[,i + 1] <- pronosticos[[i]][iNoNAs, iEstacion]
  }
  colnames(valoresEstacionI) <- c(colnames(valoresObservaciones)[iEstacion], names(pronosticos))
  
  # iClase <- unique(clasesFechas)[1]
  for (iClase in unique(clasesFechas)) {
    iFechas <- which(clasesFechas == iClase)
    fechasIClase <- fechas[iFechas]
    valoresEstacionIiClase <- valoresEstacionI[iFechas,]
    
    linePlot(x=fechasIClase, y=valoresEstacionIiClase, 
             dibujarPuntos=F, titulo=paste(colnames(valoresEstacionI)[1], '-', iClase), tituloEjeX='Tiempo', tituloEjeY='C', 
             dibujarEscala=T, dibujarEjes=T, xyLims=NULL, dibujar=interactive(), DPI=600, tamanioFuenteTextos=15, 
             escalaGraficos=1, annotateMean=F, nomArchSalida='D:/Tesis/Artigas2014.png') 
  }
}

plotValidationStatsEspaciales <- function(
    coordsObservaciones, statsEspaciales, carpetaSalida='Resultados/Validacion2/', shpBase=NULL, 
    xyLims=NULL, nColsPlots=min(3, length(ordenModelosPorColumnas)), 
    ordenModelosPorColumnas=names(statsEspaciales), tamaniosPuntos=4,
    tamanioFuentePuntos=3, tamanioFuenteEjes=15, tamanioFuenteTitulo=22) {
  statsNames <- colnames(statsEspaciales[[1]])
  
  gs <- list()
  length(gs) <- length(ordenModelosPorColumnas)
  escalaGraficos <- nColsPlots
  dir.create(carpetaSalida, showWarnings = F)
  
  i <- 1
  for (i in seq_along(statsNames)) {
    statName <- statsNames[i]
    iEstadistico <- which(startsWith(x=statName, prefix=dfInfoValidationStats$StatNames))
    if (length(iEstadistico) > 1) {
      # Keep the longest match
      iEstadistico <- iEstadistico[
        which.max(sapply(dfInfoValidationStats$StatNames[iEstadistico], nchar))]
    }
    
    infoVS <- dfInfoValidationStats[iEstadistico, ]
    vsEspaciales <- sapply(statsEspaciales[ordenModelosPorColumnas], function(x) {return(x[, i])})
    
    minEscala <- suppressWarnings(as.numeric(infoVS$minEscala)) 
    if (is.na(minEscala)) { minEscala <- get(infoVS$minEscala)(vsEspaciales, na.rm=T) }
    maxEscala <- suppressWarnings(as.numeric(infoVS$maxEscala))
    if (is.na(maxEscala)) { maxEscala <- get(infoVS$maxEscala)(vsEspaciales, na.rm=T) }
    
    if (!is.na(infoVS$medioEscala)) {
      escala <- crearEscalaTresPuntos(
        inicio=minEscala, medio=infoVS$medioEscala, fin=maxEscala, intervaloFinalCerrado=T, 
        nIntervalos=11, continuo=T, space='rgb')
    } else { 
      escala <- crearEscalaDosPuntos(
        inicio = minEscala, fin = maxEscala, brewerPal = 'RdYlGn', 
        invertirPaleta = infoVS$invertirColores, nIntervalos = 11, continuo = T) 
    }
    escala <- ajustarExtremosEscala(escala=escala, datos=vsEspaciales, nDigitos=2, redondear=TRUE)
    
    j <- 1
    for (j in seq_along(ordenModelosPorColumnas)) {
      print(paste0(statName, ': ', ordenModelosPorColumnas[j]))
      coordsObservaciones$value <- vsEspaciales[, j]
      nomArchivo <- paste0(
          carpetaSalida, sprintf("%02d", iEstadistico), '-', statName, '/', sprintf("%02d", j), '-', 
          ordenModelosPorColumnas[j], '.png')
      gs[[j]] <- mapearPuntosGGPlot(
        puntos=coordsObservaciones, shpBase=shpBase, xyLims=xyLims, 
        zcol='value', dibujarTexto=TRUE, 
        escala=escala, tamaniosPuntos=tamaniosPuntos, tamanioFuentePuntos=tamanioFuentePuntos,
        tamanioFuenteEjes=tamanioFuenteEjes, tamanioFuenteTitulo=tamanioFuenteTitulo, nDigitos=2, 
        titulo=paste0(ordenModelosPorColumnas[j], ': ', statsNames[i]), 
        nomArchResultados=nomArchivo, dibujar=F, alturaEscalaContinua=unit(x=0.65, units='in'), 
        escalaGraficos=escalaGraficos)
    }
    
    nomArchMapa <- paste0(
      carpetaSalida, sprintf("%02d", iEstadistico), '-', statName, '_Espacial.png')
    png(nomArchMapa, width = 1920 * escalaGraficos, height = 1017 * escalaGraficos)
    tryCatch(expr = print(multiplot(plotlist=gs, cols = nColsPlots)), finally = dev.off())
  }
}

plotValidationStatsTemporales <- function(
    fechas, statsTemporales, carpetaSalida='Resultados/Validacion2/', 
    nColsPlots=min(3, length(ordenModelosPorColumnas)), 
    ordenModelosPorColumnas=names(statsTemporales)) {
  statsNames <- colnames(statsTemporales[[1]])
  
  if (is.character(fechas)) {
    fechas <- lubridate::ymd_hms(fechas, truncated = 5)
  }
  
  meses <- months(fechas, abbreviate = TRUE)
  meses <- factor(meses, levels = unique(meses), ordered = T)
  gs <- list()
  length(gs) <- length(ordenModelosPorColumnas)
  gs2 <- list()
  length(gs2) <- length(ordenModelosPorColumnas)
  escalaGraficos <- nColsPlots
  dir.create(carpetaSalida, showWarnings = F)
  
  i <- 1
  for (i in seq_along(statsNames)) {
    statName <- statsNames[i]
    iEstadistico <- which(startsWith(x=statName, prefix=dfInfoValidationStats$StatNames))
    if (length(iEstadistico) > 1) {
      # Keep the longest match
      iEstadistico <- iEstadistico[
        which.max(sapply(dfInfoValidationStats$StatNames[iEstadistico], nchar))]
    }
    
    infoVS <- dfInfoValidationStats[iEstadistico, ]
    vsTemporales <- sapply(statsTemporales[ordenModelosPorColumnas], function(x) {return(x[, i])})
    
    minEscala <- suppressWarnings(as.numeric(infoVS$minEscala)) 
    if (is.na(minEscala)) { minEscala <- get(infoVS$minEscala)(vsTemporales, na.rm=T) }
    maxEscala <- suppressWarnings(as.numeric(infoVS$maxEscala))
    if (is.na(maxEscala)) { maxEscala <- get(infoVS$maxEscala)(vsTemporales, na.rm=T) }
    
    xyLimsLinePlots <- crearXYLims(min(fechas), max(fechas), minEscala, maxEscala)
    xyLimsBoxPlots <- crearXYLims(min(meses), max(meses), minEscala, maxEscala)
    
    j <- 1
    for (j in seq_along(ordenModelosPorColumnas)) {
      print(paste0(statName, ': ', ordenModelosPorColumnas[j]))

      vsTemporalesJ <- vsTemporales[, j]

      titulo <- paste0(ordenModelosPorColumnas[j], ': ', statName)      
      gs[[j]] <- linePlot(
        x= fechas, y=vsTemporalesJ, dibujarPuntos=F, titulo=titulo, xyLims=xyLimsLinePlots, 
        dibujar=F, escalaGraficos=escalaGraficos, lineaRegresion=T, formulaRegresion=y~1, 
        annotateMean = T)
      
      gs2[[j]] <- boxplot_GGPlot(
        clases=meses, y=vsTemporalesJ, titulo=titulo, escalaGraficos=escalaGraficos, 
        xyLims=xyLimsBoxPlots)
    }
    
    nomArchGraficos <- paste0(
      carpetaSalida, sprintf("%02d", iEstadistico), '-', statName, '_Temporal.png')
    png(nomArchGraficos, width = 1920 * escalaGraficos, height = 1017 * escalaGraficos)
    tryCatch(expr = print(multiplot(plotlist=gs, cols = nColsPlots)), finally = dev.off())
    
    nomArchGraficos <- paste0(
      carpetaSalida, sprintf("%02d", iEstadistico), '-', statName, '_TemporalBoxPlots.png')
    png(nomArchGraficos, width = 1920 * escalaGraficos, height = 1017 * escalaGraficos)
    tryCatch(expr = print(multiplot(plotlist=gs2, cols = nColsPlots)), finally = dev.off())
  }
}

plotCVStationByStation <- function(fechasObservaciones, valoresObservaciones, pronosticos, clasesFechas=year(fechasObservaciones),
                                   formatoFechas="%Y-%m-%d %H:%M") {
  # pronosticos <- cvs[1:2]
  
  primerYUltimoNoNA <- function(x) {
    iNoNAs <- which(apply(x, 1, FUN = function(x) { any(!is.na(x)) } ))
    return(c(iNoNAs[1], iNoNAs[length(iNoNAs)]))
  }
  primerosYUltimosNoNA <- sapply(pronosticos, FUN = primerYUltimoNoNA)
  primerNoNA <- min(primerosYUltimosNoNA[1,])
  ultimoNoNA <- max(primerosYUltimosNoNA[2,])
  
  iNoNAs <- primerNoNA:ultimoNoNA
  fechas <- fechasObservaciones[iNoNAs]
  clasesFechas <- clasesFechas[iNoNAs]
  
  iEstacion <- 1



  
} 
