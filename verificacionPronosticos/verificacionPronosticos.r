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
} else { script.dir.verificacionPronosticos <- paste(dirname(script.dir.verificacionPronosticos), '/', sep='') }

source(paste(script.dir.verificacionPronosticos, '../Graficas/graficas.r', sep=''))

# "Constantes"
StatNames <- c('ME', 'MAE', 'MAD', 'MSE', 'VarDif', 'RMSE', 'Corr', 'RankCorr', 'CorrAnom', 'RankCorrAnom', 'Cant. Datos')
dfInfoValidationStats <- data.frame(StatNames = c('ME', 'MAE', 'MAD', 'MSE', 'VarDif', 'RMSE', 'Corr', 'RankCorr', 'CorrAnom', 'RankCorrAnom', 'Cant. Datos'),
                                    LongStatNames=c('Mean Error', 'Mean Absolute Error', 'Mean Absolute Deviation', 'Mean Squared Error', 'Difference Variance', 
                                                    'Root Mean Squared Error', 'Pearson\'s Coefficient of Correlation', 'Spearman\'s Coefficient of Correlation', 'Anomaly Correlation (Pearson)', 
                                                    'Anomaly Correlation (Spearman)', 'Cantidad de Datos Utilizados'),
                                    valoresPerfectos = c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 'max'),
                                    peoresValores = c('max', 'max', 'max', 'max', 'max', 'max', 0, 0, 0, 0, 0),
                                    minEscala = c('min', 0, 0, 0, 0, 0, -1, -1, -1, -1, 0),
                                    medioEscala = c(0, NA, NA, NA, NA, NA, 0, 0, 0, 0, NA),
                                    maxEscala = c('max', 'max', 'max', 'max', 'max', 'max', 1, 1, 1, 1, 'max'),
                                    invertirColores = c(F, T, T, T, T, T, F, F, F, F, F), stringsAsFactors = F)

#cbind(StatNames, valoresPerfectos, peoresValores, minEscala, medioEscala, maxEscala)

# Funciones
calcValidationStatisticsEx <- function(pronostico, observacion, climatologia) {
  # Saco los nulos para acelerar las operaciones
  i <- !is.na(pronostico) & !is.na(observacion)
  CantDatos <- sum(i) 
  
  if (CantDatos > 0) {
    p <- pronostico[i]
    o <- observacion[i]
    
    #instant_pkgs('Metrics')
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
  } else {
    ME <- NA
    MAE <- NA
    MAD <- NA
    MSE <- NA
    RMSE <- NA
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
  
  return(c(ME, MAE, MAD, MSE, VarDif, RMSE, Corr, RankCorr, CorrAnom, RankCorrAnom, CantDatos))
}

calcValidationStatistics <- function(pronostico, observacion, climatologia) {
  return(as.data.frame(t(setNames(calcValidationStatisticsEx(pronostico, observacion, climatologia), StatNames))))
}

calcValidationStatisticsTemporal <- function(pronosticos, observaciones, climatologias) {
  res <- matrix(NA, nrow=nrow(observaciones), ncol = length(StatNames))
  for (i in 1:nrow(observaciones)) {
    res[i, ] <- calcValidationStatisticsEx(pronostico = pronosticos[i, ], observacion = observaciones[i, ], climatologia = climatologias[i, ])
  }
  colnames(res) <- StatNames
  rownames(res) <- rownames(observaciones)
  return(as.data.frame(res))
}

calcValidationStatisticsEspacial <- function(pronosticos, observaciones, climatologias) {
  res <- matrix(NA, nrow=ncol(observaciones), ncol = length(StatNames))
  for (i in 1:ncol(observaciones)) {
    res[i, ] <- calcValidationStatisticsEx(pronostico = pronosticos[,i], observacion = observaciones[,i], climatologia = climatologias[, i])
  }
  colnames(res) <- StatNames
  rownames(res) <- colnames(observaciones)
  return(as.data.frame(res))
}

calcValidationStatisticsOverall <- function(nombreModelo, pronosticos, observaciones, climatologias) {
  res <- calcValidationStatistics(pronostico=as.vector(pronosticos), observacion=as.vector(observaciones), climatologia=as.vector(climatologias))
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
  archivoStatsOverall <- paste(carpetaSalida, 'ValidationStats_', validationStats$nombreModelo, '.txt', sep='')
  write.table(x = round(validationStats$statsOverall, 3), sep = '\t', file=archivoStatsOverall)  
  
  archivoStatsTemporales <- paste(carpetaSalida, 'ValidationStatsTemporales_', validationStats$nombreModelo, '.png', sep='')
  
  png(archivoStatsTemporales, width=1920, height=1080)
  par(mfrow=c(4,3))
  for (i in 1:ncol(validationStats$statsTemporales))
    plot(x= fechas, validationStats$statsTemporales[, i], type='l', col='blue', lwd=0.5, main=colnames(validationStats$statsTemporales)[i])
  par(mfrow=c(1,1))
  dev.off()
  
  archivoStatsEspaciales <- paste(carpetaSalida, 'ValidationStatsEspaciales_', validationStats$nombreModelo, '.png', sep='')
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
    coordsObservaciones, shpBase=NULL, xyLims=NULL, nColsPlots = 3, ordenModelosPorColumnas=NULL, 
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
  
  validationStats <- list()
  length(validationStats) <- length(pronosticos)
  for (i in 1:length(pronosticos)) {
    validationStats[[i]] <- calcAllValidationStatistics(nombreModelo = names(pronosticos)[i], pronosticos = pronosticos[[i]][iNoNAs, ],
                                                        observaciones = observaciones, climatologias = climatologias)
  }
  names(validationStats) <- names(pronosticos)
  
  statsOverall <- t(sapply(validationStats, function(x) {return (x$statsOverall)}))
  rownames(statsOverall) <- names(pronosticos)
  write.table(x = statsOverall, paste(carpetaSalida, 'validationStatsOverall.tsv', sep=''), sep = '\t', dec = '.', row.names = TRUE, col.names = TRUE)
  
  spValidationStats <- SpatialPointsDataFrame(coords = coordsObservaciones, data=data.frame(value=rep(NA, nrow(coordsObservaciones))))
  gs <- list()
  length(gs) <- length(ordenModelosPorColumnas)
  escalaGraficos <- nColsPlots
  dir.create(carpetaSalida, showWarnings = F)
  
  ####### Espaciales ########
  i <- 1
  for (i in 1:ncol(statsOverall)) {
    print(i)
    iEstadistico <- colnames(statsOverall)[i]
    infoVS <- dfInfoValidationStats[i, ]
    vsEspaciales <- sapply(validationStats, function(x) {return(x$statsEspaciales[,i])})
    
    minEscala <- suppressWarnings(as.numeric(infoVS$minEscala)) 
    if (is.na(minEscala)) { minEscala <- get(infoVS$minEscala)(vsEspaciales, na.rm=T) }
    maxEscala <- suppressWarnings(as.numeric(infoVS$maxEscala))
    if (is.na(maxEscala)) { maxEscala <- get(infoVS$maxEscala)(vsEspaciales, na.rm=T) }
    
    if (!is.na(infoVS$medioEscala)) { escala <- crearEscalaTresPuntos(inicio = minEscala, medio = infoVS$medioEscala, fin = maxEscala, intervaloFinalCerrado = T, nIntervalos = 11, continuo = T, space = 'rgb')
    } else { escala <- crearEscalaDosPuntos(inicio = minEscala, fin = maxEscala, brewerPal = 'RdYlGn', invertirPaleta = infoVS$invertirColores, nIntervalos = 11, continuo = T) }
    escala <- ajustarExtremosEscala(escala = escala, datos = vsEspaciales, nDigitos = 2, redondear = TRUE)

    j <- 1
    for (j in seq_along(ordenModelosPorColumnas)) {
      print(j)
      jPron <- which(names(pronosticos) == ordenModelosPorColumnas[j])
      spValidationStats$value <- vsEspaciales[, jPron]
      nomArchivo <- paste(carpetaSalida, sprintf("%02d", i), '-', iEstadistico, '/', sprintf("%02d", j), '-', ordenModelosPorColumnas[j], '.png', sep='')
      gs[[j]] <- mapearPuntosGGPlot(puntos = spValidationStats, shpBase = shpBase, xyLims = xyLims, dibujarTexto = TRUE, escala = escala, 
                                    tamaniosPuntos = tamaniosPuntos, tamanioFuentePuntos = tamanioFuentePuntos, tamanioFuenteEjes = tamanioFuenteEjes, tamanioFuenteTitulo=tamanioFuenteTitulo,
                                    nDigitos = 2, titulo = paste(ordenModelosPorColumnas[j], ': ', iEstadistico, sep=''), 
                                    nomArchResultados = nomArchivo, dibujar = F, alturaEscalaContinua = unit(x=0.65, units = 'in'), 
                                    escalaGraficos = escalaGraficos)
    }
    
    nomArchMapa <- paste(carpetaSalida, sprintf("%02d", i), '-', iEstadistico, '_Espacial.png', sep='')
    png(nomArchMapa, width = 1920 * escalaGraficos, height = 1017 * escalaGraficos, type='cairo')
    tryCatch(expr = print(multiplot(plotlist=gs, cols = nColsPlots)), finally = dev.off())
  }
  
  ####### Temporales ########
  meses <- months(fechas)
  meses <- factor(meses, levels = unique(meses), ordered = T)
  gs2 <- list()
  length(gs2) <- length(ordenModelosPorColumnas)
  i <- 1
  for (i in 1:ncol(statsOverall)) {
    iEstadistico <- colnames(statsOverall)[i]
    infoVS <- dfInfoValidationStats[i, ]
    vsTemporales <- sapply(validationStats, function(x) {return(x$statsTemporales[, i])})
    
    minEscala <- suppressWarnings(as.numeric(infoVS$minEscala)) 
    if (is.na(minEscala)) { minEscala <- get(infoVS$minEscala)(vsTemporales, na.rm=T) }
    maxEscala <- suppressWarnings(as.numeric(infoVS$maxEscala))
    if (is.na(maxEscala)) { maxEscala <- get(infoVS$maxEscala)(vsTemporales, na.rm=T) }
    
    xyLimsLinePlots <- crearXYLims(min(fechas), max(fechas), minEscala, maxEscala)
    xyLimsBoxPlots <- crearXYLims(min(meses), max(meses), minEscala, maxEscala)

    j <- 4
    for (j in 1:length(ordenModelosPorColumnas)) {
      print(j)
      jPron <- which(names(pronosticos) == ordenModelosPorColumnas[j])
      vsTemporalesJ <- vsTemporales[, jPron]
      
      gs[[j]] <- linePlot(x= fechas, y=vsTemporalesJ, dibujarPuntos = F, titulo = paste(ordenModelosPorColumnas[j], ': ', iEstadistico, sep=''), 
                          xyLims = xyLimsLinePlots, dibujar = F, escalaGraficos=escalaGraficos, lineaRegresion = T, formulaRegresion = y~1,
                          annotateMean = T)
      
      gs2[[j]] <- boxplot_GGPlot(clases = meses, y = vsTemporalesJ, titulo = paste(ordenModelosPorColumnas[j], ': ', iEstadistico, sep=''), 
                                 escalaGraficos = escalaGraficos, xyLims = xyLimsBoxPlots)
    }
    
    which.min(validationStats$`GRK-LST_Night_Combinada_Clim_mean+x+y`$statsTemporales$ME)
    fechas[331]
    
    nomArchGraficos <- paste(carpetaSalida, sprintf("%02d", i), '-', iEstadistico, '_Temporal.png', sep='')
    png(nomArchGraficos, width = 1920 * escalaGraficos, height = 1017 * escalaGraficos, type='cairo')
    tryCatch(expr = print(multiplot(plotlist=gs, cols = nColsPlots)), finally = dev.off())
    
    nomArchGraficos <- paste(carpetaSalida, sprintf("%02d", i), '-', iEstadistico, '_TemporalBoxPlots.png', sep='')
    png(nomArchGraficos, width = 1920 * escalaGraficos, height = 1017 * escalaGraficos, type='cairo')
    tryCatch(expr = print(multiplot(plotlist=gs2, cols = nColsPlots)), finally = dev.off())
  }
  
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
