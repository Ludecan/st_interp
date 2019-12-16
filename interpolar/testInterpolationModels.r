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
if (iFrame >= 3) { script.dir.testInterpolationModels <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.testInterpolationModels <- NULL }
while ((is.null(script.dir.testInterpolationModels) || is.na(regexpr('testInterpolationModels.r', script.dir.testInterpolationModels, fixed=T)[1])) && iFrame >= 0) {
  script.dir.testInterpolationModels <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.testInterpolationModels)) { script.dir.testInterpolationModels <- ''
} else { script.dir.testInterpolationModels <- paste(dirname(script.dir.testInterpolationModels), '/', sep='') }

source(paste(script.dir.testInterpolationModels, '../verificacionPronosticos/verificacionPronosticos.r', sep=''))

getTSeqs <- function(fechasObservaciones) {
  anios <- sort(unique(year(fechasObservaciones)), decreasing = T)
  getInicioFinAnio <- function(anio, fechasObservaciones) {
    ts <- which(year(fechasObservaciones) == anio)
    tIni <- ts[1]
    tFin <- ts[length(ts)]
    return(c(tIni, tFin))
  }
  tSeqs <- sapply(anios, FUN = getInicioFinAnio, fechasObservaciones = fechasObservaciones)
  return(data.frame(anio=anios, tIni = tSeqs[1,], tFin=tSeqs[2,]))
}

testRegressors <- function(
    valoresObservaciones, pathsRegresores, pathSHPNotNUll, pathResultados='Resultados/', 
    seriesName='Rainfall', outputTableFilename=NULL) {
  carpeta <- paste(pathResultados, '1-Exploracion/', sep='')
  serie <- unlist(c(valoresObservaciones))
  
  res <- matrix(nrow=ncol(pathsRegresores), ncol = 7)
  colnames(res) <- c('Pearson', 'Spearman', 'Adj. R^2', 'CantDatos', 'CoberturaMínima', 'CoberturaMedia', 'CoberturaMáxima')
  rownames(res) <- colnames(pathsRegresores)
  
  iNoNaSerie <- !is.na(serie)
  i <- 1
  for (i in 1:ncol(pathsRegresores)) {
    print(paste(Sys.time(), ' : TestsRegresores ', i, ' - ', colnames(pathsRegresores)[i], sep=''))
    
    aux <- readGDAL(fname = pathsRegresores[which(!is.na(pathsRegresores[, i]))[1], i])
    shpMaskNoNulos <- cargarSHPYObtenerMascaraParaGrilla(pathSHP = pathSHPNotNUll, grilla=aux)
    nPixeles <- sum(shpMaskNoNulos$mask)
    
    # mapearGrillaGGPlot(grilla = aux, shpBase = shpMaskNoNulos$shp, continuo=T)
    rm(aux)
    
    regresor <- as.vector(extraerValoresRegresorSobreSP(
      objSP = coordsObservaciones, pathsRegresor = pathsRegresores[, i, drop=F], silent = F))
    
    # which(is.na(pathsRegresores[, i, drop=F]))[1]
    nrow(pathsRegresores)
    length(serie)
    length(regresor)
    
    iNoNa <- which(iNoNaSerie & !is.na(regresor))
    if (length(iNoNa) > 1) {
      s <- serie[iNoNa]
      r <- regresor[iNoNa] 
      
      res[i, 1] <- cor(s, r)
      res[i, 2] <- cor(s, r, method = 'spearman')
      m <- lm(formula = 'y~x+1', data = data.frame(x=r, y=s))
      smry <- summary(m)
      res[i, 3] <- smry$adj.r.squared
      res[i, 4] <- length(iNoNa)
      
      nNoNulosPorCuadrantes <- contarNoNulosPorCuadrantes(
        pathsGeoTiffs = pathsRegresores[, i, drop=F], shpMask=shpMaskNoNulos, nCoresAUsar = 0)
      nNoNulos <- rowSums(nNoNulosPorCuadrantes[, -ncol(nNoNulosPorCuadrantes)])
      
      # plot(nNoNulos)
      
      aux <- range(nNoNulos)
      res[i, 5] <- aux[1] / nPixeles * 100
      res[i, 6] <- mean(nNoNulos) / nPixeles * 100
      res[i, 7] <- aux[2] / nPixeles * 100
      
      # quantile(nNoNulos, probs=c(0, 0.01, 0.02, 0.5, 0.98, 0.99, 1)) / nPixeles * 100
      
      arch <- paste(carpeta, seriesName, '_vs_', colnames(pathsRegresores)[i], '.png', sep='')
      linePlot(
        x=r, y=s, tituloEjeX=colnames(pathsRegresores)[i], tituloEjeY=seriesName, 
        lineaRegresion = T, intervalosConfianza = T,  dibujarPuntos = T, 
        titulo = paste(seriesName, ' vs ', colnames(pathsRegresores)[i], sep=''), 
        dibujar = interactive(), nomArchSalida = arch)
      
      rm(s, r, m, smry, aux, nNoNulosPorCuadrantes, nNoNulos, arch)
    }
    rm(shpMaskNoNulos, nPixeles, iNoNa, regresor)
  }
  
  if (!is.null(outputTableFilename)) {
    write.table(res, file = paste(carpeta, outputTableFilename), col.names=T, row.names=T, 
                append = F, quote = F, sep = "\t", eol = "\r", dec = ".", 
                qmethod = c("escape", "double"))  
  }
  rm(serie, carpeta)
  
  return(res)
}

st_interpCrossValidation <- function(
    coordsObservaciones, fechasObservaciones, valoresObservaciones, params, pathsRegresores, 
    pathResultados='Resultados/3-GrilladoYCV/', recalcCV=FALSE) {
  nomModelo <- nombreModelo(params = params, pathsRegresores=pathsRegresores)
  archCV <- paste(pathResultados, nomModelo, '/cv/cv_', nomModelo, '.tsv', sep = '')
  if (recalcCV || !file.exists(archCV) || file.info(archCV)$size == 0) {
    dir.create(dirname(archCV), showWarnings = F, recursive = T)
    linea <- paste(Sys.time(), ': CV ', nomModelo, sep='')
    print(linea)
    write(linea, file="tiemposEjecucion.txt", append=TRUE)
    
    #eliminarSerieTemporalCompleta = FALSE
    #longitudesEnColumnas = T
    
    cv <- universalGriddingCV(
      coordsObservaciones = coordsObservaciones, fechasObservaciones = fechasObservaciones, 
      valoresObservaciones = valoresObservaciones, params = params, pathsRegresores = pathsRegresores,
      eliminarSerieTemporalCompleta = FALSE, longitudesEnColumnas = TRUE)
    
    write.table(x = cv, file = archCV, sep = '\t', na = '-99', dec = '.', row.names = T, col.names = T)
  } else {
    cv <- read.csv(archCV, sep = '\t', na.strings = '-99', dec = '.')
  }
  return(cv)
}

st_interpCrossValidations <- function(
    coordsObservaciones, fechasObservaciones, valoresObservaciones, listaParams, listaRegresores, 
    pathResultados='Resultados/3-GrilladoYCV/', recalcCV=FALSE, 
    modelosACorrer=1:length(listaParams)) {
  cvs <- vector(mode="list", length = length(modelosACorrer))
  # i <- 1
  for (i in seq_along(modelosACorrer)) {
    iModel <- modelosACorrer[i]
    params <- listaParams[[iModel]]
    # params$nCoresAUsar <- 1
    if (is.na(listaRegresores[iModel])) { pathsRegresores <- NULL
    } else { pathsRegresores <- listaRegresores[[iModel]] }
    
    names(cvs)[iModel] <- nombreModelo(params = params, pathsRegresores=pathsRegresores)
    
    try(cvs[[i]] <- st_interpCrossValidation(
      coordsObservaciones, fechasObservaciones, valoresObservaciones, params = params, 
      pathsRegresores = pathsRegresores, pathResultados=pathResultados, recalcCV=recalcCV))
  }
  return(cvs)
}

calcValidationStatistics <- function(valoresObservaciones, cvs, climatologias=NULL, 
                                     pathResultados='Resultados/4-Validacion/') {
  source(paste(pathScripts, 'verificacionPronosticos/verificacionPronosticos.r', sep=''))
  
  i <- 1
  # Estadísticos de Validación
  validationStatsOverall <- data.frame()
  validationStatsEspaciales <- list()
  validationStatsTemporales <- list()
  length(validationStatsEspaciales) <- length(cvs)
  length(validationStatsTemporales) <- length(cvs)
  
  for (i in 1:length(cvs)) {
    nomModelo <- names(cvs)[i]
    validationStatsOverall <- rbind(
      validationStatsOverall, 
      calcValidationStatisticsOverall(
        nombreModelo = nomModelo, pronosticos = cvs[[i]], observaciones = valoresObservaciones, 
        climatologias = climatologias))
    
    validationStatsEspaciales[[i]] <- calcValidationStatisticsEspacial(
      pronosticos = cvs[[i]], observaciones = valoresObservaciones, climatologias = climatologias)
    validationStatsTemporales[[i]] <- calcValidationStatisticsTemporal(
      pronosticos = cvs[[i]], observaciones = valoresObservaciones, climatologias = climatologias)
  }
  round(validationStatsOverall, 2)
  dir.create(pathResultados, showWarnings = F, recursive = T)
  write.table(x = validationStatsOverall, paste(pathResultados, 'validationStatsOverall.tsv'), 
              sep = '\t', dec = '.', row.names = TRUE, col.names = TRUE)
  
  return(list(validationStatsOverall=validationStatsOverall,
              validationStatsEspaciales=validationStatsEspaciales,
              validationStatsTemporales=validationStatsTemporales))
}
