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