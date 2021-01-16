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

iFrame <- sys.nframe()
if (iFrame >= 3) { script.dir.leerEscalas <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.leerEscalas <- NULL }
while ((is.null(script.dir.leerEscalas) || is.na(regexpr('leerEscalas.r', script.dir.leerEscalas, fixed=T)[1])) && iFrame >= 0) {
  script.dir.leerEscalas <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.leerEscalas)) { script.dir.leerEscalas <- ''
} else { script.dir.leerEscalas <- paste(dirname(script.dir.leerEscalas), '/', sep='') }

source(paste0(script.dir.leerEscalas, 'mapearEx.r'), encoding = 'WINDOWS-1252')
source(paste0(script.dir.leerEscalas, '../instalarPaquetes/instant_pkgs.r'), encoding = 'WINDOWS-1252')

instant_pkgs(c('jsonlite'))
library(jsonlite)

leerEscala <- function(nombreArchivo, datos, numDecimales=NULL) {
  if (file.exists(nombreArchivo)) {
    fEscala <- file(nombreArchivo, open='r')
    
    linea <- readLines(con=fEscala, n=1, ok=F)
    posIgual <- regexpr('=', linea, fixed=T)[1]
    
    if (posIgual > -1) {
      nomParam <- substr(x=linea, start=1, stop=posIgual - 1)
      valorParam <- substr(x=linea, start=posIgual+1, stop=nchar(linea[i]))
    } else {
      # formato viejo
      seek(fEscala, 0, rw = "r")
      valorParam <- 'TEspecificacionEscalaFija'
    }
    
    escala <- read.table(file=fEscala, sep=" ", dec=".", header=T, colClasses=c("numeric", "character"), comment.char="")
    close(fEscala)
    
    if (valorParam == 'TEspecificacionEscalaFija') {
      minDato <- trunc(min(datos, na.rm=T) * 10^numDecimales) / 10^numDecimales
      maxDato <- ceiling(max(datos, na.rm=T) * 10^numDecimales) / 10^numDecimales
      if (escala$escala[1] > minDato) { escala$escala[1] <- minDato }
      
      if (!is.null(numDecimales)) {
        if (escala$colores[nrow(escala)] == "" & escala$escala[nrow(escala)] < maxDato) { 
          escala$escala[nrow(escala)] <- maxDato
        }
      }
      rm(minDato, maxDato)
    } else if (valorParam == 'TEspecificacionEscalaRelativaAlMinimoYMaximo') {
      minDato <- trunc(min(datos, na.rm=T) * 10^numDecimales) / 10^numDecimales
      maxDato <- ceiling(max(datos, na.rm=T) * 10^numDecimales) / 10^numDecimales
      
      escala$escala <- escala$escala / 100 * (maxDato - minDato) + minDato
    } else if (valorParam == 'TEspecificacionEscalaCuantil') {
      escala$escala <- quantile(datos, probs=escala$escala/100, na.rm=T)
    }
    
    if (!is.null(numDecimales)) escala$escala <- round(escala$escala, numDecimales)
    escala <- crearEscala(escala=escala$escala, colores=escala$colores)
  } else {
    escala <- crearEscala(escala=seq(min(datos), max(datos), length.out=9),
                          colores=brewer.pal(9, 'YlOrRd'))
  }
  return(escala)
}

leerEscalaJSON <- function(nombreArchivo) {
  aux <- fromJSON(nombreArchivo)
  return(crearEscala(aux$escala, aux$colores, aux$continuo))
}

redondearEscala <- function(escala, nDigitos=1) {
  n <- length(escala)
  return(c(achicarToNDigitos(escala[1], nDigitos = nDigitos), round(escala[2:(n-1)], nDigitos), 
           agrandarToNDigitos(escala[n], nDigitos = nDigitos)))
}

leerEscalas <- function(params, datosGrillaPredicciones, datosGrillaSD) {
  # Escalas
  if (params$graficarObservaciones || params$graficarPrediccion) {
    escalaObservaciones <- leerEscalaJSON(nombreArchivo=paste(params$pathProceso, "escalaObservaciones.json", sep=''))
    escalaObservaciones <- ajustarExtremosEscala(escala = escalaObservaciones, datos = datosGrillaPredicciones, nDigitos = params$numDecimales)
  } else { escalaObservaciones <- NULL }
  
  if (params$graficarVarianza || params$graficarSD) {
    escalaSD <- leerEscalaJSON(nombreArchivo=paste(params$pathProceso, "escalaSD.json", sep=''))
    escalaSD <- ajustarExtremosEscala(escala = escalaSD, datos = datosGrillaSD, nDigitos = params$numDecimales)
  } else { escalaSD <- NULL }
  
  return (list(escalaObservaciones=escalaObservaciones, escalaSD=escalaSD))
}

leerEspecificacionEscala <- function(jsonEspecificacionEscala) {
  return(fromJSON(jsonEspecificacionEscala))
}

leerEspecificacionesEscalas <- function(nEspecificaciones, baseNomArchEscalas, extensionArchsEscalas='.json') {
  especificaciones <- vector(mode = "list", length = nEspecificaciones)
  if (nEspecificaciones > 0) {
    for (i in 1:nEspecificaciones) {
      especificaciones[[i]] <- fromJSON(paste(baseNomArchEscalas, i - 1, extensionArchsEscalas, sep=''))
    }
  }
  return(especificaciones)
}

darEscala <- function(especificacion, valores, ajustarExtremos=T) {
  if (especificacion$Clase == 'TEspecificacionEscalaFija') {
    # TEspecificacionEscalaFija
    escala <- crearEscala(escala = redondearEscala(especificacion$iniciosIntervalos, especificacion$nDigitos), colores = especificacion$colores, 
                          brewerPal = especificacion$brewerPal, continuo=especificacion$continuo,
                          iniciosIntervalosIsoLineas = especificacion$iniciosIntervalosIsoLineas)
  } else if (especificacion$Clase == 'TEspecificacionEscalaRelativaAlMinimoYMaximo') {
    # TEspecificacionEscalaRelativaAlMinimoYMaximo
    rango <- range(valores, na.rm = T)
    min <- rango[1]
    max <- rango[2]
    escala <- crearEscala(
      escala = redondearEscala(((max - min) * 0.01) * especificacion$iniciosIntervalos  + min, nDigitos = especificacion$nDigitos),
      colores = especificacion$colores, brewerPal = especificacion$brewerPal, 
      continuo=especificacion$continuo)
  } else if (especificacion$Clase == 'TEspecificacionEscalaRelativaAlMinimoYMaximoDistinguir0') {
    # TEspecificacionEscalaRelativaAlMinimoYMaximoDistinguir0
    minVal <- 10^(-especificacion$nDigitos)
    i <- which(!is.na(valores) & valores >= minVal)
    if (length(i) > 0) {
      rango <- range(valores[i])
      min <- rango[1]
      max <- rango[2]
      escala <- redondearEscala(((max - min) * 0.01) * especificacion$iniciosIntervalos  + min, nDigitos = especificacion$nDigitos)
      escala <- crearEscala(
        escala = escala, colores = especificacion$colores, brewerPal = especificacion$brewerPal, 
        continuo=especificacion$continuo)
      escala$escala <- c(0, escala$escala)
      escala$colores <- c(especificacion$colorCero, escala$colores)
    } else {
      escala <- c(0, 1)
      escala <- crearEscala(
        escala = escala, colores = especificacion$colores, brewerPal = especificacion$brewerPal, 
        continuo=especificacion$continuo)
      escala$colores <- especificacion$colorCero
    }
  } else if (especificacion$Clase == 'TEspecificacionEscalaCuantil') {
    # TEspecificacionEscalaCuantil
    escala <- crearEscala(escala = redondearEscala(escala = quantile(x = valores, probs = especificacion$iniciosIntervalos/100, na.rm = T), nDigitos = especificacion$nDigitos),
                          colores = especificacion$colores, brewerPal = especificacion$brewerPal, continuo=especificacion$continuo) 
  } else if (especificacion$Clase == 'TEspecificacionEscalaRelativaAlMinimoYMaximo_PorNeutro') {
    # TEspecificacionEscalaRelativaAlMinimoYMaximo_PorNeutro
    rango <- range(valores, na.rm = T)
    min <- rango[1]
    max <- rango[2]
    
    cotaInf <- ifelse(especificacion$cotaInf == NULL_DOUBLE, yes = NA, no = especificacion$cotaInf)
    valorNeutro <- especificacion$valorNeutro
    cotaSup <- ifelse(especificacion$cotaSup == NULL_DOUBLE, yes = NA, no = especificacion$cotaSup)
    colorInf <- especificacion$colorInf
    colorNeutro <- especificacion$colorNeutro
    colorSup <- especificacion$colorSup
    espacioColores <- especificacion$espacioColores
    nIntervalosSiTodosAUnLadoDelNeutro <- especificacion$nIntervalosSiTodosAUnLadoDelNeutro
    nIntervalos <- especificacion$nIntervalos
    
    if (max - min < 1E-3) {
      escala <- min
      if (min > valorNeutro) {
        if (!is.na(cotaSup)) {
          colores <- colorRamp(c(colorNeutro, colorSup), space=espacioColores)((min - valorNeutro) / (cotaSup - valorNeutro))
        } else { colores <- colorSup }
      } else if (min < valorNeutro) {
        if (!is.na(cotaInf)) {
          colores <- colorRamp(c(colorInf, colorNeutro), space=espacioColores)((min - cotaInf) / (valorNeutro - cotaInf))
        } else { colores <- TColorDynArray.Create(colorInf) }
      } else { colores <- TColorDynArray.Create(colorNeutro) }
    } else {
      if (min >= valorNeutro) {
        escala <- seq(from=min, to = max, length.out = nIntervalosSiTodosAUnLadoDelNeutro + 1)
        colores <- colorRamp(c(colorNeutro, colorSup), space=espacioColores)(seq(from=0, to = 1, length.out = nIntervalosSiTodosAUnLadoDelNeutro))
      } else if (max <= valorNeutro) {
        escala <- seq(from=min, to = max, length.out = nIntervalosSiTodosAUnLadoDelNeutro + 1)
        colores <- colorRamp(c(colorInf, colorNeutro), space=espacioColores)(seq(from=0, to = 1, length.out = nIntervalosSiTodosAUnLadoDelNeutro))
      } else {
        escala <- seq(from=min, to = max, length.out = nIntervalos + 1)
        # TO-DO: Esto no está bien, los colores no deberían repartirse según 0, 0.5, 1, sino según cotaInf, valorNeutro, cotaSup
        colores <- colorRamp(c(colorInf, colorNeutro, colorSup), space=espacioColores)(seq(from=0, to = 1, length.out = nIntervalos))
      }
    }

    escala <- crearEscala(escala = redondearEscala(escala, nDigitos = especificacion$nDigitos), colores = colores, continuo=especificacion$continuo)
  } else {
    stop(paste('leerEscalas.darEscala: clase de especificación de escala no implementada ', especificacion$Clase, sep=''))
  }
  if (ajustarExtremos) {
    escala <- ajustarExtremosEscala(escala, datos=valores, nDigitos = especificacion$nDigitos, redondear = F)
  }
  return(escala)
}

darEscalas <- function(especificaciones, valores) {
  menosInf <- -.Machine$double.xmax
  masInf <- .Machine$double.xmax
  escalas <- vector(mode = "list", length = length(especificaciones))
  if (length(especificaciones) > 0) {
    i<-1
    for (i in seq.int(from = 1, to = length(especificaciones))) {
      escalas[[i]] <- darEscala(especificacion = especificaciones[[i]], valores = valores)
      iMenosInf <- which(escalas[[i]]$escala == menosInf)
      if (length(iMenosInf) > 0) escalas[[i]]$escala[iMenosInf] <- -Inf
      iMasInf <- which(escalas[[i]]$escala == masInf)
      if (length(iMasInf) > 0) escalas[[i]]$escala[iMasInf] <- Inf
    }
  }
  return(escalas)    
}

crearEspecificacionEscalaFija <- function(iniciosIntervalos, colores=NULL, brewerPal='Spectral', continuo=F, nDigitos=1) {
  return(list(Clase='TEspecificacionEscalaFija', iniciosIntervalos=iniciosIntervalos, colores=colores, brewerPal=brewerPal, continuo=continuo, nDigitos=nDigitos))
}

crearEspecificacionEscalaRelativaAlMinimoYMaximo <- function(iniciosIntervalos=seq(from=0, to=100, length.out = 11), colores=NULL, brewerPal='Spectral', continuo=F, nDigitos=1) {
  return(list(Clase='TEspecificacionEscalaRelativaAlMinimoYMaximo', iniciosIntervalos=iniciosIntervalos, colores=colores, brewerPal=brewerPal, continuo=continuo, nDigitos=nDigitos))
}

crearEspecificacionEscalaRelativaAlMinimoYMaximoDistinguir0 <- function(iniciosIntervalos=seq(from=0, to=100, length.out = 11), colores=NULL, colorCero=rgb(245L, 245L, 249L, alpha = 0, maxColorValue = 255L), brewerPal='Spectral', continuo=F, nDigitos=1) {
  return(list(Clase='TEspecificacionEscalaRelativaAlMinimoYMaximoDistinguir0', iniciosIntervalos=iniciosIntervalos, colores=colores, colorCero=colorCero, brewerPal=brewerPal, continuo=continuo, nDigitos=nDigitos))
}

crearEspecificacionEscalaCuantil <- function(iniciosIntervalos, colores=NULL, brewerPal='Spectral', continuo=F, nDigitos=1) {
  return(list(Clase='TEspecificacionEscalaCuantil', iniciosIntervalos=iniciosIntervalos, colores=colores, brewerPal=brewerPal, continuo=continuo, nDigitos=nDigitos))
}

crearEspecificacionEscalaBoletinPluviometrico <- function() {
  iniciosIntervalos <- c(0, 0.1, 5, 25, 75, 150)
  #c(rgb(red = 240, green = 240, blue = 240, maxColorValue = 255),
  #  colorRampPalette(colors = c(rgb(red = 204, green = 255, blue = 255, maxColorValue = 255),
  #                              rgb(red = 0, green = 0, blue = 102, maxColorValue = 255)), 
  #                   space = 'rgb')(length(iniciosIntervalos)-1))
  colores <- c('#F0F0F0', '#CCFFFF', '#99BFD9', '#6680B2', '#33408C', '#000066')
  return(crearEspecificacionEscalaFija(iniciosIntervalos = iniciosIntervalos, colores = colores, continuo = F, nDigitos = 1))
}

crearEspecificacionEscalaCuantilParaLluvia <- function(datos, nIntervalos=9, continuo=F, nDigitos=0) {
  iniciosIntervalos <- c(0, quantile(datos[datos > 0], probs=seq(from=0, to = 1, length.out = nIntervalos)))
  iniciosIntervalos[length(iniciosIntervalos)] <- agrandarToNDigitos(x = iniciosIntervalos[length(iniciosIntervalos)], nDigitos = nDigitos)
  iniciosIntervalos <- round(iniciosIntervalos, digits = nDigitos)
  colores <- c(rgb(red = 240, green = 240, blue = 240, maxColorValue = 255),
               colorRampPalette(colors = c(rgb(red = 204, green = 255, blue = 255, maxColorValue = 255),
                                           rgb(red = 0, green = 0, blue = 102, maxColorValue = 255)), 
                                space = 'rgb')(length(iniciosIntervalos)-2))
  
  #colores <- c('#F0F0F0', brewer.pal(nIntervalos - 1, 'Blues'))
  return(crearEspecificacionEscalaFija(iniciosIntervalos = iniciosIntervalos, colores = colores, continuo = continuo, nDigitos = nDigitos))
}

crearEspecificacionEscalaEquiespaciadaParaLluvia <- function(datos, nIntervalos=9, continuo=F, nDigitos=0) {
  rango <- range(datos[datos>0], na.rm = T)
  iniciosIntervalos <- c(0, seq(from=rango[1], to = rango[2], length.out = nIntervalos))
  iniciosIntervalos[length(iniciosIntervalos)] <- agrandarToNDigitos(x = iniciosIntervalos[length(iniciosIntervalos)], nDigitos = nDigitos)
  iniciosIntervalos <- round(iniciosIntervalos, digits = nDigitos)
  colores <- c(rgb(red = 240, green = 240, blue = 240, maxColorValue = 255),
               colorRampPalette(colors = c(rgb(red = 204, green = 255, blue = 255, maxColorValue = 255),
                                           rgb(red = 0, green = 0, blue = 102, maxColorValue = 255)), 
                                space = 'rgb')(length(iniciosIntervalos)-2))
  
  #colores <- c('#F0F0F0', brewer.pal(nIntervalos - 1, 'Blues'))
  return(crearEspecificacionEscalaFija(iniciosIntervalos = iniciosIntervalos, colores = colores, continuo = continuo, nDigitos = nDigitos))
}
