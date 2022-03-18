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
if (iFrame >= 3) { script.dir.graficas <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.graficas <- NULL }
while ((is.null(script.dir.graficas) || is.na(regexpr('graficas.r', script.dir.graficas, fixed=T)[1])) && iFrame >= 0) {
  script.dir.graficas <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.graficas)) { script.dir.graficas <- ''
} else { script.dir.graficas <- paste0(dirname(script.dir.graficas), '/') }

set1 <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')
paired <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
set3 <- c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
unikn_pair <- c("#3E5496","#8290BB","#008ECE","#59C7EB","#077187","#6AAAB7","#0A9086","#54BFB7","#8E2043","#BC7A8F","#E0607E","#ECA0B2","#FEA090","#FECFC7","#B8BCC1","#E1E2E5")
colores64 <- c('#000000','#00FF00','#0000FF','#FF0000','#01FFFE','#FFA6FE','#FFDB66','#006401','#010067','#95003A',
               '#007DB5','#FF00F6','#FFEEE8','#774D00','#90FB92','#0076FF','#D5FF00','#FF937E','#6A826C','#FF029D',
               '#FE8900','#7A4782','#7E2DD2','#85A900','#FF0056','#A42400','#00AE7E','#683D3B','#BDC6FF','#263400',
               '#BDD393','#00B917','#9E008E','#001544','#C28C9F','#FF74A3','#01D0FF','#004754','#E56FFE','#788231',
               '#0E4CA1','#91D0CB','#BE9970','#968AE8','#BB8800','#43002C','#DEFF74','#00FFC6','#FFE502','#620E00',
               '#008F9C','#98FF52','#7544B1','#B500FF','#00FF78','#FF6E41','#005F39','#6B6882','#5FAD4E','#A75740',
               '#A5FFD2','#FFB167','#009BFF','#E85EBE')

colores269 <- c('#000000','#FFFF00','#1CE6FF','#FF34FF','#FF4A46',
                '#008941','#006FA6','#A30059','#FFDBE5','#7A4900','#0000A6','#63FFAC','#B79762','#004D43','#8FB0FF',
                '#997D87','#5A0007','#809693','#FEFFE6','#1B4400','#4FC601','#3B5DFF','#4A3B53','#FF2F80','#61615A',
                '#BA0900','#6B7900','#00C2A0','#FFAA92','#FF90C9','#B903AA','#D16100','#DDEFFF','#000035','#7B4F4B',
                '#A1C299','#300018','#0AA6D8','#013349','#00846F','#372101','#FFB500','#C2FFED','#A079BF','#CC0744',
                '#C0B9B2','#C2FF99','#001E09','#00489C','#6F0062','#0CBD66','#EEC3FF','#456D75','#B77B68','#7A87A1',
                '#788D66','#885578','#FAD09F','#FF8A9A','#D157A0','#BEC459','#456648','#0086ED','#886F4C','#34362D',
                '#B4A8BD','#00A6AA','#452C2C','#636375','#A3C8C9','#FF913F','#938A81','#575329','#00FECF','#B05B6F',
                '#8CD0FF','#3B9700','#04F757','#C8A1A1','#1E6E00','#7900D7','#A77500','#6367A9','#A05837','#6B002C',
                '#772600','#D790FF','#9B9700','#549E79','#FFF69F','#201625','#72418F','#BC23FF','#99ADC0','#3A2465',
                '#922329','#5B4534','#FDE8DC','#404E55','#0089A3','#CB7E98','#A4E804','#324E72','#6A3A4C','#83AB58',
                '#001C1E','#D1F7CE','#004B28','#C8D0F6','#A3A489','#806C66','#222800','#BF5650','#E83000','#66796D',
                '#DA007C','#FF1A59','#8ADBB4','#1E0200','#5B4E51','#C895C5','#320033','#FF6832','#66E1D3','#CFCDAC',
                '#D0AC94','#7ED379','#012C58','#7A7BFF','#D68E01','#353339','#78AFA1','#FEB2C6','#75797C','#837393',
                '#943A4D','#B5F4FF','#D2DCD5','#9556BD','#6A714A','#001325','#02525F','#0AA3F7','#E98176','#DBD5DD',
                '#5EBCD1','#3D4F44','#7E6405','#02684E','#962B75','#8D8546','#9695C5','#E773CE','#D86A78','#3E89BE',
                '#CA834E','#518A87','#5B113C','#55813B','#E704C4','#00005F','#A97399','#4B8160','#59738A','#FF5DA7',
                '#F7C9BF','#643127','#513A01','#6B94AA','#51A058','#A45B02','#1D1702','#E20027','#E7AB63','#4C6001',
                '#9C6966','#64547B','#97979E','#006A66','#391406','#F4D749','#0045D2','#006C31','#DDB6D0','#7C6571',
                '#9FB2A4','#00D891','#15A08A','#BC65E9','#FFFFFE','#C6DC99','#203B3C','#671190','#6B3A64','#F5E1FF',
                '#FFA0F2','#CCAA35','#374527','#8BB400','#797868','#C6005A','#3B000A','#C86240','#29607C','#402334',
                '#7D5A44','#CCB87C','#B88183','#AA5199','#B5D6C3','#A38469','#9F94F0','#A74571','#B894A6','#71BB8C',
                '#00B433','#789EC9','#6D80BA','#953F00','#5EFF03','#E4FFFC','#1BE177','#BCB1E5','#76912F','#003109',
                '#0060CD','#D20096','#895563','#29201D','#5B3213','#A76F42','#89412E','#1A3A2A','#494B5A','#A88C85',
                '#F4ABAA','#A3F3AB','#00C6C8','#EA8B66','#958A9F','#BDC9D2','#9FA064','#BE4700','#658188','#83A485',
                '#453C23','#47675D','#3A3F00','#061203','#DFFB71','#868E7E','#98D058','#6C8F7D','#D7BFC2','#3C3E6E',
                '#D83D66','#2F5D9B','#6C5E46','#D25B88','#5B656C','#00B57F','#545C46','#866097','#365D25','#252F99',
                '#00CCFF','#674E60','#FC009C','#92896B')

source(paste0(script.dir.graficas, '../instalarPaquetes/instant_pkgs.r'), encoding = 'WINDOWS-1252')
instant_pkgs(c("colorspace", "ggplot2", "ragg", 'reshape'))

crearXYLims <- function(xMin, xMax, yMin, yMax, expand=0) {
  res <- list(xLim=c(xMin, xMax), yLim=c(yMin, yMax), expand=expand)
  if (is.POSIXt(xMin)) res$xLim <- with_tz(time=res$xLim, tzone=tz(xMin))
  return (res)
}

getColoresSetParaN <- function(n) {
  if (n == 1) { return(set1[2])
  } else if(n <= 9) { return(set1[1:n])
  } else if(n <= 12) { return(paired[1:n])
  #} else if(n <= 16) { return(unikn_pair[1:n])
  } else if(n <= 64) { return(colores64[1:n])
  } else { return(colores269[1:n]) }
  #} else if(n <= 269) { return(colores269[1:n]) 
  #} else { stop("graficas.getColoresSetParaN: Can't pick more than 269 colors") }
}

getColoresParaSeries <- function(y, seriesAdicionalesReservadas=0) {
  if (!is.matrix(y)) { n <- 1
  } else { n <- ncol(y) }
  
  res <- getColoresSetParaN(n)
  
  return(c(res, rep('#000000', seriesAdicionalesReservadas)))
}

getTiposDeLineaParaNSeries <- function(n, tiposDistintos=FALSE, seriesAdicionalesReservadas=0) {
  # Imagen con los tipos de línea disponibles en:
  # http://www.sthda.com/english/wiki/ggplot2-line-types-how-to-change-line-types-of-a-graph-in-r-software
  if (tiposDistintos) { res <- c('solid', 'longdash', 'dotdash', 'dotted', 'twodash', 'dashed')[1:n]
  } else { res <- rep('solid', n) }
  return(c(res, rep('blank', seriesAdicionalesReservadas)))
}

getTiposDeLineaParaSeries <- function(y, tiposDistintos=FALSE, seriesAdicionalesReservadas=0) {
  if (!is.matrix(y)) { n <- 1
  } else { n <- ncol(y) }
  
  return(getTiposDeLineaParaNSeries(n = n, tiposDistintos = tiposDistintos, seriesAdicionalesReservadas = seriesAdicionalesReservadas))
}

aplicarOpcionesAGrafico <- function(
    p, titulo='', tituloEjeX='', tituloEjeY='', dibujarEjes=T, xyLims=NULL, tamanioFuenteTextos=15,
    escalaGraficos=1, sinFondo=T) {
  if (sinFondo) p <- p + theme_bw()
  #if (sinFondo) p <- p + theme(panel.background=element_blank())
  p <- p + 
       labs(x=tituloEjeX, y=tituloEjeY, title=titulo) + 
       theme(plot.title = element_text(size=tamanioFuenteTextos * escalaGraficos, hjust = 0.5),
             text = element_text(size=tamanioFuenteTextos * escalaGraficos)) 
  if (!is.null(xyLims)) {
    if (is.null(xyLims$expand)) {
      expand <- 0
    } else {
      expand <- xyLims$expand
    }
    p <- p + coord_cartesian(xlim=xyLims$xLim, ylim=xyLims$yLim, expand=expand)
  }
  if (!dibujarEjes) {
    p <- p + 
      theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())  
  }
  
  return(p)
}

guardarGrafico <- function(p, nomArchSalida, DPI=120, widthPx=800, heightPx=800) {
  path <- dirname(nomArchSalida)
  if (!file.exists(path)) dir.create(path, showWarnings=F, recursive=T)
  ggsave(p, file=nomArchSalida, dpi=DPI, width = widthPx / DPI, height = heightPx / DPI, units = 'in')
}

linePlot <- function(
    x, y, lineaRegresion=F, intervalosConfianza=lineaRegresion, metodoRegresion=lm, 
    formulaRegresion=y~x, dibujarPuntos=F, dibujarLineas=T, titulo='', tituloEjeX='', tituloEjeY='', 
    dibujarEscala=T, dibujarEjes=T, xyLims=NULL, dibujar=interactive(), nomArchSalida=NULL, DPI=120, 
    widthPx=800, heightPx=800, tamanioFuenteTextos=15, grosorLineas=1, escalaGraficos=1, 
    annotateMean=F, colores=getColoresParaSeries(y), tiposDeLinea=getTiposDeLineaParaSeries(y), 
    sinFondo=TRUE, nombresSeries=NULL, nombreVariableGrupo=NULL) {
  df <- data.frame(x, y)
  
  if (is.matrix(y) && ncol(y) > 1) {
    df <- melt(df, id.vars = 'x')
    if (!is.null(nombresSeries)) { df$variable <- factor(nombresSeries[as.integer(df$variable)]) }
    p <- ggplot(df, aes(x=x, y=value)) +
         scale_linetype_manual(values = tiposDeLinea) +
         scale_color_manual(values = colores)
    if (!is.null(nombreVariableGrupo)) p <- p + labs(colour=nombreVariableGrupo, linetype=nombreVariableGrupo)
  } else {
    p <- ggplot(df, aes(x=x, y=y))
  }

  if (dibujarLineas){
    if (is.matrix(y) && ncol(y) > 1) {
      p <- p + geom_line(size=grosorLineas, aes(colour=variable, linetype=variable))
    } else {
      p <- p + geom_line(color=colores[1], linetype=tiposDeLinea[1], size=grosorLineas)
    }
  }
    
  if (dibujarPuntos) p <- p + geom_point(color=colores)
  if (lineaRegresion) p <- p + stat_smooth(method=metodoRegresion, se=intervalosConfianza, formula=formulaRegresion)
  
  if (annotateMean) {
    mediaY <- mean(y, na.rm=T)
    yPos <- max(y, xyLims$yLim[2], na.rm = T)
    p <- p + annotate('text', x = min(x, na.rm = T), y = yPos, label=round(mediaY, 2))
  }
  
  p <- aplicarOpcionesAGrafico(
    p=p, titulo=titulo, tituloEjeX=tituloEjeX, tituloEjeY=tituloEjeY, xyLims=xyLims, 
    dibujarEjes=dibujarEjes, tamanioFuenteTextos=tamanioFuenteTextos, escalaGraficos=escalaGraficos, 
    sinFondo=sinFondo)
  
  if (dibujar) print(p)
  if (!is.null(nomArchSalida)) guardarGrafico(p, nomArchSalida, DPI, widthPx=widthPx, heightPx=heightPx)
  return (p)
}

scatterPlot <- function(
  x, y, lineaRegresion=F, intervalosConfianza=lineaRegresion, metodoRegresion=lm, 
  formulaRegresion=y~x, colorLineaRegresion='red', lineaXIgualY=F, 
  colorLineaXIgualY='black', titulo='', tituloEjeX='', tituloEjeY='', dibujarEscala=T,
  dibujarEjes=T, xyLims=NULL, dibujar=interactive(), nomArchSalida=NULL, DPI=120, 
  widthPx=800, heightPx=800, figurasPuntos=1,
  coloresPuntos=rep('red', length(as.vector(x))), colorMap=NULL,
  tamaniosPuntos=3, tamanioFuenteTextos=15, escalaGraficos=1, sinFondo=TRUE
) {
  if (is.matrix(x) && is.matrix(y) && ncol(x) > 1 && ncol(y) > 1) {
    df <- data.frame(
      variable=c(sapply(colnames(x), function(x, n) rep(x, n), n=nrow(x))), 
      x=as.vector(x), y=as.vector(y), colores=coloresPuntos
    )
    i <- complete.cases(df)
    df <- df[i,]
    multiScatterPlot <- TRUE
  } else {
    auxX <- as.vector(x)
    auxY <- as.vector(y)
    
    i <- !is.na(auxX) & !is.na(auxY)
    df <- data.frame(
      variable=1, x=auxX[i], y=auxY[i], color=coloresPuntos[i]
    )
    
    multiScatterPlot <- FALSE
  }

  if (length(figurasPuntos) > 1) figurasPuntos <- figurasPuntos[i]
  if (length(tamaniosPuntos) > 1) tamaniosPuntos <- tamaniosPuntos[i]
  
  if (multiScatterPlot) {
    colores <- getColoresParaSeries(x)
    
    if (lineaXIgualY) {
      p <- ggplot(df, aes(x=x, y=y, colour=variable)) + 
        geom_abline(slope=1, intercept=0, colour=colorLineaXIgualY, size=escalaGraficos * 0.75) + 
        geom_point(shape=figurasPuntos, size=tamaniosPuntos) + 
        scale_colour_manual(values=colores)
    } else {
      p <- ggplot(df, aes(x=x, y=y, colour=variable)) + 
        geom_point(shape=figurasPuntos, size=tamaniosPuntos) + 
        scale_colour_manual(values=colores)
    }
  } else {
    p <- ggplot(df, aes(x=x, y=y))
    
    if (lineaXIgualY) {
       p <- p + geom_abline(slope=1, intercept=0, colour=colorLineaXIgualY, size=escalaGraficos * 0.75)
    }
      
    if (!is.null(colorMap)) {
      idx <- colorMap$label %in% df$color
      
      values <- setNames(colorMap$color[idx], colorMap$label[idx])
      p <- p + 
        geom_point(aes(colour=color), shape=figurasPuntos, size=tamaniosPuntos) + 
        scale_colour_manual(values=values)
    } else {
      p <- p + geom_point(shape=figurasPuntos, colour=df$color, size=tamaniosPuntos)
    }
  }

  if (lineaRegresion) p <- p + stat_smooth(method=metodoRegresion, se=intervalosConfianza, formula=formulaRegresion, colour=colorLineaRegresion, size=escalaGraficos)
  p <- aplicarOpcionesAGrafico(p = p, titulo = titulo, tituloEjeX = tituloEjeX, tituloEjeY = tituloEjeY, dibujarEjes = dibujarEjes, 
                               xyLims = xyLims, tamanioFuenteTextos = tamanioFuenteTextos, escalaGraficos = escalaGraficos, sinFondo=sinFondo)
  
  if (FALSE) {
    if (is.null(xyLims)) {
      rangoX <- range(df$x)
      rangoY <- range(df$y)
    } else {
      rangoX <- xyLims$xLim
      rangoY <- xyLims$yLim
    }
    
    p <- p + scale_x_continuous(breaks=seq(from = rangoX[1], to = rangoX[2], length.out = 11)) +
             scale_y_continuous(breaks=seq(from = rangoY[1], to = rangoY[2], length.out = 11))
  }
  
  if (dibujar) print(p)
  if (!is.null(nomArchSalida)) guardarGrafico(p, nomArchSalida, DPI, widthPx=widthPx, heightPx=heightPx)
  return (p)
}

boxplot_GGPlot <- function(
    clases, y, titulo='', tituloEjeX='', tituloEjeY='', dibujarEscala=T, dibujarEjes=T, xyLims=NULL, 
    dibujar=interactive(), nomArchSalida=NULL, DPI=120, widthPx=800, heightPx=800, 
    tamanioFuenteTextos=15, escalaGraficos=1, notch=F, dibujarPuntos=F) {
  i <- !is.na(clases) & !is.na(y)
  
  df <- data.frame(classes=factor(clases[i]), values=y[i])
  p <- ggplot(df, aes(classes, values)) + geom_boxplot(aes(group=classes), notch=notch)
  if (dibujarPuntos) p <- p + geom_jitter(width = 0.2)
  p <- aplicarOpcionesAGrafico(
    p=p, titulo=titulo, tituloEjeX=tituloEjeX, tituloEjeY=tituloEjeY, xyLims=xyLims, 
    tamanioFuenteTextos=tamanioFuenteTextos, escalaGraficos=escalaGraficos)
  if (dibujar) print(p)
  if (!is.null(nomArchSalida)) guardarGrafico(p, nomArchSalida, DPI, widthPx=widthPx, heightPx=heightPx)
  return (p)
}

graficoCorrVsDistancia <- function(
  dist, corr, clasesEstaciones=rep(1, nrow(corr)), defaultClass=NULL, nomArchSalida=NULL,
  tamaniosPuntos=2, figurasPuntos=20, widthPx=1600, heightPx=900, DPI=120
) {
  upperTri <- upper.tri(dist, diag=T)
  
  df <- data.frame(x=dist[upperTri], y=corr[upperTri])
  formulaModelo <- y~Matern(d=x, range=range, smoothness=smoothness)
  modeloMatern <- nls(formulaModelo, data=df, start=list(smoothness=0.3, range=400))
  params <- coef(modeloMatern)
  smoothness <- params[1]
  range <- params[2]
  formulaModeloAjustada <- as.formula(formulaModelo, env = new.env(range, smoothness))
  
  clasesUnique <- unique(clasesEstaciones)
  if (length(clasesUnique) > 1) {
    n <- (length(clasesUnique) * (length(clasesUnique) + 1)) / 2
    colorMap <- data.frame(label=rep(NA_character_, n), color=rep(NA_character_, n))
    coloresPuntos <- matrix(NA_character_, nrow=nrow(corr), ncol=ncol(corr))
    n <- 1
    for (i in seq_along(clasesUnique)) {
      for (j in seq.int(i, length(clasesUnique))) {
        idxClase1 <- clasesEstaciones == clasesUnique[i]
        idxClase2 <- clasesEstaciones == clasesUnique[j]
        
        if (!is.null(defaultClass) && 
            (clasesUnique[i] == clasesUnique[j] & sum(idxClase1) == 1)
        ) {
          coloresPuntos[idxClase1, idxClase2] <- defaultClass
          coloresPuntos[idxClase2, idxClase1] <- defaultClass
        } else {
          if (i == j) {
            colorMap[n, 'label'] <- clasesUnique[i]
          } else {
            colorMap[n, 'label'] <- paste0(clasesUnique[i], ',', clasesUnique[j])  
          }
                    
          coloresPuntos[idxClase1, idxClase2] <- colorMap[n, 'label']
          coloresPuntos[idxClase2, idxClase1] <- colorMap[n, 'label']
    }
        n <- n + 1
      }
    }
    colorMap <- colorMap[!is.na(colorMap$label), ]
    colorMap$color <- getColoresSetParaN(n=nrow(colorMap))
    rownames(colorMap) <- 1:nrow(colorMap)
    
    # which(clasesEstaciones!='Ok')
    # clasesEstaciones[clasesEstaciones!='Ok']
    # unique(coloresPuntos[upperTri]) %in% colorMap$label
  } else {
    coloresPuntos <- rep('red', prod(dim(dist)))
    colorMap <- NULL
  }
  
  xyLimsScatterPlot <- crearXYLims(
    yMin=achicarToNDigitos(min(corr[upperTri], na.rm=T), 2), yMax=1,
    xMin=0, xMax=agrandarToNDigitos(max(dist[upperTri]), 1)
  )
  
  # TODO: agregar labels de los puntos. Idea interesante acá:
  # https://stackoverflow.com/questions/52009545/r-ggplot-for-hover-text
  return(scatterPlot(
    x=dist[upperTri], y=corr[upperTri], lineaRegresion=T, intervalosConfianza=F, 
                     formulaRegresion=formulaModeloAjustada, xyLims=xyLimsScatterPlot,
    coloresPuntos=coloresPuntos[upperTri], colorMap=colorMap, 
    tituloEjeX='Distancia[km]', tituloEjeY='Correlación',
    dibujar=F, nomArchSalida=nomArchSalida, tamaniosPuntos=tamaniosPuntos,
    figurasPuntos=figurasPuntos, widthPx=widthPx, heightPx=heightPx, DPI=DPI
  ))
}
