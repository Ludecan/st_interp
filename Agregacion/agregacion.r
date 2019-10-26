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
if (iFrame >= 3) { script.dir.agregacion <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.agregacion <- NULL }
while ((is.null(script.dir.agregacion) || is.na(regexpr('agregacion.r', script.dir.agregacion, fixed=T)[1])) && iFrame >= 0) {
  script.dir.agregacion <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.agregacion)) { script.dir.agregacion <- ''
} else { script.dir.agregacion <- paste(dirname(script.dir.agregacion), '/', sep='') }

source(paste(script.dir.agregacion, '../instalarPaquetes/instant_pkgs.r', sep=''))
instant_pkgs(c('stats', 'sp'))

naSiTodosNAFuncSiNo <- function(x, func, ...) {
  x <- x[!is.na(x)]
  if (length(x)==0) { return (NA)
  } else { return (func(x, ...)) }
}

naSiPorcNAsFuncSiNo <- function(x, func, porcDatosNoNulos=0.75, ...) {
  iNoNA <- !is.na(x)
  if (sum(iNoNA) > length(x) * porcDatosNoNulos) { return(func(x[iNoNA], ...))
  } else { return(NA) }
}

ocurrenciaUmbral <- function(x, umbral=0, comparador='=') {
  # Retorna 1 si existe un valor que cumpla con umbral la relacion definida por comparador 
  # Retorna 0 sino
  i <- !is.na(x)
  if (any(i)) { 
    if (comparador == '=') { return (any(x[i] == umbral))
    } else if (comparador == '>') { return (any(x[i] > umbral))
    } else if (comparador == '<') { return (any(x[i] < umbral))
    } else if (comparador == '>=') { return (any(x[i] >= umbral))
    } else if (comparador == '<=') { return (any(x[i] <= umbral))
    } else { stop(paste('Comparador desconocido', comparador))
    }
  } else { return(FALSE) }
}

contarUmbral <- function(x, umbral=0, comparador='=') {
  # Retorna la cantidad de valores que cumplan con umbral la relacion definida por comparador 
  i <- !is.na(x)
  if (any(i)) { 
    if (comparador == '=') { return (sum(x[i] == umbral))
    } else if (comparador == '>') { return (sum(x[i] > umbral))
    } else if (comparador == '<') { return (sum(x[i] < umbral))
    } else if (comparador == '>=') { return (sum(x[i] >= umbral))
    } else if (comparador == '<=') { return (sum(x[i] <= umbral))
    } else { stop(paste('Comparador desconocido', comparador))
    }
  } else { return(0) }
}

porcentajeUmbral <- function(x, umbral=0, comparador='=') {
  # Retorna el porcentaje de valores que cumplan con umbral la relacion definida por comparador 
  i <- !is.na(x)
  if (any(i)) { 
    if (comparador == '=') { return (sum(x[i] == umbral) / sum(i) * 100)
    } else if (comparador == '>') { return (sum(x[i] > umbral) / sum(i) * 100)
    } else if (comparador == '<') { return (sum(x[i] < umbral) / sum(i) * 100)
    } else if (comparador == '>=') { return (sum(x[i] >= umbral) / sum(i) * 100)
    } else if (comparador == '<=') { return (sum(x[i] <= umbral) / sum(i) * 100)
    } else { stop(paste('Comparador desconocido', comparador))
    }
  } else { return(FALSE) }
}

parseFuncionAgregacion <- function(funcionAgregacion) {
  #x <- seq(1, 10, length.out=10)
  #funcionAgregacion <- 'max'
  #funcionAgregacion <- 'ocurrenciaUmbral;umbral=1;comparador=='
  #lala <- parseFuncionAgregacion(funcionAgregacion)
  #formals(lala)
  #lala(x)  
  
  paramsSplit <- unlist(strsplit(funcionAgregacion, ";"))
  funcion <- get(paramsSplit[1])
  paramsFuncion <- formals(funcion)
  if (!is.null(paramsFuncion)) {
    paramsSplit <- paramsSplit[-1]
    for (i in 1:length(paramsSplit)) {
      posIgual <- regexpr('=', paramsSplit[i], fixed=T)[1]
      nomParam <- substr(x=paramsSplit[i], start=1, stop=posIgual - 1)
      valParam <- type.convert(substr(x=paramsSplit[i], start=posIgual+1, stop=nchar(paramsSplit[i])), as.is=T)
      
      formals(funcion)[nomParam] <- valParam
    }
  }
  return (funcion)
}

agregar <- function(x, funcionAgregacion, claseIndiceI=rep(1, nrow(x)), ordenarPorClases=F, reducirDimensionesUnitarias=T, na.rm=T, ...) {
  # reduce la primer dimension de x de acuerdo a las clases en claseIndiceI usando funcionAgregacion
  # retorna una matriz length(unique(claseIndiceI)) x dim(x)[2:length(dim(x))] o length(unique(claseIndiceI)) x 1 si x es un vector
  # por ejemplo si x = (1, 2, 3, 4, 5, 6), claseIndiceI=(1,1, 2,2, 3,3) y funcionAgregacion=max
  # retorna una matriz de 3x1 = (2, 4, 6)
  # si x = (1, 2, 
  #         3, 4, 
  #         5, 6), claseIndiceI=(1, 1, 2) y funcionAgregacion=max 
  # retorna una matriz de 2x2 = (3, 4, 
  #                              5, 6)
  # el valor por defecto de claseIndiceI reduce la primer dimension a un único valor

  #funcionAgregacion <- max
  #claseIndiceI <- c(2,2, 1, 1, 1)
  if (is.character(funcionAgregacion))
    funcionAgregacion <- parseFuncionAgregacion(funcionAgregacion)
  
  if (na.rm) { res <- aggregate(x, by=list(claseIndiceI), FUN=naSiTodosNAFuncSiNo, func=funcionAgregacion, ...)
  } else { res <- aggregate(x, by=list(claseIndiceI), FUN=funcionAgregacion, ...) }
  res <- res[,-1]
  if (ordenarPorClases) {
    clases <- unique(claseIndiceI)
    res[order(clases),]
  }

  if (!is.null(dim(x))) {
    dimRes <- c(length(unique(claseIndiceI)), dim(x)[2:length(dim(x))])
    if (reducirDimensionesUnitarias) dimRes <- dimRes[dimRes > 1]

    if (length(dimRes) == 2) { res <- matrix(unlist(res), nrow=dimRes[1], ncol=dimRes[2])
    } else { res <- array(unlist(res), dim=dimRes) }
    
    if (!is.null(colnames(x)) & ncol(x) == ncol(res)) colnames(res) <- colnames(x)
  }
  
  return(res)
}

agregacionEspacialAPoligonos <- function(spObj, shpPoligonos, funcionAgregacion, zcol=1, na.rm=T) {
  if (is.character(funcionAgregacion))
    funcionAgregacion <- parseFuncionAgregacion(funcionAgregacion)

  if (na.rm) {
    return (as.numeric(over(y=spObj[,zcol], x=shpPoligonos, fn=naSiTodosNAFuncSiNo, func=funcionAgregacion)[,1]))
  } else {
    return(as.numeric(over(y=spObj[,zcol], x=shpPoligonos, fn=funcionAgregacion)[,1]))
  }
}

testAgregacion <- function() {
  funcionAgregacion <- max
  
  x1 <- c(1, 2, NA, NA, 5, 6)
  claseIndiceI1 <- c(1,1, 2,2, 3,3)
  x1
  claseIndiceI1
  agregar(x1, funcionAgregacion, claseIndiceI1, na.rm=T)
  
  x2 <- matrix(data=c(1, 2, NA, 4, NA, 6), nrow=3, ncol=2, byrow=T)
  claseIndiceI2 <- c(1, 1, 2)
  x2
  claseIndiceI2
  agregar(x2, funcionAgregacion, claseIndiceI2, na.rm=T)
  
  x3 <- array(c(rnorm(30)), dim=c(3, 5, 2), dimnames=c('T', 'Y', 'X'))
  x3[1,1,1] <- NA
  x3[2,1,1] <- NA
  x3[3,1,1] <- NA
  claseIndiceI3 <- c(1, 1, 2)
  x3[1,,]
  x3[2,,]
  x3[3,,]
  claseIndiceI3
  res <- agregar(x3, funcionAgregacion, claseIndiceI3, na.rm=F)
  res[1,,]
  res[2,,] == x3[3,,]
  
  xdf <- data.frame(x=x1, y=x1+1, z=x1*2)
  claseIndiceIdf <- c(1, 2, 1, 2, 3, 3)
  xdf
  claseIndiceIdf
  agregar(xdf, funcionAgregacion, claseIndiceIdf)
}