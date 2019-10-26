script.dir.PathUtils <- dirname((function() { attr(body(sys.function()), "srcfile") })()$filename)

# source(paste(script.dir.PathUtils, '/../instalarPaquetes/instant_pkgs.r', sep=''))
# instant_pkgs("tools")
# funciones útiles que ya vienen en R
# basename, dirname, dir.create

nombreArchSinExtension <- function(filename) {
  #             1234567890123456
  #filename <- 'c:/lele.lala.txt'
  #filename <- 'c:/lele'
  #filename <- '//192.168.0.11/wrfout_d02_2014-02-15_00%3A00%3A00'
  iUltimaBarra <- max(gregexpr(pattern = '/', filename)[[1]])
  iUltimoPunto <- max(gregexpr(pattern = '\\.', filename)[[1]])
  # lo que hay después del último punto es una extensión solo si el último punto
  # está después de la última barra
  if (iUltimoPunto >= 1 & iUltimaBarra < iUltimoPunto) { 
    return(substr(filename, start=1, stop=iUltimoPunto-1))
  } else { return(filename) }
}

nombreArchSinPathNiExtension <- function(filename) {
  return (nombreArchSinExtension(basename(filename)))
}

# devuelve la extensión sin el punto
getFileExt <- function(filename) {
  iUltimaBarra <- sapply(gregexpr(pattern = '/', filename), FUN = max)
  iUltimoPunto <- sapply(gregexpr(pattern = '\\.', filename), FUN = max)
  
  # lo que hay después del último punto es una extensión solo si el último punto
  # está después de la última barra
  res <- character(length(filename))
  for (i in 1:length(filename)) {
    if (!is.na(iUltimoPunto[i]) & iUltimoPunto[i] >= 1 & iUltimaBarra[i] < iUltimoPunto[i]) {
      res[i] <- substr(filename[i], start=iUltimoPunto[i]+1, stop=nchar(filename[i]))
    } else { res[i] <- '' }
  }
  return(res)
}

changeFileExt <- function(filename, nuevaExtensionConPunto) {
  return(paste(nombreArchSinExtension(filename), nuevaExtensionConPunto, sep=''))
}

changeFileDrive <- function(filename, nuevoDriveSinDosPuntos) {
  posDosPuntos <- min(gregexpr(pattern = ':', filename)[[1]])
  if (posDosPuntos == 2) {
    return(paste(nuevoDriveSinDosPuntos, substr(filename, start=2, stop=nchar(filename)), sep=''))
  } else {
    return(filename)
  }
}

prependToFileName <- function(filename, preFijo) {
  return(sapply(filename, function(x) {
    nomArch <- basename(x)
    paste(gsub(x = x, pattern = nomArch, replacement = ''), preFijo, nomArch, sep='')  
  }))
}

appendToFileName <- function(filename, postFijo) {
  ext <- getFileExt(filename)
  return(ifelse(ext != "", 
                yes = paste(substr(filename, start=1, stop=nchar(filename) - (nchar(ext) + 1)), postFijo, '.', ext, sep=''),
                no = paste(filename, postFijo, sep='')))
}

agregarCarpetaAlFinal <- function(filename, carpeta) {
  if ((carpeta != '') && (carpeta != './')) {
    if (grepl(pattern = '/$', x = carpeta)) {
      return(paste(dirname(filename), '/', carpeta, basename(filename), sep=''))  
    } else { return(paste(dirname(filename), '/', carpeta, '/', basename(filename), sep=''))
    }
  } else { return(filename) }
}

crearDirectoriosSiNoExisten <- function(directorio) {
  if (!file.existsEx(directorio))
    dir.create(directorio, recursive=T, showWarnings=F)
}

esRutaAbsoluta <- function(filename) {
  segundoChar <- substr(filename, start=2, stop=2)
  return (ifelse(segundoChar == ':', TRUE, ifelse(segundoChar == '/', substr(filename, start=1, stop=1) == '/', FALSE)))
}

pathParaGuardadoDeArchivos <- function(filename, pathProceso) {
  return (ifelse(esRutaAbsoluta(filename), yes = filename, agregarCarpetaAlFinal(filename, carpeta = pathProceso)))
}

removeTrailingSlash <- function(paths) {
  ultimochar <- substr(paths, start=nchar(paths), stop=nchar(paths))
  res <- sapply(paths, 
                FUN = function(x) {
                  ultimochar <- substr(x, start=nchar(x), stop=nchar(x))
                  if (ultimochar == '/') { return(substr(x, start=1, stop=nchar(x)-1))
                  } else { return(x) }
                })
  return(res)
}

file.existsEx <- function(path) {
  return(file.exists(removeTrailingSlash(path)))
}

subirDirectorio <- function(path) {
  iUltimaBarra <- max(gregexpr(pattern = '/', removeTrailingSlash(path))[[1]])
  if (length(iUltimaBarra) > 0) { return(substr(path, start=1, stop=iUltimaBarra))
  } else { stop(paste('En el path "', path, '" no hay carpeta anterior', sep ='')) }
}

ultimaCarpeta <- function(paths, removerTrailingSlash=T) {
  if (removerTrailingSlash) { iBarras <- gregexpr(pattern = '/', removeTrailingSlash(paths))
  } else { iBarras <- gregexpr(pattern = '/', paths) }
  res <- lapply(1:length(paths), 
                FUN = function(x, paths, iBarras) {
                  if (length(iBarras[x]) > 0) { 
                    return(substr(paths[x], start=iBarras[[x]][length(iBarras[[x]])-1]+1, stop=iBarras[[x]][length(iBarras[[x]])]-1))
                  } else { 
                    return(NA)
                  }
                }, paths, iBarras)
  return(res)
}

pathDesdeUltimaCarpeta <- function(paths) {
  iBarras <- gregexpr(pattern = '/', removeTrailingSlash(paths))
  res <- sapply(1:length(paths), 
                FUN = function(x, paths, iBarras) {
                        if (length(iBarras[x]) > 0) { iIni <- iBarras[[x]][length(iBarras[[x]])-1]+1
                        } else { iIni <- 1 }
                        return(substr(paths[x], start=iIni, stop=nchar(paths[x])))
                }, paths, iBarras)
  return(res)
}
