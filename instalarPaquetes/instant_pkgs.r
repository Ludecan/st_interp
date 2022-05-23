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
if (iFrame >= 3) { script.dir.instantPkgs <- sys.frame(iFrame - 3)$ofile 
} else { script.dir.instantPkgs <- NULL }
while ((is.null(script.dir.instantPkgs) || is.na(regexpr('instant_pkgs.r', script.dir.instantPkgs, fixed=T)[1])) && iFrame >= 0) {
  script.dir.instantPkgs <- sys.frame(iFrame)$ofile
  iFrame <- iFrame - 1
}
if (is.null(script.dir.instantPkgs)) { script.dir.instantPkgs <- ''
} else { script.dir.instantPkgs <- paste0(dirname(script.dir.instantPkgs), '/') }

cargarPaquetes <- function(pkgs, silent=T, nRetries=25) {
  if (length(pkgs) > 0) {
    attached <- search()
    attached_pkgs <- attached[grepl("package", attached)]
    need_to_attach <- pkgs[!pkgs %in% gsub("package:", "", attached_pkgs)]
    
    if (silent) {
      for (i in seq_along(need_to_attach)) { 
        nIntentos <- 0
        while (!try(suppressPackageStartupMessages(library(need_to_attach[i], character.only=T, logical.return = T, quietly = silent))) &
               nIntentos < nRetries) {
          Sys.sleep(1)
          nIntentos <- nIntentos + 1
        }
      }
    } else {
      for (i in seq_along(need_to_attach)) { 
        nIntentos <- 0
        while (!try(library(need_to_attach[i], character.only=T, logical.return = T, quietly = silent)) &
               nIntentos < nRetries) {
          Sys.sleep(1)
          nIntentos <- nIntentos + 1
        }
      }
    }
  }
  return(NULL)
}

checkInstallPackage <- function(i, pkgs, minVersions, paquetesInstalados) {
  return(
    !pkgs[i] %in% rownames(paquetesInstalados)
    || (
      !is.na(minVersions[i])
      && utils::compareVersion(paquetesInstalados[pkgs[i], 'Version'], minVersions[i]) < 0
    )
  )
}

checkInstallPackages <- function(pkgs, minVersions) {
  paquetesInstalados <- installed.packages()
  paquetesInstalados <- paquetesInstalados[rownames(paquetesInstalados) %in% pkgs, , drop=F]
  bPaquetesAInstalar <- sapply(
    seq_along(pkgs), checkInstallPackage, pkgs, minVersions, paquetesInstalados)
  return(bPaquetesAInstalar)
}

instant_pkgs <- function(pkgs, minVersions=rep(NA_character_, length(pkgs)), silent=TRUE, doCargarPaquetes=TRUE) {
  if (length(pkgs) > 0) {
    bPaquetesAInstalar <- checkInstallPackages(pkgs, minVersions)
    paquetesAInstalar <- pkgs[bPaquetesAInstalar]

    if (length(paquetesAInstalar) > 0) {
      Ncpus <- parallel::detectCores(T, T)
      
      tryCatch(
        expr = {
          if (.Platform$OS.type == "windows") {
            utils::install.packages(
              paquetesAInstalar, repos='https://cran.rstudio.com/', dependencies=TRUE, type='binary', Ncpus=Ncpus)
          } else {
            utils::install.packages(
              paquetesAInstalar, repos='https://cran.rstudio.com/', dependencies=TRUE, Ncpus=Ncpus)
          }
        },
        except={
          if (.Platform$OS.type == "windows") {
            utils::install.packages(
              paquetesAInstalar, repos='https://cran.rstudio.com/', dependencies=TRUE, type='binary')
          } else {
            utils::install.packages(
              paquetesAInstalar, repos='https://cran.rstudio.com/', dependencies=TRUE)
          }
        }
      )
    }
    if (doCargarPaquetes) cargarPaquetes(pkgs, silent = silent)
  }
  return(NULL)
}

instant_pkgs_github <- function(reposgithub, pkgs=basename(reposgithub), minVersions, silent=T) {
  if (length(pkgs) > 0) {
    bPaquetesAInstalar <- checkInstallPackages(pkgs, minVersions)
    paquetesAInstalar <- pkgs[bPaquetesAInstalar]
    reposgithub <- reposgithub[bPaquetesAInstalar]
    if (length(reposgithub) > 0) {
      # devtools tiene algún problema con la carga en los threads. Trato de limitar su carga lo más posible
      instant_pkgs('devtools')
      install_github(repo = reposgithub, quiet = FALSE)
    }
    
    cargarPaquetes(pkgs, silent = silent)
  }
  return(NULL)
}

updateAllPkgs <- function() {
  update.packages(lib.loc = .libPaths()[1])
}

detachAllPackages <- function() {
  #Util para debugear. Descarga todos los paquetes  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

checkInstalledPackages <- function(minPkgVersionsPath) {
  minPkgVersions <- read.table(
    minPkgVersionsPath, header=TRUE, sep='\t', encoding='UTF-8', stringsAsFactors=FALSE
  )
  installedPackages <- installed.packages()
  
  toUpdatePkgs <- list()
  toUpdatePkgVersions <- list()
  
  i_pkg <- 1
  i_pkg <- which(installedPackages[, 'Package'] == 'spatial')
  for (i_pkg in seq.int(from=1, to=nrow(installedPackages))) {
    pkg <- installedPackages[i_pkg, 'Package']
    priority <- installedPackages[i_pkg, 'Priority']
    if (is.na(priority) || priority != 'base') {
      version <- installedPackages[i_pkg, 'Version']
      
      i_minPkg = which(minPkgVersions[, 'Package'] == pkg)
      
      if (length(i_minPkg) > 0) {
        minVersion <- minPkgVersions[i_minPkg[1], 'Version']
        
        if (utils::compareVersion(version, minVersion) == -1) {
          print(paste0('Updating ', pkg, ' from ', version, ' to ', minVersion))
          
          n <- length(toUpdatePkgs) + 1
          toUpdatePkgs[n] <- pkg
          toUpdatePkgVersions[n] <- minVersion
        } # else {
        #print(paste0('Skipping ', pkg, ' as installed version ', version, ' >= ', minVersion))
        #}
      }      
    }
  }

  if (length(toUpdatePkgs)) {
    toUpdatePkgs <- unlist(toUpdatePkgs)
    toUpdatePkgVersions <- unlist(toUpdatePkgVersions)
    idx <- order(toUpdatePkgs)
    toUpdatePkgs <- toUpdatePkgs[idx]
    toUpdatePkgVersions <- toUpdatePkgVersions[idx]
    
    toUpdatePkgVersions <- toUpdatePkgVersions[!duplicated(toUpdatePkgs)]
    toUpdatePkgs <- toUpdatePkgs[!duplicated(toUpdatePkgs)]
    
  
    instant_pkgs(pkgs=toUpdatePkgs, minVersions=toUpdatePkgVersions, doCargarPaquetes=FALSE)
  }
}

minPkgVersionsPath <- paste0(script.dir.instantPkgs, 'min_pkg_versions.tsv')
if (!exists('installedPackagesChecked') && file.exists(minPkgVersionsPath)) {
  checkInstalledPackages(minPkgVersionsPath)
}
installedPackagesChecked <- TRUE

createMinPackageVersions <- function() {
  packages <-installed.packages()
  write.table(x = packages, file = minPkgVersionsPath, sep='\t', row.names = T, col.names = T, 
              fileEncoding = 'UTF-8')
}
#renv::snapshot()
