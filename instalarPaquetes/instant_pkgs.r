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

# Los writes comentados son para debuggear el bug de library en R 3.4.0 cuando se usan procesos
#write(getwd(), paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)

cargarPaquetes <- function(pkgs, silent=T, nRetries=25) {
  #write('CP1', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
  if (length(pkgs) > 0) {
    #write('CP2', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    attached <- search()
    attached_pkgs <- attached[grepl("package", attached)]
    need_to_attach <- pkgs[!pkgs %in% gsub("package:", "", attached_pkgs)]
    #write('CP3', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    
    #write(need_to_attach, paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    if (silent) {
      # i <- seq_along(need_to_attach)[1]
      for (i in seq_along(need_to_attach)) { 
        nIntentos <- 0
        while (!try(suppressPackageStartupMessages(library(need_to_attach[i], character.only=T, logical.return = T, quietly = silent))) &
               nIntentos < nRetries) {
          Sys.sleep(1)
          nIntentos <- nIntentos + 1
        }
        #if (nIntentos >= nRetries)
        #  write(need_to_attach[i], paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
      }
    } else {
      for (i in seq_along(need_to_attach)) { 
        nIntentos <- 0
        while (!try(library(need_to_attach[i], character.only=T, logical.return = T, quietly = silent)) &
               nIntentos < nRetries) {
          Sys.sleep(1)
          nIntentos <- nIntentos + 1
        }
        #if (nIntentos >= nRetries)
        #  write(need_to_attach[i], paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
      }
    }
    # write('CP4', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
  }
  return(NULL)
}

instant_pkgs <- function(pkgs, minVersions=rep(NA_character_, length(pkgs)), silent=TRUE, doCargarPaquetes=TRUE) {
  if (length(pkgs) > 0) {
    # write('IP1', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    # write(pkgs, paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    
    paquetesInstalados <- installed.packages()
    iInstalados <- pkgs %in% rownames(paquetesInstalados)
    iInstalados[iInstalados] <- is.na(minVersions[iInstalados]) | 
      paquetesInstalados[pkgs[iInstalados], 'Version'] >= minVersions[iInstalados]
    paquetesAInstalar <- pkgs[!iInstalados]
    # write('IP2', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)

    if (length(paquetesAInstalar) > 0) {
      # write('IP3', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
      if (.Platform$OS.type == "windows") {
        utils::install.packages(paquetesAInstalar, repos = 'https://cran.rstudio.com/', dependencies = TRUE, type='binary')
      } else {
        utils::install.packages(paquetesAInstalar, repos = 'https://cran.rstudio.com/', dependencies = TRUE)
      }
      # write('IP4', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    }
    if (doCargarPaquetes) cargarPaquetes(pkgs, silent = silent)
    # write('IP5', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
  }
  return(NULL)
}

instant_pkgs_github <- function(reposgithub, pkgs=basename(reposgithub), minVersions, silent=T) {
  if (length(pkgs) > 0) {
    # write('IPG1', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    # write(pkgs, paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    paquetesInstalados <- installed.packages()
    iInstalados <- pkgs %in% rownames(paquetesInstalados)
    iInstalados[iInstalados] <- is.na(minVersions[iInstalados]) | 
      paquetesInstalados[pkgs[iInstalados], 'Version'] >= minVersions[iInstalados]
    paquetesAInstalar <- pkgs[!iInstalados]
    reposgithub <- reposgithub[!iInstalados]
    # write('IPG2', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    # i <- seq_along(reposgithub)[1]
    if (length(reposgithub) > 0) {
      # write('IPG3', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
      # devtools tiene algún problema con la carga en los threads. Trato de limitar su carga lo más posible
      instant_pkgs('devtools')
      # write('IPG4', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
      install_github(repo = reposgithub, quiet = FALSE)
      # write('IPG5', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    }
    
    cargarPaquetes(pkgs, silent = silent)
    # write('IPG6', paste('D:/testsMCH/instantPkgs/', Sys.getpid(), '.txt', sep=''), append = T)
    # flock::unlock(file.lock = lock)
  }
  return(NULL)
}

updateAllPkgs <- function() {
  update.packages(lib.loc = .libPaths()[1])
}

#instant_pkgs(pkgs = 'flock', silent = T)

detachAllPackages <- function() {
  #Util para debugear. Descarga todos los paquetes  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}