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

require(rgdal)
mod11a1 <- readGDAL('//192.168.1.223/mch/datosEnGrilla/fdg/MOD11A1/MOD11A1_LST_Day/MOD11A1_2011-12-25.LST_Day_1km.tif')
# nX
length(unique(sp::coordinates(mod11a1)[,1]))
# nY
length(unique(sp::coordinates(mod11a1)[,2]))
length(unique(sp::coordinates(mod11a1)[,1]))*length(unique(coordinates(mod11a1)[,2]))
bbox(mod11a1)

grillaUy <- readGDAL('//192.168.1.223/mch/ReportesWeb/Mapas/TempAire/2016_03_18_18_00_00_TempAire_Automaticas_Horaria.tif')

# nX
length(unique(sp::coordinates(grillaUy)[,1]))
# nY
length(unique(sp::coordinates(grillaUy)[,2]))
length(unique(sp::coordinates(grillaUy)[,1]))*length(unique(coordinates(grillaUy)[,2]))

source('D:/Workspace/MCH2/MCH/MCHLib/R/Scripts/interpolar/interpolarEx.r', encoding = 'WINDOWS-1252')
source('D:/Workspace/MCH2/MCH/MCHLib/R/Scripts/interpolar/mapearEx.r', encoding = 'WINDOWS-1252')

shp <- cargarSHP('C:/mch/ArchivosProcesosLocales/MapaUruguayVacio/uruguay_departamentos.shp', '+proj=longlat +datum=WGS84')

proj4string(shp)
proj4string(mod11a1)
proj4string(grillaUy)

bbox(shp)
bbox(mod11a1)
bbox(grillaUy)

bbox(spTransform(shp, mod11a1@proj4string))
lele <- spTransform(shp, mod11a1@proj4string)

bbox(spTransform(mod11a1, shp@proj4string))
bbox(spTransform(grillaUy, shp@proj4string))

LatitudMaxPais <-	-30
LatitudMinPais <-	-35
LongitudMaxPais <- -53.12
LongitudMinPais <- -58.5

xLim <- bbox(shp)[1, ]
yLim <- bbox(shp)[2, ]
options(digits = 15)
xLim
yLim
nDigitos <- 2
c(achicarToNDigitos(xLim[1], nDigitos), agrandarToNDigitos(xLim[2], nDigitos))
c(achicarToNDigitos(yLim[1], nDigitos), agrandarToNDigitos(yLim[2], nDigitos))

# Coordenadas del país
xLim <- c(-58.45, -53.15)
yLim <- c(-35, -30)
p4strLatLong <- '+proj=longlat +datum=WGS84 +no_defs'
SRS_stringLatLong <- "EPSG:4326"
p4strPlana <- '+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs'
SRS_stringPlana <- "EPSG:32721"

# Creo los 4 vertices que encierran al país
coordsBB <- matrix(data=c(xLim[1], yLim[1], #abajo, izquierda
                          xLim[1], yLim[2], #arriba, izquierda
                          xLim[2], yLim[2], #arriba, derecha
                          xLim[2], yLim[1]), #abajo, derecha
                   ncol=2, byrow=T)
sps <- SpatialPoints(coords = coordsBB, proj4string = CRS(projargs = p4strLatLong, SRS_string = SRS_stringLatLong))
# Los proyecto para hallar el rango de coordenadas en la proyección plana
sps <- spTransform(sps, CRS(projargs = p4strPlana, SRS_string = SRS_stringPlana))
# Hallo los rangos de coordenadas en la proyección plana
xLim <- bbox(sps)[1,]
yLim <- bbox(sps)[2,]

xDim <- 1000
yDim <- 1000
grid <- GridTopology(c(xLim[1], yLim[1]), c(xDim, yDim), cells.dim = c(ceiling((xLim[2]-xLim[1])/xDim), ceiling((yLim[2]-yLim[1])/yDim)))
grid <- SpatialGrid(grid = grid, proj4string = CRS(projargs = p4strPlana, SRS_string = SRS_stringPlana))

(xLim[2]-xLim[1])/ceiling((xLim[2]-xLim[1])/xDim)
(yLim[2]-yLim[1])/ceiling((yLim[2]-yLim[1])/yDim)

data <- data.frame(as.integer(!is.na(over(grid, geometry(mod11a1)))))
grid <- SpatialGridDataFrame(grid, data, proj4string = CRS(projargs = p4strPlana, SRS_string = SRS_stringPlana))

shp <- spTransform(shp, CRS(projargs = p4strPlana, SRS_string = SRS_stringPlana))
mapearGrillaGGPlot(grilla = grid, shpBase = shp)

geometry(mod11a1)@grid
geometry(grillaUy)@grid
geometry(grid)@grid
nrow(mod11a1)
nrow(grid)
nrow(grillaUy)
