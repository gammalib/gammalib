#############################################################################
# Makefile.am for model module                                              #
# ------------------------------------------------------------------------- #
# Copyright (C) 2016-2021 Juergen Knoedlseder                               #
# ------------------------------------------------------------------------- #
#                                                                           #
#  This program is free software: you can redistribute it and/or modify     #
#  it under the terms of the GNU General Public License as published by     #
#  the Free Software Foundation, either version 3 of the License, or        #
#  (at your option) any later version.                                      #
#                                                                           #
#  This program is distributed in the hope that it will be useful,          #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#  GNU General Public License for more details.                             #
#                                                                           #
#  You should have received a copy of the GNU General Public License        #
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                           #
#############################################################################

# Set include directories
INCLUDES = -I$(top_srcdir)/include \
           -I$(top_builddir)/src \
           -I$(top_srcdir)/src

# Define sources for this directory
sources = GModelPar.cpp \
          GModelAssociation.cpp \
          GModelAssociations.cpp \
          GModels.cpp \
          GModel.cpp \
          GModelRegistry.cpp \
          GModelSky.cpp \
          GModelData.cpp \
          GModelSpatial.cpp \
          GModelSpatialRegistry.cpp \
          GModelSpatialPointSource.cpp \
          GModelSpatialRadial.cpp \
          GModelSpatialRadialDisk.cpp \
          GModelSpatialRadialRing.cpp \
          GModelSpatialRadialGauss.cpp \
          GModelSpatialRadialShell.cpp \
          GModelSpatialRadialProfile.cpp \
          GModelSpatialRadialProfileGauss.cpp \
          GModelSpatialRadialProfileDMZhao.cpp \
          GModelSpatialRadialProfileDMBurkert.cpp \
          GModelSpatialRadialProfileDMEinasto.cpp \
          GModelSpatialElliptical.cpp \
          GModelSpatialEllipticalDisk.cpp \
          GModelSpatialEllipticalGauss.cpp \
          GModelSpatialDiffuse.cpp \
          GModelSpatialDiffuseConst.cpp \
          GModelSpatialDiffuseCube.cpp \
          GModelSpatialDiffuseMap.cpp \
          GModelSpatialComposite.cpp \
          GModelSpectral.cpp \
          GModelSpectralRegistry.cpp \
          GModelSpectralBrokenPlaw.cpp \
          GModelSpectralConst.cpp \
          GModelSpectralComposite.cpp \
          GModelSpectralExpInvPlaw.cpp \
          GModelSpectralExpPlaw.cpp \
          GModelSpectralExponential.cpp \
          GModelSpectralSuperExpPlaw.cpp \
          GModelSpectralFunc.cpp \
          GModelSpectralGauss.cpp \
          GModelSpectralLogParabola.cpp \
          GModelSpectralMultiplicative.cpp \
          GModelSpectralBins.cpp \
          GModelSpectralNodes.cpp \
          GModelSpectralPlaw.cpp \
          GModelSpectralPlawPhotonFlux.cpp \
          GModelSpectralPlawEnergyFlux.cpp \
          GModelSpectralSmoothBrokenPlaw.cpp \
          GModelSpectralTable.cpp \
          GModelSpectralTablePar.cpp \
          GModelSpectralTablePars.cpp \
          GModelTemporal.cpp \
          GModelTemporalRegistry.cpp \
          GModelTemporalConst.cpp \
          GModelTemporalLightCurve.cpp \
          GModelTemporalPhaseCurve.cpp \
          GException_model.cpp

# Build libtool library
noinst_LTLIBRARIES = libmodel.la
libmodel_la_SOURCES = $(sources)
libmodel_la_LDFLAGS =
