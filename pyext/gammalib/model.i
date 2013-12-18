/***************************************************************************
 *                          model.i - Model module                         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall model.i                                         *
 ***************************************************************************/
/**
 * @file model.i
 * @brief Model module
 * @author Juergen Knoedlseder
 */
%module model
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
#include "GException.hpp"
#include "GTools.hpp"
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";
%import(module="gammalib.base") "GContainer.i";
%import(module="gammalib.base") "GRegistry.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i";

/* __ Inform about base classes __________________________________________ */
%import(module="gammalib.opt") "GOptimizerPar.i";
%import(module="gammalib.opt") "GOptimizerPars.i";

/* __ Model handling _____________________________________________________ */
%include "GModelPar.i"
%include "GModels.i"
%include "GModel.i"
%include "GModelRegistry.i"
%include "GModelSky.i"
%include "GModelData.i"
%include "GModelSpatial.i"
%include "GModelSpatialRegistry.i"
%include "GModelSpatialPointSource.i"
%include "GModelSpatialRadial.i"
%include "GModelSpatialRadialDisk.i"
%include "GModelSpatialRadialGauss.i"
%include "GModelSpatialRadialShell.i"
%include "GModelSpatialElliptical.i"
%include "GModelSpatialEllipticalDisk.i"
%include "GModelSpatialDiffuse.i"
%include "GModelSpatialDiffuseConst.i"
%include "GModelSpatialDiffuseCube.i"
%include "GModelSpatialDiffuseMap.i"
%include "GModelSpectral.i"
%include "GModelSpectralRegistry.i"
%include "GModelSpectralConst.i"
%include "GModelSpectralExpPlaw.i"
%include "GModelSpectralFunc.i"
%include "GModelSpectralNodes.i"
%include "GModelSpectralPlaw.i"
%include "GModelSpectralPlaw2.i"
%include "GModelSpectralLogParabola.i"
%include "GModelTemporal.i"
%include "GModelTemporalRegistry.i"
%include "GModelTemporalConst.i"
