/***************************************************************************
 *                cta.i - Cherenkov Telescope Array module                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
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
 * swig -c++ -python -Wall cta.i                                           *
 ***************************************************************************/
/**
 * @file cat.i
 * @brief Cherenkov Telescope Array module
 * @author Juergen Knoedlseder
 */
%module cta
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
%import(module="gammalib.obs") "GObservation.i";
%import(module="gammalib.obs") "GEvent.i";
%import(module="gammalib.obs") "GEventBin.i";
%import(module="gammalib.obs") "GEventAtom.i";
%import(module="gammalib.obs") "GEvents.i";
%import(module="gammalib.obs") "GEventCube.i";
%import(module="gammalib.obs") "GEventList.i";
%import(module="gammalib.obs") "GResponse.i";
%import(module="gammalib.obs") "GInstDir.i";
%import(module="gammalib.obs") "GRoi.i";
%import(module="gammalib.model") "GModel.i";
%import(module="gammalib.model") "GModelData.i";

/* __ CTA ________________________________________________________________ */
%include "GCTAObservation.i"
%include "GCTAOnOffObservation.i"
%include "GCTAOnOffObservations.i"
%include "GCTAEventCube.i"
%include "GCTAEventList.i"
%include "GCTAEventBin.i"
%include "GCTAEventAtom.i"
%include "GCTAPointing.i"
%include "GCTAResponse.i"
%include "GCTAResponseTable.i"
%include "GCTAAeff.i"
%include "GCTAAeffPerfTable.i"
%include "GCTAAeffArf.i"
%include "GCTAAeff2D.i"
%include "GCTAPsf.i"
%include "GCTAPsfPerfTable.i"
%include "GCTAPsfVector.i"
%include "GCTAPsf2D.i"
%include "GCTAPsfKing.i"
%include "GCTAEdisp.i"
%include "GCTAEdispPerfTable.i"
%include "GCTABackground.i"
%include "GCTABackground3D.i"
%include "GCTAInstDir.i"
%include "GCTARoi.i"
%include "GCTAModelBackground.i"
%include "GCTAModelInstBackground.i"
%include "GCTAModelRadial.i"
%include "GCTAModelRadialRegistry.i"
%include "GCTAModelRadialGauss.i"
%include "GCTAModelRadialPolynom.i"
%include "GCTAModelRadialProfile.i"
%include "GCTAModelRadialAcceptance.i"



