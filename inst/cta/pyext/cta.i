/***************************************************************************
 *        cta - Cherenkov Telescope Array support python bindings          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
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
%module cta
%feature("autodoc", "1");

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="support") "GException.i";

/* __ Inform about base classes __________________________________________ */
%import(module="obs") "GObservation.i";
%import(module="obs") "GEvent.i";
%import(module="obs") "GEventBin.i";
%import(module="obs") "GEventAtom.i";
%import(module="obs") "GEvents.i";
%import(module="obs") "GEventCube.i";
%import(module="obs") "GEventList.i";
%import(module="obs") "GPointing.i";
%import(module="obs") "GResponse.i";
%import(module="obs") "GInstDir.i";
%import(module="obs") "GRoi.i";
%import(module="model") "GModel.i";
%import(module="model") "GModelData.i";

/* __ CTA ________________________________________________________________ */
%include "GCTAObservation.i"
%include "GCTAEventCube.i"
%include "GCTAEventList.i"
%include "GCTAEventBin.i"
%include "GCTAEventAtom.i"
%include "GCTAPointing.i"
%include "GCTAResponse.i"
%include "GCTAInstDir.i"
%include "GCTARoi.i"
%include "GCTAModelRadial.i"
%include "GCTAModelRadialRegistry.i"
%include "GCTAModelRadialGauss.i"
%include "GCTAModelRadialPolynom.i"
%include "GCTAModelRadialProfile.i"
%include "GCTAModelRadialAcceptance.i"
%include "GCTADir.i"



