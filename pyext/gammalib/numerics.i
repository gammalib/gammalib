/***************************************************************************
 *                      numerics.i - Numerics module                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2020 by Juergen Knoedlseder                         *
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
 * swig -c++ -python -Wall numerics.i                                      *
 ***************************************************************************/
/**
 * @file numerics.i
 * @brief Numerics module
 * @author Juergen Knoedlseder
 */
%module numerics
%feature("autodoc", "1");

/* __ Headers needed for compilation _____________________________________ */
%{
#include <stddef.h>
#include "GException.hpp"
#include "GTools.hpp"
%}

/* __ Include standard typemaps for vectors and strings __________________ */
%include stl.i
%include std_vector.i
%include std_complex.i
%template(vectori) std::vector<int>;
%template(vectorc) std::vector<std::complex<double> >;

/* __ Include GammaLib typemaps __________________________________________ */
%include typemap_GChatter.i
%include typemap_GTuple.i
%include typemap_GFilename.i

/* __ Include interface classes __________________________________________ */
%import(module="gammalib.base") "GBase.i";

/* __ Make sure that exceptions are catched ______________________________ */
%import(module="gammalib.support") "GException.i"; 

/* __ Numerics module ____________________________________________________ */
%include "GDerivative.i"
%include "GFunction.i"
%include "GFunctions.i"
%include "GIntegral.i"
%include "GMath.i"
%include "GNdarray.i"
%include "GFft.i"
%include "GFftWavetable.i"
