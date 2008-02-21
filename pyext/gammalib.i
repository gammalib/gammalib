/***************************************************************************
 *                         gammalib - SWIG file                            *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall -includeall gammalib.i                          *
 ***************************************************************************/
%module gammalib

/* Linear Algebra */
#include "GVector.i"
#include "GMatrixBase.i"
#include "GMatrix.i"
#include "GSparseMatrix.i"
#include "GSymMatrix.i"

/* Model fitting */
#include "GFitPar.i"
