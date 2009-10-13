/***************************************************************************
 *                     GMatrixTools.hpp  -  Matrix tools                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GMATRIXTOOLS_HPP
#define GMATRIXTOOLS_HPP

/* __ Includes ___________________________________________________________ */
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Prototypes ________________________________________________________ */
GMatrix convert(const GSymMatrix& m);
GMatrix convert(const GSparseMatrix& m);

#endif /* GMATRIXTOOLS_HPP */
