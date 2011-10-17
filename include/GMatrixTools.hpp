/***************************************************************************
 *                     GMatrixTools.hpp  -  Matrix tools                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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
 ***************************************************************************/

#ifndef GMATRIXTOOLS_HPP
#define GMATRIXTOOLS_HPP

/* __ Includes ___________________________________________________________ */
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Prototypes _________________________________________________________ */
GMatrix       matrix(const GSymMatrix& m);
GMatrix       matrix(const GSparseMatrix& m);
GSymMatrix    sym_matrix(const GMatrix& m);
GSymMatrix    sym_matrix(const GSparseMatrix& m);
GSparseMatrix sparse_matrix(const GMatrix& m);
GSparseMatrix sparse_matrix(const GSymMatrix& m);


#endif /* GMATRIXTOOLS_HPP */
