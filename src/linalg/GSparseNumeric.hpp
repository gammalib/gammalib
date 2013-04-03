/***************************************************************************
 *       GSparseNumeric.hpp - Sparse matrix numeric analysis class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2013 by Juergen Knoedlseder                         *
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
/**
 * @file GSparseNumeric.hpp
 * @brief Sparse matrix numeric analysis class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPARSENUMERIC_HPP
#define GSPARSENUMERIC_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GMatrixSparse.hpp"
#include "GSparseSymbolic.hpp"

/* __ Definitions ________________________________________________________ */

/* __ Enumerators ________________________________________________________ */

/* __ Structures _________________________________________________________ */

/* __ Prototypes _________________________________________________________ */


/***********************************************************************//**
 * @class GSparseNumeric
 *
 * @brief Sparse matrix numeric analysis class
 ***************************************************************************/
class GSparseNumeric {

  // Friend classes
  friend class GMatrixSparse;

  // I/O friends
  friend std::ostream& operator<< (std::ostream& os, const GSparseNumeric& n);

public:
  // Constructors and destructors
  GSparseNumeric(void);
  virtual ~GSparseNumeric(void);

  // Assignment operator
  GSparseNumeric& operator=(const GSparseNumeric& n);

  // Functions
  void cholesky_numeric_analysis(const GMatrixSparse& m, const GSparseSymbolic& s);

private:
  // Functions
  int cs_ereach(const GMatrixSparse* A, int k, const int* parent, int* s, int* w);

  // Data
  GMatrixSparse* m_L;        // L for LU and Cholesky, V for QR
  GMatrixSparse* m_U;        // U for LU, R for QR, not used for Cholesky
  int*           m_pinv;     // partial pivoting for LU
  double*        m_B;        // beta [0..n-1] for QR
  int            m_n_pinv;   // Number of elements in m_pinv
  int            m_n_B;      // Number of elements in m_B
};

#endif /* GSPARSENUMERIC_HPP */
