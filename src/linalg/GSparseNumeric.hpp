/***************************************************************************
 *      GSparseNumeric.hpp  -  sparse matrix numeric analysis class        *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2006 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GSPARSENUMERIC_HPP
#define GSPARSENUMERIC_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GSparseMatrix.hpp"
#include "GSparseSymbolic.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Definitions ________________________________________________________ */

/* __ Enumerators ________________________________________________________ */

/* __ Structures _________________________________________________________ */

/* __ Prototypes _________________________________________________________ */


/***************************************************************************
 *                      GSparseNumeric class definition                    *
 ***************************************************************************/
class GSparseNumeric {
  // Friend classes
  friend class GSparseMatrix;

  // I/O friends
  friend ostream& operator<< (ostream& os, const GSparseNumeric& n);

public:
  // Constructors and destructors
  GSparseNumeric();
  ~GSparseNumeric();

  // Assignment operator
  GSparseNumeric& operator= (const GSparseNumeric& n);

  // Functions
  void cholesky_numeric_analysis(const GSparseMatrix& m, const GSparseSymbolic& s);

private:
  // Functions
  int cs_ereach(const GSparseMatrix* A, int k, const int* parent, int* s, int* w);

  // Data
  GSparseMatrix* m_L;        // L for LU and Cholesky, V for QR
  GSparseMatrix* m_U;        // U for LU, R for QR, not used for Cholesky
  int*           m_pinv;     // partial pivoting for LU
  double*        m_B;        // beta [0..n-1] for QR
  int            m_n_pinv;   // Number of elements in m_pinv
  int            m_n_B;      // Number of elements in m_B
};


/***************************************************************************
 *                              Inline members                             *
 ***************************************************************************/


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/

#endif /* GSPARSENUMERIC_HPP */
