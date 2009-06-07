/***************************************************************************
 *     GSparseSymbolic.hpp  -  sparse matrix symbolic analysis class       *
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

#ifndef GSPARSESYMBOLIC_HPP
#define GSPARSESYMBOLIC_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GSparseMatrix.hpp"

/* __ Definitions ________________________________________________________ */

/* __ Enumerators ________________________________________________________ */

/* __ Structures _________________________________________________________ */

/* __ Prototypes _________________________________________________________ */


/***************************************************************************
 *                     GSparseSymbolic class definition                    *
 ***************************************************************************/
class GSparseSymbolic {
  // Friend classes
  friend class GSparseMatrix;
  friend class GSparseNumeric;

  // Binary operator friends

  // I/O friends
  friend std::ostream& operator<< (std::ostream& os, const GSparseSymbolic& s);

  // Friend functions

public:
  // Constructors and destructors
  GSparseSymbolic();
//  GSparseSymbolic(const GSparseSymbolic& s) { cout << "would sym to" << endl;}
  ~GSparseSymbolic();

  // Assignment operator
  GSparseSymbolic& operator= (const GSparseSymbolic& s);

  // Binary operators

  // Unary operators

  // Functions
  void cholesky_symbolic_analysis(int order, const GSparseMatrix& m);

private:
  // Functions
  int*        cs_amd (int order, const GSparseMatrix* A);
  int*        cs_counts(const GSparseMatrix* A, const int* parent, const int* post, int ata);
  int*        cs_etree(const GSparseMatrix* A, int ata);
  int         cs_fkeep(GSparseMatrix* A, int(*fkeep)(int, int, double, void*), void* other);
  int         cs_leaf(int i, int j, const int* first, int* maxfirst, int* prevleaf, int* ancestor, int* jleaf);
  int*        cs_pinv(int const* p, int n);
  int*        cs_post(const int* parent, int n);
  int         cs_tdfs(int j, int k, int* head, const int* next, int* post, int* stack);
  static void init_ata(const GSparseMatrix* AT, const int* post, int* wrk_int, int** head, int** next);
  static int  cs_diag(int i, int j, double aij, void* other);
  static int  cs_wclear(int mark, int lemax, int* w, int n);

  // Data
  int*   m_pinv;        // Inverse row permutation for QR, fill reduce permutation for Cholesky
  int*   m_q;           // Fill-reducing column permutation for LU and QR
  int*   m_parent;      // elimination tree for Cholesky and QR
  int*   m_cp;	        // column pointers for Cholesky, row counts for QR
  int*   m_leftmost;    // leftmost[i] = min(find(A(i,:))), for QR
  int    m_m2;	        // # of rows for QR, after adding fictitious rows
  double m_lnz;         // # entries in L for LU or Cholesky; in V for QR
  double m_unz;         // # entries in U for LU; in R for QR
  int    m_n_pinv;      // Number of elements in m_pinv
  int    m_n_q;         // Number of elements in m_q
  int    m_n_parent;    // Number of elements in m_parent
  int    m_n_cp;        // Number of elements in m_cp
  int    m_n_leftmost;  // Number of elements in m_leftmost
};


/***************************************************************************
 *                              Inline members                             *
 ***************************************************************************/


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/

#endif /* GSPARSESYMBOLIC_HPP */
