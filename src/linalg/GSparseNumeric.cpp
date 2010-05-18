/***************************************************************************
 *      GSparseNumeric.cpp  -  sparse matrix numeric analysis class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSparseNumeric.hpp"

/* __ Macros _____________________________________________________________ */
#define CS_FLIP(i) (-(i)-2)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_MARKED(w,j) (w [j] < 0)


/*==========================================================================
 =                                                                         =
 =                   GSparseNumeric constructors/destructors               =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                        GSparseNumeric constructor                       *
 ***************************************************************************/
GSparseNumeric::GSparseNumeric()
{
  // Initialise private members for clean destruction
  m_L      = NULL;
  m_U      = NULL;
  m_pinv   = NULL;
  m_B      = NULL;
  m_n_pinv = 0;
  m_n_B    = 0;
  
  // Return
  return;
}


/***************************************************************************
 *                        GSparseNumeric destructor                       *
 ***************************************************************************/
GSparseNumeric::~GSparseNumeric()
{
  // De-allocate only if memory has indeed been allocated
  if (m_L    != NULL) delete m_L;
  if (m_U    != NULL) delete m_U;
  if (m_pinv != NULL) delete [] m_pinv;
  if (m_B    != NULL) delete [] m_B;
  
  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                         GSparseNumeric operators                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                     GSparseNumeric assignment operator                  *
 ***************************************************************************/
GSparseNumeric& GSparseNumeric::operator= (const GSparseNumeric& n)
{ 
    // Execute only if object is not identical
    if (this != &n) {

        // De-allocate only if memory has indeed been allocated
        if (m_L    != NULL) delete m_L;
        if (m_U    != NULL) delete m_U;
        if (m_pinv != NULL) delete [] m_pinv;
        if (m_B    != NULL) delete [] m_B;

        // Initialise private members for clean destruction
        m_L      = NULL;
        m_U      = NULL;
        m_pinv   = NULL;
        m_B      = NULL;
        m_n_pinv = 0;
        m_n_B    = 0;

	    // Copy m_L if it exists
	    if (n.m_L != NULL)
	        m_L = new GSparseMatrix(*n.m_L);

	    // Copy m_U if it exists
	    if (n.m_U != NULL)
	        m_U = new GSparseMatrix(*n.m_U);
	
	    // Copy m_pinv array if it exists
	    if (n.m_pinv != NULL && n.m_n_pinv > 0) {
	        m_pinv = new int[n.m_n_pinv];
            for (int i = 0; i < n.m_n_pinv; ++i)
                m_pinv[i] = n.m_pinv[i];
	        m_n_pinv = n.m_n_pinv;
	    }

	    // Copy m_B array if it exists
	    if (n.m_B != NULL && n.m_n_B > 0) {
	        m_B = new double[n.m_n_B];
            for (int i = 0; i < n.m_n_B; ++i)
                m_B[i] = n.m_B[i];
	        m_n_B = n.m_n_B;
	    }

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                      GSparseNumeric member functions                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                        cholesky_numeric_analysis                        *
 * ----------------------------------------------------------------------- *
 * L = chol (A, [pinv parent cp]), pinv is optional. This function         *
 * allocates the memory to hold the information.                           *
 * ----------------------------------------------------------------------- *
 * Input:   A                    Sparse matrix                             *
 *          S                    Symbolic analysis of sparse matrix        *
 ***************************************************************************/
void GSparseNumeric::cholesky_numeric_analysis(const GSparseMatrix& A, 
                                               const GSparseSymbolic& S)
{
  // De-allocate memory that has indeed been previously allocated
  if (m_L    != NULL) delete m_L;
  if (m_U    != NULL) delete m_U;
  if (m_pinv != NULL) delete [] m_pinv;
  if (m_B    != NULL) delete [] m_B;

  // Initialise members
  m_L      = NULL;
  m_U      = NULL;
  m_pinv   = NULL;
  m_B      = NULL;
  m_n_pinv = 0;
  m_n_B    = 0;

  // Return if arrays in the symbolic analysis have not been allocated
  if (!S.m_cp || !S.m_parent) return;
  
  // Declare
  double lki;
  int top, i, p, k;

  // Assign input matrix attributes
  int n = A.m_cols;

    // Allocate int workspace
    int  wrk_size = 2*n;
    int* wrk_int  = new int[wrk_size];

    // Allocate double workspace
    wrk_size = n;
    double* wrk_double = new double[wrk_size];

  // Assign pointers
  int* cp     = S.m_cp;
  int* pinv   = S.m_pinv;
  int* parent = S.m_parent;

  // Assign C = A(p,p) where A and C are symmetric and the upper part stored
  GSparseMatrix C = (pinv) ? cs_symperm(A, pinv) : (A);

  // Assign workspace pointer
  int*    c = wrk_int;
  int*    s = wrk_int + n;
  double* x = wrk_double;
  
  // Assign C matrix pointers 
  int*    Cp = C.m_colstart;
  int*    Ci = C.m_rowinx; 
  double* Cx = C.m_data;
  
    // Allocate L matrix
    m_L = new GSparseMatrix(n, n, cp[n]);

  // Assign L matrix pointers 
  int*    Lp = m_L->m_colstart;
  int*    Li = m_L->m_rowinx; 
  double* Lx = m_L->m_data;

  // Initialise column pointers of L and c
  for (k = 0; k < n; ++k) 
    Lp[k] = c[k] = cp[k];
	
  // Compute L(:,k) for L*L' = C
  for (k = 0; k < n; k++) {
  
	// Nonzero pattern of L(k,:). 
	// Returns -1 if parent = NULL, s = NULL or c = NULL
    top  = cs_ereach(&C, k, parent, s, c);      // find pattern of L(k,:)
	x[k] = 0;                                   // x (0:k) is now zero
	
	// x = full(triu(C(:,k)))
	for (p = Cp[k]; p < Cp[k+1]; p++) {
      if (Ci[p] <= k) 
	    x[Ci[p]] = Cx[p];
	}
	
	// d = C(k,k)
	double d = x[k];
	
	// Clear x for k+1st iteration
	x[k] = 0;
	
	// Triangular solve: Solve L(0:k-1,0:k-1) * x = C(:,k)
	for ( ; top < n; top++) {
      i    = s[top];                            // s [top..n-1] is pattern of L(k,:)
      lki  = x[i]/Lx[Lp[i]];                    // L(k,i) = x (i) / L(i,i)
      x[i] = 0;                                 // clear x for k+1st iteration
	  for (p = Lp[i]+1; p < c[i]; p++)
		x[Li[p]] -= Lx[p] * lki;
	  d    -= lki * lki;                        // d = d - L(k,i)*L(k,i)
	  p     = c[i]++;
	  Li[p] = k;                                // store L(k,i) in column i
	  Lx[p] = lki;
	}
	
	// Compute L(k,k)
	if (d <= 0)
	  throw GException::matrix_not_pos_definite(
	        "GSparseNumeric::cholesky_numeric_analysis(GSparseMatrix&, const GSparseSymbolic&)",
			k, d);

    // Store L(k,k) = sqrt (d) in column k
	p     = c[k]++;
	Li[p] = k;
	Lx[p] = sqrt(d);
  }
  
  // Finalize L
  Lp[n] = cp[n];
  
  // Free workspace
  delete [] wrk_int;
  delete [] wrk_double;

  // Return void
  return; 
}


/*==========================================================================
 =                                                                         =
 =                      GSparseNumeric private functions                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                 cs_ereach                               *
 * ----------------------------------------------------------------------- *
 * Find nonzero pattern of Cholesky L(k,1:k-1) using etree and             *
 * triu(A(:,k))                                                            *
 * ----------------------------------------------------------------------- *
 * Input:   A                    Sparse matrix                             *
 *          k                    Node                                      *
 *          parent               ...                                       *
 *          s                    ...                                       *
 *          w                    ...                                       *
 * Input:   top                  Index, where s[top..n-1] contains pattern *
 *                               of L(k,:); -1 if arrays have not been     *
 *                               allocated                                 *
 ***************************************************************************/
int GSparseNumeric::cs_ereach(const GSparseMatrix* A, int k, 
                                          const int* parent, int* s, int* w)
{
  // Return -1 if arrays have not been allocated
  if (!parent || !s || !w) return (-1);

  // Assign A matrix attributes and pointers 
  int  n   = A->m_cols;
  int  top = n;
  int* Ap  = A->m_colstart;
  int* Ai  = A->m_rowinx; 
    
  // Mark node k as visited
  CS_MARK(w, k);
  
  // Loop over elements of node
  int i;
  for (int p = Ap[k]; p < Ap[k+1]; ++p) {
	i = Ai[p];                  // A(i,k) is nonzero
	if (i > k) continue;        // only use upper triangular part of A

    // Traverse up etree
	int len;
	for (len = 0; !CS_MARKED(w,i); i = parent[i]) {
      s[len++] = i;             // L(k,i) is nonzero
	  CS_MARK(w, i);            // mark i as visited
	}
	
	// Push path onto stack
	while (len > 0) 
	  s[--top] = s[--len];
  }
  
  // Unmark all nodes
  for (int p = top; p < n; ++p) 
    CS_MARK(w, s[p]);
	
  // Unmark node k
  CS_MARK(w, k);
  
  // Return index top, where s[top..n-1] contains pattern of L(k,:)
  return top;
}


/*==========================================================================
 =                                                                         =
 =                   GSparseNumeric static member functions                =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                           GSparseNumeric friends                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSparseNumeric& n)
{
  // Put header is stream
  os << "=== GSparseNumeric ===";

  // Show L or V matrix
  if (n.m_L != NULL) {
    os << std::endl <<
	   " === L or V ";
	os << *n.m_L << std::endl;
  }

  // Show U or R matrix
  if (n.m_U != NULL) {
    os << std::endl <<
       " === U or R matrix ";
	os << *n.m_U << std::endl;
  }

  // Show inverse permutation if it exists
  if (n.m_pinv != NULL) {
    os << std::endl << " Inverse permutation (of " << n.m_n_pinv << " elements)" << std::endl;
	for (int i = 0; i < n.m_n_pinv; ++i)
	  os << " " << n.m_pinv[i];
  } 

  // Show beta array if it exists
  if (n.m_B != NULL) {
    os << std::endl << " Beta[0.." << n.m_n_B << "]" << std::endl;
	for (int i = 0; i < n.m_n_B; ++i)
	  os << " " << n.m_B[i];
  } 
  
  // Return output stream
  return os;
}

/*==========================================================================
 =                                                                         =
 =                      GSparseNumeric exception classes                   =
 =                                                                         =
 ==========================================================================*/
