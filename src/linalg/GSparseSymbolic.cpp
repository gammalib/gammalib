/***************************************************************************
 *     GSparseSymbolic.cpp  -  sparse matrix symbolic analysis class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2012 by Juergen Knoedlseder                         *
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
 * @file GSparseSymbolic.cpp
 * @brief Sparse matrix symbolic analysis class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GSparseMatrix.hpp"
#include "GSparseSymbolic.hpp"

/* __ Macros _____________________________________________________________ */
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)

/* __ Definitions ________________________________________________________ */
//#define G_DEBUG_SPARSE_SYM_AMD_ORDERING               // Debug AMD ordering
//#define G_DEBUG_SPARSE_CHOLESKY           // Analyse Cholesky decomposition


/*==========================================================================
 =                                                                         =
 =                  GSparseSymbolic constructors/destructors               =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                       GSparseSymbolic constructor                       *
 ***************************************************************************/
GSparseSymbolic::GSparseSymbolic()
{
  // Initialise private members for clean destruction
  m_n_pinv     = 0;
  m_n_q        = 0;
  m_n_parent   = 0;
  m_n_cp       = 0;
  m_n_leftmost = 0;
  m_pinv       = NULL;
  m_q          = NULL;
  m_parent     = NULL;
  m_cp         = NULL;
  m_leftmost   = NULL;
  m_m2         = 0;
  m_lnz        = 0.0;
  m_unz        = 0.0;

  // Return
  return;
}


/***************************************************************************
 *                        GSparseSymbolic destructor                       *
 ***************************************************************************/
GSparseSymbolic::~GSparseSymbolic()
{
  // De-allocate only if memory has indeed been allocated
  if (m_pinv     != NULL) delete [] m_pinv;
  if (m_q        != NULL) delete [] m_q;
  if (m_parent   != NULL) delete [] m_parent;
  if (m_cp       != NULL) delete [] m_cp;
  if (m_leftmost != NULL) delete [] m_leftmost;

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                        GSparseSymbolic operators                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                    GSparseSymbolic assignment operator                  *
 ***************************************************************************/
GSparseSymbolic& GSparseSymbolic::operator= (const GSparseSymbolic& s)
{ 
    // Execute only if object is not identical
    if (this != &s) {

      // De-allocate only if memory has indeed been allocated
      if (m_pinv     != NULL) delete [] m_pinv;
      if (m_q        != NULL) delete [] m_q;
      if (m_parent   != NULL) delete [] m_parent;
      if (m_cp       != NULL) delete [] m_cp;
      if (m_leftmost != NULL) delete [] m_leftmost;

      // Initialise private members for clean destruction
      m_n_pinv     = 0;
      m_n_q        = 0;
      m_n_parent   = 0;
      m_n_cp       = 0;
      m_n_leftmost = 0;
      m_pinv       = NULL;
      m_q          = NULL;
      m_parent     = NULL;
      m_cp         = NULL;
      m_leftmost   = NULL;
      m_m2         = 0;
      m_lnz        = 0.0;
      m_unz        = 0.0;

      // Copy data members
      m_m2  = s.m_m2;
      m_lnz = s.m_lnz;
      m_unz = s.m_unz;
	
	  // Copy m_pinv array if it exists
	  if (s.m_pinv != NULL && s.m_n_pinv > 0) {
	      m_pinv = new int[s.m_n_pinv];
          for (int i = 0; i < s.m_n_pinv; ++i)
              m_pinv[i] = s.m_pinv[i];
	      m_n_pinv = s.m_n_pinv;
	  }

	  // Copy m_q array if it exists
	  if (s.m_q != NULL && s.m_n_q > 0) {
	      m_q = new int[s.m_n_q];
          for (int i = 0; i < s.m_n_q; ++i)
              m_q[i] = s.m_q[i];
	      m_n_q = s.m_n_q;
      }

      // Copy m_parent array if it exists
	  if (s.m_parent != NULL && s.m_n_parent > 0) {
	      m_parent = new int[s.m_n_parent];
          for (int i = 0; i < s.m_n_parent; ++i)
              m_parent[i] = s.m_parent[i];
	      m_n_parent = s.m_n_parent;
	  }

	  // Copy m_cp array if it exists
	  if (s.m_cp != NULL && s.m_n_cp > 0) {
	      m_cp = new int[s.m_n_cp];
          for (int i = 0; i < s.m_n_cp; ++i)
              m_cp[i] = s.m_cp[i];
	      m_n_cp = s.m_n_cp;
	  }

	  // Copy m_leftmost array if it exists
	  if (s.m_leftmost != NULL && s.m_n_leftmost > 0) {
	      m_leftmost = new int[s.m_n_leftmost];
          for (int i = 0; i < s.m_n_leftmost; ++i)
              m_leftmost[i] = s.m_leftmost[i];
	      m_n_leftmost = s.m_n_leftmost;
	  }

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                     GSparseSymbolic member functions                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                       cholesky_symbolic_analysis                        *
 * ----------------------------------------------------------------------- *
 * Ordering and symbolic analysis for a Cholesky factorization. This       *
 * function allocates the memory to hold the information.                  *
 * ----------------------------------------------------------------------- *
 * Input:   order                0: natural                                *
 *                               1: Cholesky decomposition                 *
 *          m                    Sparse matrix                             *
 ***************************************************************************/
void GSparseSymbolic::cholesky_symbolic_analysis(int order, 
												 const GSparseMatrix& m)
{
  // Debug
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << endl << "GSparseSymbolic::cholesky_symbolic_analysis entered" << endl;
  cout << " order ......................: " << order << endl;
  cout << " number of rows .............: " << m.m_rows << endl;
  cout << " number of columns ..........: " << m.m_cols << endl;
  cout << " number of non-zero elements : " << m.m_elements << endl;
  #endif

  // De-allocate memory that has indeed been previously allocated
  if (m_pinv     != NULL) delete [] m_pinv;
  if (m_q        != NULL) delete [] m_q;
  if (m_parent   != NULL) delete [] m_parent;
  if (m_cp       != NULL) delete [] m_cp;
  if (m_leftmost != NULL) delete [] m_leftmost;
  
  // Initialise members
  m_n_pinv     = 0;
  m_n_q        = 0;
  m_n_parent   = 0;
  m_n_cp       = 0;
  m_n_leftmost = 0;
  m_pinv       = NULL;
  m_q          = NULL;
  m_parent     = NULL;
  m_cp         = NULL;
  m_leftmost   = NULL;
  m_m2         = 0;
  m_lnz        = 0.0;
  m_unz        = 0.0;

  // Check if order type is valid
  if (order < 0 || order > 1)
    throw GException::invalid_order(
		  "GSparseSymbolic::cholesky_symbolic_analysis(int, GSparseMatrix*)",
	      order, 0, 1);

  // Assign input matrix attributes
  int n = m.m_cols;

  // P = amd(A+A'), or natural
  int* P = cs_amd(order, &m);
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << " AMD ordering permutation (" << P << ")" << endl;
  if (P != NULL) {
    for (int i = 0; i < n; ++i)
      cout << " " << P[i];
    cout << endl;
  }
  #endif
  
  // Find inverse permutation and store it in class member 'm_pinv'.
  // Note that if P = NULL or n < 1 this function returns NULL.
  m_pinv   = cs_pinv(P, n);
  m_n_pinv = n;
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << " Inverse permutation (" << m_pinv << ")" << endl;
  if (m_pinv != NULL) {
    for (int i = 0; i < n; ++i)
      cout << " " << m_pinv[i];
    cout << endl;
  }
  #endif
  
  // Delete workspace
  if (P != NULL) delete [] P;
  
  // C = spones(triu(A(P,P)))
  GSparseMatrix C = cs_symperm(m, m_pinv);
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << " C = spones(triu(A(P,P))) " << C << endl;
  #endif
  
  // Find etree of C and store it in class member 'm_parent'
  m_parent   = cs_etree(&C, 0);
  m_n_parent = n;
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << " Elimination tree (" << m_parent << ")" << endl;
  if (m_parent != NULL) {
    for (int i = 0; i < n; ++i)
      cout << " " << m_parent[i];
    cout << endl;
  }
  #endif

  // Post order the etree
  // Note that if m_parent = NULL or n < 1 this function returns NULL.
  int* post = cs_post(m_parent, n);
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << " Post ordered elimination tree (" << post << ")" << endl;
  if (post != NULL) {
    for (int i = 0; i < n; ++i)
      cout << " " << post[i];
    cout << endl;
  }
  #endif
  
  // Find column counts of chol(C)
  // Note that if m_parent = NULL or post = NULL this function returns NULL.
  int* c = cs_counts(&C, m_parent, post, 0);
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << " Column counts (" << c << ")" << endl;
  if (c != NULL) {
    for (int i = 0; i < n; ++i)
      cout << " " << c[i];
    cout << endl;
  }
  #endif

  // Delete workspace
  if (post != NULL) delete [] post;
  
    // Allocate column pointers for Cholesky decomposition
    m_n_cp = n+1;
    m_cp   = new int[m_n_cp];

  // Find column pointers for L
  m_unz = m_lnz = cs_cumsum(m_cp, c, n);
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << " Column pointers for L (" << m_cp << ")" << endl;
  if (m_cp != NULL) {
    for (int i = 0; i < m_n_cp; ++i)
      cout << " " << m_cp[i];
    cout << endl;
  }
  cout << " Number of non-zero elements in L: " << m_lnz << endl;
  #endif

  // Delete workspace
  if (c != NULL) delete [] c;
  
  // Delete symbolic analysis if it is not valid
  if (m_lnz < 0) {

    // De-allocate memory that has indeed been previously allocated
    if (m_pinv     != NULL) delete [] m_pinv;
    if (m_q        != NULL) delete [] m_q;
    if (m_parent   != NULL) delete [] m_parent;
    if (m_cp       != NULL) delete [] m_cp;
    if (m_leftmost != NULL) delete [] m_leftmost;
  
    // Initialise members
    m_n_pinv     = 0;
    m_n_q        = 0;
    m_n_parent   = 0;
    m_n_cp       = 0;
    m_n_leftmost = 0;
    m_pinv       = NULL;
    m_q          = NULL;
    m_parent     = NULL;
    m_cp         = NULL;
    m_leftmost   = NULL;
    m_m2         = 0;
    m_lnz        = 0.0;
    m_unz        = 0.0;
  }

  // Debug
  #if defined(G_DEBUG_SPARSE_CHOLESKY)
  cout << "GSparseSymbolic::cholesky_symbolic_analysis finished" << endl;
  #endif
  
  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                     GSparseSymbolic private functions                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                 cs_amd                                  *
 * ----------------------------------------------------------------------- *
 * Applies an approximate minimum degree ordering algorithm to derive the  *
 * permutations that minimise the fill-in.                                 *
 * p = amd(A+A') if symmetric is true, or amd(A'A) otherwise               *
 * ----------------------------------------------------------------------- *
 * Input:   order                0: natural                                *
 *                               1: Cholesky decomposition                 *
 *                               2: LU decomposition                       *
 *                               3: QR decomposition                       *
 *          A                    Sparse matrix                             *
 * Output:  p[0..n]              integer array of n+1 elements             *
 ***************************************************************************/
int* GSparseSymbolic::cs_amd(int order, const GSparseMatrix* A)
{
  // Throw exception if order is not a decomposition
  if (order < 0 || order > 3)
    throw GException::invalid_order("GSparseSymbolic::cs_amd(int, const GSparseMatrix*)",
	                    order, 0, 3);
  
  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << endl << "GSparseSymbolic::cs_amd entered" << endl;
  cout << " order ......................: " << order << endl;
  cout << " number of rows .............: " << A->m_rows << endl;
  cout << " number of columns ..........: " << A->m_cols << endl;
  cout << " number of non-zero elements : " << A->m_elements << endl;
  #endif

  // Declare
  int dense;
  int  d, dk, dext, lemax = 0, e, elenk, eln, i, j, k, k1,
	k2, k3, jlast, ln, mindeg = 0, nvi, nvj, nvk, mark, wnvi,
	ok, nel = 0, p, p1, p2, p3, p4, pj, pk, pk1, pk2, pn, q;
  unsigned int h ;

  // Assign input matrix attributes
  int m = A->m_rows;
  int n = A->m_cols;


  // ======================================================================  
  // Step 1: Construct matrix C
  // ======================================================================  

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Step 1: Construct matrix C" << endl;
  #endif

  // Initialise matrix C. We do not care at this point about the actual
  // dimension since we attribute the matrix later anyways. It's just
  // that there is no empty constructor of a matrix, so we have to put
  // something by default ...
  GSparseMatrix C(m,n);

  // Get (logical) transpose of A: AT = A'
// NOTE: WE ONLY NEED THE LOGICAL TRANSPOSE HERE. HOWEVER, WE HAVE NO
// LOGICAL ADDITION SO FAR ...
//  GSparseMatrix AT = cs_transpose(*A, 0);
  GSparseMatrix AT = cs_transpose(*A, 1);

  // Find dense threshold
  dense = (int)CS_MAX(16, 10 * sqrt((double)n));
  dense = CS_MIN(n-2, dense);
  
  // Case A: Cholesky decomposition of a symmetric matrix: C = A + A'
  if (order == 1 && n == m)
//  NOTE: HERE WE ONLY NEED LOGICAL ADDITION
//	C = cs_add (A, AT, 0, 0);
	C = *A + AT;
  
  // Case B: LU decomposition: C = A' * A with no dense rows
  else if (order == 2) {
  
    // Assign A' matrix attributes
    int* ATp = AT.m_colstart;
    int* ATi = AT.m_rowinx; 
    
	// Drop dense columns from AT
	for (p2 = 0, j = 0; j < m; j++) {
	  p      = ATp[j];                       // column j of AT starts here
	  ATp[j] = p2;                           // new column j starts here
	  if (ATp[j+1] - p > dense) continue;    // skip dense col j
	  for ( ; p < ATp[j+1]; p++) 
	    ATi[p2++] = ATi[p] ;
	}
	
	// Finalise AT
	ATp[m] = p2;
	
	// Get (logical) transpose of AT: A2 = AT'
// NOTE: WE ONLY NEED THE LOGICAL TRANSPOSE HERE. HOWEVER, WE HAVE NO
// LOGICAL MULTIPLICATION SO FAR ...
//    GSparseMatrix A2 = cs_transpose(AT, 0);
    GSparseMatrix A2 = cs_transpose(AT, 1);
	
	// C = A' * A with no dense rows
//  NOTE: cs_multiply NOT YET IMPLEMENTED
//	C = A2 ? cs_multiply(AT, A2) : NULL;
	C = AT * A2;
  }
  
  // Case C: QR decomposition: C = A' * A
  else
//  NOTE: cs_multiply NOT YET IMPLEMENTED
//	C = cs_multiply(AT, A);
	C = AT * *A;

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Matrix C " << C << endl;
  #endif
  
  // Drop diagonal entries from C
  cs_fkeep(&C, &cs_diag, NULL);

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Dropped diagonal entries from matrix C " << C << endl;
  #endif
  
  // Assign C matrix attributes
  int* Cp  = C.m_colstart;
  int  cnz = Cp[n];

    // Allocate result array
    int* P = new int[n+1];

    // Allocate workspace
    int  wrk_size = 8*(n+1);
    int* wrk_int  = new int[wrk_size];

  // Add elbow room to C  
  int elbow_room = cnz/5 + 2*n;           // Request additional # of elements
  C.alloc_elements(cnz, elbow_room);      // Appand elements
  int nzmax = C.m_elements;               // Save the maximum # of entries for garbage coll.
  C.m_elements -= elbow_room;             // Reverse element increase of alloc_elements

  // (Re-)assign C matrix attributes
  Cp      = C.m_colstart;
  int* Ci = C.m_rowinx;

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Added elbow room to matrix C " << C << endl;
  #endif

  // Set workspace pointers
  int* len    = wrk_int; 
  int* nv     = wrk_int +   (n+1); 
  int* next   = wrk_int + 2*(n+1);
  int* head   = wrk_int + 3*(n+1); 
  int* elen   = wrk_int + 4*(n+1); 
  int* degree = wrk_int + 5*(n+1);
  int* w      = wrk_int + 6*(n+1); 
  int* hhead  = wrk_int + 7*(n+1);
  int* last   = P;                        // use P as workspace for last


  // ======================================================================  
  // Step 2: Initialize quotient graph
  // ======================================================================  

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Step 2: Initialize quotient graph" << endl;
  #endif

  // Setup array that contains for each column the number of non-zero elements
  for (k = 0 ; k < n ; k++) 
    len[k] = Cp[k+1] - Cp[k];
  len[n] = 0;

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Non-zero elements per column: ";
  for (k = 0 ; k <= n ; k++)
    cout << " " << len[k];
  cout << endl;
  #endif
	
  // Loop over all nodes (i.e. matrix columns)
  for (i = 0 ; i <= n ; i++) {
	head[i]   = -1;                       // degree list i is empty
	last[i]   = -1;
	next[i]   = -1;
	hhead[i]  = -1;                       // hash list i is empty
	nv[i]     = 1;                        // node i is just one node
	w[i]      = 1;                        // node i is alive
	elen[i]   = 0;                        // Ek of node i is empty
	degree[i] = len[i];                   // degree of node i
  }
  
  // Clear w
  mark = cs_wclear(0, 0, w, n);
  
  // End-point initialisations
  elen[n] = -2;                           // n is a dead element
  Cp[n]   = -1;                           // n is a root of assembly tree
  w[n]    = 0;                            // n is a dead element
  
  
  // ======================================================================  
  // Step 3: Initialize degree lists
  // ======================================================================  

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Step 3: Initialize degree lists" << endl;
  #endif

  // Loop over all nodes
  for (i = 0 ; i < n ; i++) {
  
    // Get degree of node
	d = degree[i];
	
	// Case A: Node i is empty
	if (d == 0) {
      elen[i] = -2;                       // element i is dead
      nel++ ;
      Cp[i]   = -1;                       // i is a root of assembly tree
      w[i]    = 0;                        // i is a dead element
	}
	
	// Case B: Node is dense
	else if (d > dense) {
      nv[i]   = 0;                        // absorb i into element n
      elen[i] = -1;                       // node i is dead
      nel++;
      Cp[i] = CS_FLIP(n);
      nv[n]++;
	}
	
	// Case C: Node is sparse
	else {
      if (head[d] != -1) last[head[d]] = i;
      next[i] = head[d];                  // put node i in degree list d
      head[d] = i;
	}
  } // endfor: looped over all columns


  // ======================================================================  
  // Step 4: Build permutations
  // ======================================================================  

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Step 4: Build permutations" << endl;
  #endif
  
  // While (selecting pivots) do
  while (nel < n) {
  
    // Select node of minimum approximate degree
	for (k = -1; mindeg < n && (k = head [mindeg]) == -1; mindeg++) ;
	if (next[k] != -1) last [next[k]] = -1;
	head[mindeg]  = next[k];              // remove k from degree list
	elenk         = elen[k];              // elenk = |Ek|
	nvk           = nv[k];                // # of nodes k represents
	nel          += nvk;                  // nv[k] nodes of A eliminated
	
	// Garbage collection. This frees the elbow rooms in matrix C
	if (elenk > 0 && cnz + mindeg >= nzmax) {
	
	  // Loop over all columns
      for (j = 0; j < n; j++) {
	    // If j is a live node or element, then ...
		if ((p = Cp[j]) >= 0) {
          Cp[j] = Ci[p];                  // save first entry of object
          Ci[p] = CS_FLIP(j);             // first entry is now CS_FLIP(j)
		}
	  }
	  
	  // Scan all memory
      for (q = 0, p = 0; p < cnz; ) {
	    // If we found object j, then ...
		if ((j = CS_FLIP(Ci[p++])) >= 0) {
          Ci[q] = Cp[j];                  // restore first entry of object
          Cp[j] = q++;                    // new pointer to object j
          for (k3 = 0; k3 < len[j]-1; k3++) 
		    Ci[q++] = Ci[p++];
		}
	  }
	  cnz = q;                            // Ci[cnz...nzmax-1] now free
	  
	} // endif: garbage collection done
	
	// Construct new element ---------------------------------------- */
	dk    = 0;
	nv[k] = -nvk;                         // flag k as in Lk
	p     = Cp[k];
	pk1   = (elenk == 0) ? p : cnz;       // do in place if elen[k] == 0
	pk2   = pk1;
	for (k1 = 1; k1 <= elenk + 1; k1++) {
	
	  // ...
      if (k1 > elenk) {
		e  = k;                           // search the nodes in k
		pj = p;                           // list of nodes starts at Ci[pj]
		ln = len[k] - elenk;              // length of list of nodes in k
	  }
	  else {
		e  = Ci[p++];                     // search the nodes in e
		pj = Cp[e];
		ln = len[e];                      // length of list of nodes in e
	  }
	  
	  // ...
	  for (k2 = 1 ; k2 <= ln ; k2++) {
		i = Ci[pj++];
		if ((nvi = nv [i]) <= 0)          // Skip of node i is dead or seen
		  continue;
		dk        += nvi;                 // degree[Lk] += size of node i
		nv[i]      = -nvi;                // negate nv[i] to denote i in Lk
		Ci[pk2++]  = i;                   // place i in Lk
		if (next[i] != -1) 
		  last[next[i]] = last[i];
		
		// remove i from degree list
		if (last[i] != -1)
          next[last[i]] = next[i];
		else
          head[degree[i]] = next[i];
	  } // endfor: k2
	  
	  // ...
	  if (e != k) {
		Cp[e] = CS_FLIP(k);               // absorb e into k
		w[e]  = 0;                        // e is now a dead element
	  }
	} // endfor: k1
	
	if (elenk != 0) cnz = pk2;            // Ci [cnz...nzmax] is free
	degree[k] = dk;                       // external degree of k - |Lk\i|
	Cp[k]     = pk1;                      // element k is in Ci[pk1..pk2-1]
	len[k]    = pk2 - pk1;
	elen[k]   = -2;                       // k is now an element
	
	// Clear w if necessary
	mark = cs_wclear(mark, lemax, w, n);
	
	// Scan 1: Find set differences |Le\Lk|
	for (pk = pk1 ; pk < pk2 ; pk++) {
      i = Ci [pk];
      if ((eln = elen[i]) <= 0)           // skip if elen[i] empty
	    continue;
	  nvi  = -nv[i];                      // nv [i] was negated
	  wnvi = mark - nvi;
	  
	  // Scan Ei
	  for (p = Cp[i]; p <= Cp[i] + eln - 1; p++) {
		e = Ci[p];
		if (w[e] >= mark)
		  w[e] -= nvi;                    // decrement |Le\Lk|
		else if (w[e] != 0)               // ensure e is a live element
		  w[e] = degree[e] + wnvi;        // 1st time e seen in scan 1
	  }
	} // endfor: scan 1 finished
	
	// Scan 2: Degree update
	for (pk = pk1; pk < pk2; pk++) {
      i  = Ci[pk];                        // consider node i in Lk
	  p1 = Cp[i];
	  p2 = p1 + elen [i] - 1;
	  pn = p1;
	  
	  // Scan Ei
	  for (h = 0, d = 0, p = p1 ; p <= p2 ; p++) {
		e = Ci[p];
		// If e is an unabsorbed element, then ...
		if (w[e] != 0) {
          dext = w[e] - mark;             // dext = |Le\Lk|
          if (dext > 0) {
			d += dext;                    // sum up the set differences
			Ci[pn++] = e;                 // keep e in Ei
			h += e;                       // compute the hash of node i
		  }
		  else {
			Cp[e] = CS_FLIP(k);            // aggressive absorb. e->k
			w[e]  = 0;                     // e is a dead element
		  }
		}
	  } // endfor: scanned Ei
	  
	  // elen[i] = |Ei|
	  elen[i] = pn - p1 + 1;
	  p3 = pn ;
	  p4 = p1 + len[i];
	  
	  // Prune edges in Ai
	  for (p = p2 + 1 ; p < p4 ; p++) {
		j = Ci[p];
		if ((nvj = nv[j]) <= 0)            // node j dead or in Lk
		  continue;
		d += nvj;                          // degree(i) += |j|
		Ci[pn++] = j;                      // place j in node list of i
		h += j;                            // compute hash for node i
	  }
	  
	  // Check for mass elimination
	  if (d == 0) {
		Cp[i]   = CS_FLIP(k);              // absorb i into k
		nvi     = -nv[i];
		dk     -= nvi;                     // |Lk| -= |i|
		nvk    += nvi;                     // |k| += nv[i]
		nel    += nvi;
		nv[i]   = 0;
		elen[i] = -1;                      // node i is dead
	  }
	  else {
		degree[i] = CS_MIN(degree[i], d);   // update degree(i)
		Ci[pn] = Ci[p3];                    // move first node to end
		Ci[p3] = Ci[p1];                    // move 1st el. to end of Ei
		Ci[p1] = k;                         // add k as 1st element in of Ei
		len[i] = pn - p1 + 1;               // new len of adj. list of node i
		h %= n;                             // finalize hash of i
		next[i]  = hhead[h];                // place i in hash bucket
		hhead[h] = i;
		last[i]  = h;                       // save hash of i in last[i]
	  }
	} // endfor: scan 2 is done
	
	// finalize |Lk|
	degree[k] = dk;
	lemax     = CS_MAX(lemax, dk);
	
	// Clear w
	mark = cs_wclear(mark+lemax, lemax, w, n);
	
	// Supernode detection
	for (pk = pk1 ; pk < pk2 ; pk++) {
      i = Ci[pk];
	  if (nv[i] >= 0) continue;             // skip if i is dead
	  h = last[i];                          // scan hash bucket of node i
	  i = hhead[h];
	  hhead[h] = -1;                        // hash bucket will be empty
	  
	  // ...
	  for ( ; i != -1 && next[i] != -1 ; i = next[i], mark++) {
		ln  = len[i];
		eln = elen[i];
		for (p = Cp[i]+1 ; p <= Cp[i] + ln-1 ; p++) 
		  w[Ci[p]] = mark;
		jlast = i;
		
		// Compare i with all j
		for (j = next[i]; j != -1; )	 {
          ok = (len[j] == ln) && (elen[j] == eln) ;
		  for (p = Cp[j] + 1; ok && p <= Cp[j] + ln - 1; p++) {
			if (w[Ci[p]] != mark)           // compare i and j
			  ok = 0;
		  }
		  
		  // i and j are identical
		  if (ok) {
			Cp[j]       = CS_FLIP(i);       // absorb j into i
			nv[i]      += nv[j];
			nv[j]       = 0;
			elen[j]     = -1;               // node j is dead
			j           = next[j];          // delete j from hash bucket
			next[jlast] = j;
		  }
		  // j and i are different
		  else {
			jlast = j;
			j     = next[j];
		  }
		} // endfor: comparison of i with all j
	  } //
	} // endfor: supernode detection
	
	// Finalize new element Lk
	for (p = pk1, pk = pk1; pk < pk2; pk++) {
      i = Ci[pk];
      if ((nvi = -nv[i]) <= 0)              // skip if i is dead
	    continue;
      nv[i] = nvi;                          // restore nv[i]
      d = degree[i] + dk - nvi;             // compute external degree(i)
	  d = CS_MIN(d, n - nel - nvi);
      if (head[d] != -1) last[head[d]] = i;
      next [i] = head [d];                  // put i back in degree list
      last [i] = -1 ;
      head [d] = i ;
      mindeg = CS_MIN (mindeg, d);          // find new minimum degree
      degree [i] = d ;
      Ci [p++] = i ;                        // place i in Lk
	}
	
	// Number of nodes absorbed into k
	nv[k] = nvk;
	
	// ...
	if ((len [k] = p-pk1) == 0) {           // length of adj list of element k
	    Cp [k] = -1 ;                       // k is a root of the tree
	    w [k] = 0 ;                         // k is now a dead element
	}
	if (elenk != 0) cnz = p ;               // free unused space in Lk

  } // endwhile


  // ======================================================================  
  // Step 5: Postordering
  // ======================================================================  

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Step 5: Postordering" << endl;
  #endif
	
  // Fix assembly tree
  for (i = 0; i < n; i++) 
    Cp[i] = CS_FLIP(Cp[i]);
  for (j = 0 ; j <= n ; j++) 
    head[j] = -1;
	
  // Place unordered nodes in lists
  for (j = n ; j >= 0 ; j--) {
    if (nv[j] > 0) continue;                // skip if j is an element
	next[j] = head [Cp[j]];                 // place j in list of its parent
	head[Cp[j]] = j;
  }
	
  // Place elements in lists
  for (e = n ; e >= 0 ; e--) {
	if (nv[e] <= 0) continue;               // skip unless e is an element
	if (Cp[e] != -1) {
	  next[e] = head[Cp[e]];                // place e in list of its parent
	  head[Cp[e]] = e;
	}
  }


  // ======================================================================  
  // Step 6: Postorder the assembly tree
  // ======================================================================  

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Step 6: Postorder the assembly tree" << endl;
  #endif
  
  // Postorder the assembly tree
  for (k = 0, i = 0 ; i <= n ; i++) {
	if (Cp [i] == -1) 
	  k = cs_tdfs(i, k, head, next, P, w);
  }
  
  // Free workspace
  delete [] wrk_int;

  // Debug
  #if defined(G_DEBUG_SPARSE_SYM_AMD_ORDERING)
  cout << " Permutation:" << endl;
  for (i = 0; i < n; ++i)
    cout << " " << P[i];
  cout << endl;
  cout << "GSparseSymbolic::cs_amd finished" << endl;
  #endif
  
  // Return result
  return P;
}


/***************************************************************************
 *                                  cs_counts                              *
 * ----------------------------------------------------------------------- *
 * Column counts of LL'=A or LL'=A'A, given parent & post ordering         *
 * ----------------------------------------------------------------------- *
 * Input:   A                    Sparse matrix A                           *
 *          parent[0..n-1]       parent array                              *
 *          post[0..n-1]         post ordering array                       *
 *          ata                  0: count LL'=A, 1: count LL'=A'A          *
 * Output:  p[0..n]              integer array of n+1 elements             *
 ***************************************************************************/
#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)
int* GSparseSymbolic::cs_counts(const GSparseMatrix* A, const int* parent, 
                                const int* post, int ata)
{
    // Declare loop variables
    int i, j, k, J, p, q, jleaf;

    // Return NULL pointer if input arrays are NULL or sparse matrix is
    // empty
    if (!parent || !post) return NULL;

    // Assign input matrix attributes
    int m = A->m_rows;
    int n = A->m_cols;
  
    // Allocate result
    int* colcount = new int[n];
    int* delta = colcount;

    // Allocate workspace
    int  wrk_size = 4*n + (ata ? (n+m+1) : 0);
    int* wrk_int  = new int[wrk_size];

  // Get (logical) transpose of A: AT = A'
  GSparseMatrix AT = cs_transpose(*A, 0);

  // Set-up pointers to workspace
  int* ancestor = wrk_int;
  int* maxfirst = wrk_int + n;
  int* prevleaf = wrk_int + 2*n;
  int* first    = wrk_int + 3*n;

  // Clear workspace
  for (k = 0; k < wrk_size; k++) 
    wrk_int[k] = -1;

  // Find first j
  for (k = 0; k < n; k++) {
	j        = post[k];
	delta[j] = (first[j] == -1) ? 1 : 0;    // delta[j]=1 if j is a leaf
	for ( ; j != -1 && first[j] == -1; j = parent[j]) first[j] = k;
  }

  // Assign output matrix attributes
  int* ATp = AT.m_colstart;
  int* ATi = AT.m_rowinx; 

  // Initialise. Note that init_ata wires the working array into the *head
  // and *next pointers, so both arrays are just aliases of the wrk_int
  // array.
  int* head = NULL;
  int* next = NULL;
  if (ata) init_ata(&AT, post, wrk_int, &head, &next);

  // Each node in its own set  
  for (i = 0; i < n; i++) ancestor[i] = i;
  
  // Loop over all columns of matrix A (or rows of matrix AT)
  for (k = 0; k < n; k++) {
	j = post[k];		// j is the kth node in postordered etree
	if (parent[j] != -1) delta[parent[j]]--;        // j is not a root
	for (J = HEAD(k,j); J != -1; J = NEXT(J)) {     // J=j for LL'=A case
	
	  // Loop over all non-zero elements in column J
      for (p = ATp[J]; p < ATp[J+1]; p++) {
		i = ATi [p] ;
		q = cs_leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
		if (jleaf >= 1) delta [j]++ ;   /* A(i,j) is in skeleton */
		if (jleaf == 2) delta [q]-- ;	/* account for overlap in q */
	  }
	}
	if (parent[j] != -1) ancestor[j] = parent[j];
  }
  
  // Sum up delta's of each child
  for (j = 0; j < n; j++) {
	if (parent[j] != -1) colcount[parent[j]] += colcount[j];
  }
  
  // Free workspace
  delete [] wrk_int;
  
  // Return result
  return colcount;
} 


/***************************************************************************
 *                                  cs_etree                               *
 * ----------------------------------------------------------------------- *
 * Compute the etree of A using triu(A), or A'A without forming A'A        *
 * ----------------------------------------------------------------------- *
 * Input:   A              Sparse matrix A                                 *
 *          ata            flag                                            *
 * Output:  etree(A)       Evaluation tree of matrix A                     *
 ***************************************************************************/
int* GSparseSymbolic::cs_etree(const GSparseMatrix* A, int ata)
{
    // Declare loop variables
    int i, k, p;

    // Assign matrix attributes
    int m = A->m_rows;
    int n = A->m_cols;

    // Allocate result array
    int* parent = new int[n];

    // Allocate workspace
    int  wrk_size = n + (ata ? m : 0);
    int* wrk_int  = new int[wrk_size];

  // Set-up pointers to workspace. If 'ata=1' then we use also a prev array
  int* ancestor = wrk_int;
  int* prev     = wrk_int + n;
  
  // If 'prev' is requested, initialise array with -1
  if (ata) {
    for (i = 0; i < m ; ++i) 
	  prev[i] = -1;
  }

  // Loop over all nodes (columns)
  for (k = 0; k < n; ++k) {
  
    // Initialise node k to having no parent and no ancestor
	parent[k]   = -1;		    // node k has no parent yet
	ancestor[k] = -1;		    // nor does k have an ancestor
	
	// Loop over all elements in column k
	for (p = A->m_colstart[k]; p < A->m_colstart[k+1]; ++p) {
	
	  // Get element index
	  i = ata ? (prev[A->m_rowinx[p]]) : (A->m_rowinx[p]);
	  
	  // Traverse from i to k
	  int inext;
	  for ( ; i != -1 && i < k; i = inext) {
		inext       = ancestor[i];        // inext = ancestor of i
		ancestor[i] = k;                  // path compression
		if (inext == -1) parent[i] = k;   // no ancestor, parent is k
	  }
	  if (ata) prev[A->m_rowinx[p]] = k;
	}
  }

  // Free workspace
  delete [] wrk_int;
  
  // Return result
  return parent;

}


/***********************************************************************//**
 * @brief cs_fkeep
 *
 * @param[in] A Sparse matrix.
 * @param[in] fkeep Screening function.
 * @param[in] other Other function.
 *
 * Drop entries for which fkeep(A(i,j)) is false. Return the number of
 * non-zero elements if ok, otherwise return -1.
 ***************************************************************************/
int GSparseSymbolic::cs_fkeep(GSparseMatrix* A, 
                              int(*fkeep)(int, int, double, void*), 
                              void* other)
{
    // Return error if some of the input pointers is invalid
    if (!A || !fkeep) return (-1);
  
    // Initialise counter for non-zero elements 
    int nz = 0;

    // Assign matrix attributes
    int     n  = A->m_cols;
    int*    Ap = A->m_colstart;
    int*    Ai = A->m_rowinx; 
    double* Ax = A->m_data;

    // Operate only if we have some data
    if (Ap != NULL && Ai != NULL && Ax != NULL) {

        // Loop over all columns
        for (int j = 0; j < n; j++) {
            int p = Ap[j];                     // get current location of col j
            Ap[j] = nz;                        // record new location of col j
            for ( ; p < Ap[j+1] ; p++) {
                if (fkeep(Ai[p], j, Ax ? Ax[p] : 1, other)) {
                    if (Ax) Ax[nz] = Ax[p];    // keep A(i,j)
                    Ai[nz++] = Ai[p] ;
                }
            }
        }
  
        // Finalise A
        Ap[n]         = nz;
        A->m_elements = nz;

    } // endif: we had data
  
    // Remove extra space from A
//  cs_sprealloc(A, 0);
    A->free_elements(nz, (A->m_elements-nz));

    // Return number of non-zero elements
    return nz;
}


/***************************************************************************
 *                                   cs_leaf                               *
 * ----------------------------------------------------------------------- *
 * Consider A(i,j), node j in ith row subtree and return lca(jprev,j)      *
 * ----------------------------------------------------------------------- *
 * Input:   i                    row to consider                           *
 *          j                    column to consider                        *
 *          first[0..n-1]        ...                                       *
 *          maxfirst[0..n-1]     ...                                       *
 *          prevleaf[0..n-1]     ...                                       *
 *          ancestor[0..n-1]     ...                                       *
 * Output:  q                    ...                                       *
 *          jleaf                ...                                       *
 ***************************************************************************/
int GSparseSymbolic::cs_leaf(int i, int j, const int *first, int *maxfirst, 
                             int *prevleaf, int *ancestor, int *jleaf)
{
  // If one of the input pointers is invalid then return -1
  if (!first || !maxfirst || !prevleaf || !ancestor || !jleaf) return (-1);

  // Declare loop variables
  int q, s, sparent;

  // ...    
  *jleaf = 0 ;
  
  // If j is not a leaf then return -1
  if (i <= j || first[j] <= maxfirst[i]) return (-1);

  // Update max first[j] seen so far
  maxfirst[i] = first[j];
  
  // jprev = previous leaf of ith subtree
  int jprev = prevleaf[i];
  prevleaf[i] = j;
  
  // j is first or subsequent leaf
  *jleaf = (jprev == -1) ? 1: 2;
  
  // If j is 1st leaf then q is the root of ith subtree
  if (*jleaf == 1) return (i);
  
  // ...
  for (q = jprev; q != ancestor[q]; q = ancestor[q]) ;
  
  // ...
  for (s = jprev; s != q; s = sparent) {
	sparent     = ancestor[s];     // path compression
	ancestor[s] = q;
  }
  
  // Return least common ancester (jprev,j)
  return (q);
}


/***************************************************************************
 *                                   cs_pinv                               *
 * ----------------------------------------------------------------------- *
 * Inverts the permutation p[0..n-1]                                       *
 * ----------------------------------------------------------------------- *
 * Input:   p[0..n-1]      integer array of n elements (NULL denotres ID)  *
 * Output:  pinv[0..n-1]   pinv[0..n-1] = p[n-1..0]                        *
 ***************************************************************************/
int* GSparseSymbolic::cs_pinv(int const *p, int n)
{
    // Return NULL pointer if input pointer is NULL. This denotes identity
    if (!p || n < 1) return NULL;
  
    // Allocate result array
    int* pinv = new int[n];
  
    // Invert the permutation
    for (int k = 0; k < n; ++k) 
        pinv[p[k]] = k;
  
    // Return result
    return pinv;
}


/***************************************************************************
 *                                   cs_post                               *
 * ----------------------------------------------------------------------- *
 * Post order a forest                                                     *
 * ----------------------------------------------------------------------- *
 * Input:   parent[0..n-1]    integer array of n elements                  *
 * Output:  post[0..n-1]      ...                                          *
 ***************************************************************************/
int* GSparseSymbolic::cs_post(const int* parent, int n)
{
    // Return NULL pointer if input pointer is NULL
    if (!parent || n < 1) return (NULL);

    // Allocate result array
    int* post = new int[n];

    // Allocate workspace
    int  wrk_size = 3 * n;
    int* wrk_int  = new int[wrk_size];

  // Set-up pointers to workspace
  int* head  = wrk_int;
  int* next  = wrk_int + n;
  int* stack = wrk_int + 2*n;  // Working array for 'cs_tdfs' function

  // Declare and initialise loop variables
  int j, k = 0;
  
  // Empty linked list
  for (j = 0 ; j < n ; j++) head[j] = -1;
  
  // Traverse nodes in reverse order to build a linked list 'head'
  for (j = n-1; j >= 0; j--) {
  
    // Skip if j is a root
	if (parent[j] == -1) continue;
	
	// Add j to list of its parent
	next[j]         = head[parent[j]];
	head[parent[j]] = j;
  }
  
  // Traverse nodes
  for (j = 0 ; j < n ; j++) {

    // Skip if j is not a root
	if (parent[j] != -1) continue;
	
	// Depth-first search and postorder of a tree rooted at node j
	k = cs_tdfs(j, k, head, next, post, stack);
  }
  
  // Free workspace
  delete [] wrk_int;
  
  // Return result
  return post;

}


/***************************************************************************
 *                                   cs_tdfs                               *
 * ----------------------------------------------------------------------- *
 * Depth-first search and postorder of a tree rooted at node j             *
 * ----------------------------------------------------------------------- *
 * Input:   j                 Node of tree root                            *
 *          k                 ...                                          *
 *          head              head[0..n-1] linked list (updated)           *
 *          next              next[0..n-1] linked list info                *
 *          post              post[0..n-1] array (updated)                 *
 *          stack             stack[0..n-1] working array (updated)        *
 * Output:  k                 Update of k                                  *
 ***************************************************************************/
int GSparseSymbolic::cs_tdfs(int j, int k, int* head, const int* next, 
                             int* post, int* stack)
{
  // Return -1 if one of the array pointers is NULL
  if (!head || !next || !post || !stack) return (-1);

  // Decrare
  int i, p;

  // Place j on the first position of the stack
  stack[0] = j;
  
  // Loop while stack is not empty
  int top = 0;
  while (top >= 0) {
	p = stack[top];              // p = top of stack
	i = head[p];                 // i = youngest child of p
	if (i == -1) {
      top--;                     // p has no unordered children left
      post[k++] = p;             // node p is the kth postordered node
	}
	else {
      head[p]      = next[i];    // remove i from children of p
	  stack[++top] = i;          // start depth-first search on child node i
	}
  }
  
  // Return result
  return k;
}


/*==========================================================================
 =                                                                         =
 =                  GSparseSymbolic static member functions                =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                  init_ata                               *
 * ----------------------------------------------------------------------- *
 * Initialise A'A.                                                         *
 * Note that the returned pointers 'head' and 'info' point into the        *
 * 'wrk_int' workspace array.                                              *
 * ----------------------------------------------------------------------- *
 * Input:   AT                Sparse matrix A'                             *
 *          post              ...                                          *
 *          wrk_int           Workspace                                    *
 *          head              pointer to head[0..n-1] linked list          *
 *          next              pointer to next[0..n-1] linked list info     *
 * Output:  -                 -                                            *
 ***************************************************************************/
void GSparseSymbolic::init_ata(const GSparseMatrix* AT, const int* post, 
                               int* wrk_int, int** head, int** next)
{
  // Declare loop variables
  int i, k, p;

  // Assign matrix attributes
  int  m   = AT->m_cols;
  int  n   = AT->m_rows;
  int* ATp = AT->m_colstart;
  int* ATi = AT->m_rowinx; 
  
  // Set-up pointers to workspace
  *head = wrk_int + 4*n;
  *next = wrk_int + 5*n+1;

  // Invert post
  for (k = 0 ; k < n ; k++) 
    wrk_int[post[k]] = k;
	
  // Loop over columns of AT matrix
  for (i = 0 ; i < m ; i++) {
  
    // Loop over elements in column i of AT matrix
	for (k = n, p = ATp[i]; p < ATp[i+1]; p++) 
	  k = CS_MIN(k,wrk_int[ATi[p]]);
	  
	// Place row i in linked list k
	(*next)[i] = (*head)[k];
	(*head)[k] = i;
  }
  
  // Return
  return;
}


/***************************************************************************
 *                                  cs_diag                                *
 * ----------------------------------------------------------------------- *
 * Keep off-diagonal entries and drop diagonal entries.                    *
 * ----------------------------------------------------------------------- *
 * Input:   i                 Row index                                    *
 *          j                 Column index                                 *
 *          aij               NOT USED                                     *
 *          other             NOT USED                                     *
 * Output:  diag              0: off-diagonal, 1: diagonal element         *
 ***************************************************************************/
int GSparseSymbolic::cs_diag(int i, int j, double aij, void* other) 
{ 
  return (i != j); 
}


/***************************************************************************
 *                                 cs_wclear                               *
 * ----------------------------------------------------------------------- *
 * Clear w array                                                           *
 * ----------------------------------------------------------------------- *
 * Input:   mark              ...                                          *
 *          lemax             ...                                          *
 *          w                 w[0..n-1] working array                      *
 *          n                 Size of the working array                    *
 * Output:  mark              Update of mark                               *
 *          w                 Cleared w[0..n-1] working array              *
 ***************************************************************************/
int GSparseSymbolic::cs_wclear(int mark, int lemax, int* w, int n)
{
    // ...
    if (mark < 2 || (mark + lemax < 0)) {
	    for (int k = 0; k < n; k++) if (w[k] != 0) w[k] = 1;
        mark = 2;
    }
  
    // Return mark. At this point, w[0..n-1] < mark holds
    return mark;
}


/*==========================================================================
 =                                                                         =
 =                          GSparseSymbolic friends                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GSparseSymbolic& s)
{
  // Put header is stream
  os << "=== GSparseSymbolic ===";
  
  // Show inverse permutation if it exists
  if (s.m_pinv != NULL) {
    os << std::endl << " Inverse permutation (of " << s.m_n_pinv << " elements)" << 
	   std::endl;
	for (int i = 0; i < s.m_n_pinv; ++i)
	  os << " " << s.m_pinv[i];
  } 

  // Show fill reducing column permutation if it exists
  if (s.m_q != NULL) {
    os << std::endl << " Fill reducing column permutation (" << s.m_n_cp << 
	   " elements)" << std::endl;
	for (int i = 0; i < s.m_n_cp; ++i)
	  os << " " << s.m_q[i];
  } 

  // Show elimination tree for Cholesky and QR if it exists
  if (s.m_parent != NULL) {
    os << std::endl << " Elimination tree (" << s.m_n_parent << " elements)" << 
	   std::endl;
	for (int i = 0; i < s.m_n_parent; ++i)
	  os << " " << s.m_parent[i];
  } 

  // Show column pointers for Cholesky row counts for QR if it exists
  if (s.m_cp != NULL) {
    os << std::endl << " Cholesky column pointers/QR row counts (" << s.m_n_cp << 
	   " elements)" << std::endl;
	for (int i = 0; i < s.m_n_cp; ++i)
	  os << " " << s.m_cp[i];
  } 

  // Show leftmost[i] = min(find(A(i,:))) if it exists
  if (s.m_leftmost != NULL) {
    os << std::endl << " leftmost[i] = min(find(A(i,:))) (" << s.m_n_leftmost << 
	   " elements)" << std::endl;
	for (int i = 0; i < s.m_n_leftmost; ++i)
	  os << " " << s.m_leftmost[i];
  } 
  
  // Return output stream
  return os;
}
