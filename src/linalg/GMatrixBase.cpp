/***************************************************************************
 *               GMatrixBase.cpp  -  matrix abstract base class            *
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

/* __ Includes ___________________________________________________________ */
#include "GMatrixBase.hpp"
#include "GException.hpp"
#include "GVector.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Method name definitions ____________________________________________ */
#define G_COPY_MEMBERS    "GMatrixBase::copy_members(const GMatrixBase&)"
#define G_SELECT_NON_ZERO "GMatrixBase::select_non_zero()"


/*==========================================================================
 =                                                                         =
 =                    GMatrixBase constructors/destructors                 =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                               Constructor                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GMatrixBase::GMatrixBase()
{
  // Initialise class members for clean destruction
  init_members();

  // Return
  return;
}


/***************************************************************************
 *                              Copy constructor                           *
 * ----------------------------------------------------------------------- *
 * The copy constructor is sufficiently general to provide the base        *
 * constructor for all derived classes, including sparse matrices.         *
 ***************************************************************************/
GMatrixBase::GMatrixBase(const GMatrixBase& m)
{
  // Initialise class members for clean destruction
  init_members();

  // Copy members
  copy_members(m);

  // Return
  return;
}


/***************************************************************************
 *                                Destructor                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GMatrixBase::~GMatrixBase()
{
  // Free class members
  free_members();

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                           GMatrixBase operators                         =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GMatrixBase& GMatrixBase::operator= (const GMatrixBase& m)
{
  // Execute only if object is not identical
  if (this != &m) {
  
    // Free members
    free_members();
  
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(m);
	
  } // endif: object was not identical

  // Return this object
  return *this;
}


/***************************************************************************
 *                              Equalty operator                           *
 * ----------------------------------------------------------------------- *
 * Two matrixes are considered equal if they have the same dimensions and  *
 * identicial elements.                                                    *
 ***************************************************************************/
int GMatrixBase::operator== (const GMatrixBase& m) const
{
  // Initalise the result to 'equal matrices'
  int result = 1;
  
  // Perform comparison (only if matrix dimensions and number of physical
  // elements are identical; otherwise set result to 'false')
  if (m_rows == m.m_rows && m_cols == m.m_cols && m_elements == m.m_elements) {
    for (int i = 0; i < m_elements; ++i) {
	  if (m_data[i] != m.m_data[i]) {
	    result = 0;
		break;
	  }
	}
  }
  else
    result = 0;

  // Return result
  return result;
}


/***************************************************************************
 *                           Non-equality operator                         *
 * ----------------------------------------------------------------------- *
 * Two matrixes are considered as non-equal (or different) if the differ   *
 * in their dimensions or if at least one element differs.                 *
 ***************************************************************************/
int GMatrixBase::operator!= (const GMatrixBase& m) const
{
  // Initalise the result to 'equal matrices'
  int result = 0;

  // Perform comparison (only if matrix dimensions and number of physical
  // elements are identical; otherwise set result to 'true')
  if (m_rows == m.m_rows && m_cols == m.m_cols && m_elements == m.m_elements) {
    for (int i = 0; i < m_elements; ++i) {
	  if (m_data[i] != m.m_data[i]) {
	    result = 1;
		break;
	  }
	}
  }
  else
    result = 1;
	
  // Return result
  return result;
}


/*==========================================================================
 =                                                                         =
 =                       GMatrixBase protected methods                     =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Negate matrix elements                          *
 * ----------------------------------------------------------------------- *
 * Inverts the sign of all matrix elements.                                *
 ***************************************************************************/
void GMatrixBase::negation(void)
{
  // Inverts the sign for all matrix elements
  for (int i = 0; i < m_elements; ++i)
    m_data[i] = -m_data[i];
	
  // Return
  return;
}


/***************************************************************************
 *                           Add matrix elements                           *
 * ----------------------------------------------------------------------- *
 * Add all matrix elements.                                                *
 ***************************************************************************/
void GMatrixBase::addition(const GMatrixBase& m)
{
  // Add all elements of matrix
  for (int i = 0; i < m_elements; ++i)
    m_data[i] += m.m_data[i];
	
  // Return
  return;
}


/***************************************************************************
 *                         Subtract matrix elements                        *
 * ----------------------------------------------------------------------- *
 * Subtract all matrix elements.                                           *
 ***************************************************************************/
void GMatrixBase::subtraction(const GMatrixBase& m)
{
  // Add all elements of matrix
  for (int i = 0; i < m_elements; ++i)
    m_data[i] -= m.m_data[i];
	
  // Return
  return;
}


/***************************************************************************
 *                    Multiply matrix elements with scalar                 *
 * ----------------------------------------------------------------------- *
 * Multiply all matrix elements with a scalar. There are three cases:      *
 * (A) the multiplier is 0, then reset all elements to 0                   *
 * (B) the multiplier is +/-1, then do nothing or negate                   *
 * (C) in any other case, multiply by multiplier                           *
 ***************************************************************************/
void GMatrixBase::multiplication(const double& s)
{
  // Case A: If multiplicator is 0 then set entire matrix to 0
  if (s == 0.0) {
    for (int i = 0; i < m_elements; ++i)
      m_data[i] = 0.0;
  }
  
  // Case C: If multiplicator is not +/- 1 then perform multiplication
  else if (fabs(s) != 1.0) {
    for (int i = 0; i < m_elements; ++i)
      m_data[i] *= s;
  }
  
  // Case B: If multiplication is -1 then negate
  else if (s == -1.0) {
    for (int i = 0; i < m_elements; ++i)
      m_data[i] = -m_data[i];
  }
	
  // Return
  return;
}


/***************************************************************************
 *                    Set all elements to a specific value                 *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrixBase::set_all_elements(const double& s)
{
  // Add all elements of matrix
  for (int i = 0; i < m_elements; ++i)
    m_data[i] = s;
	
  // Return
  return;
}


/***************************************************************************
 *                      Returns minimum matrix element                     *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GMatrixBase::get_min_element() const
{
  // Initialise minimum with first element
  double result = m_data[0];
  
  // Search all elements for the smallest one
  for (int i = 1; i < m_elements; ++i) {
    if (m_data[i] < result)
	  result = m_data[i];
  }
  
  // Return result
  return result;
}


/***************************************************************************
 *                      Returns maximum matrix element                     *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GMatrixBase::get_max_element() const
{
  // Initialise minimum with first element
  double result = m_data[0];
  
  // Search all elements for the smallest one
  for (int i = 1; i < m_elements; ++i) {
    if (m_data[i] > result)
	  result = m_data[i];
  }
  
  // Return result
  return result;
}


/***************************************************************************
 *                         Sum all matrix elements                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GMatrixBase::get_element_sum() const
{
  // Initialise matrix sum
  double result = 0.0;
  
  // Add all elements  
  for (int i = 0; i < m_elements; ++i)
    result += m_data[i];
  
  // Return result
  return result;
}


/***************************************************************************
 *                     Dump all matrix elements in ostream                 *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrixBase::dump_elements(ostream& os) const
{
  // Dump matrix elements row by row using the access function
  for (int row = 0; row < m_rows; ++row) {
    os << " ";
    for (int col = 0; col < m_cols; ++col) {
      os << (*this)(row,col);
	  if (col != m_cols-1)
	    os << ", ";
	}
	if (row != m_rows-1)
	  os << endl;
  }
  
  // Return
  return;
}


/***************************************************************************
 *           Dump row compression scheme in ostream if it exists           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrixBase::dump_row_comp(ostream& os) const
{
  // If there is a row compression the show scheme
  if (m_rowsel != NULL) {
    os << endl << " Row selection ..:";
    for (int row = 0; row < m_num_rowsel; ++row)
      os << " " << m_rowsel[row];
  }
  
  // Return
  return;
}


/***************************************************************************
 *           Dump column compression scheme in ostream if it exists        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrixBase::dump_col_comp(ostream& os) const
{
  // If there is a column compression the show scheme
  if (m_colsel != NULL) {
    os << endl << " Column selection:";
    for (int col = 0; col < m_num_colsel; ++col)
      os << " " << m_colsel[col];
  }
  
  // Return
  return;
}


/***************************************************************************
 *                              select_non_zero                            *
 * ----------------------------------------------------------------------- *
 * Determines the non-zero rows and columns in matrix and set up index     *
 * arrays that points to these rows/columns. These arraya are used for     *
 * compressed matrix factorisations.                                       *
 ***************************************************************************/
void GMatrixBase::select_non_zero(void)
{
  // Free existing selection arrays
  if (m_rowsel != NULL) delete [] m_rowsel;
  if (m_colsel != NULL) delete [] m_colsel;
  
  // Allocate selection arrays
  m_rowsel = new int[m_rows];
  if (m_rowsel == NULL)
    throw GException::mem_alloc(G_SELECT_NON_ZERO, m_rows);
  m_colsel = new int[m_cols];
  if (m_colsel == NULL)
    throw GException::mem_alloc(G_SELECT_NON_ZERO, m_cols);

  // Initialise non-zero row and column counters
  m_num_rowsel = 0;
  m_num_colsel = 0;
  
  // Declare loop variables
  int row;
  int col;
  
  // Find all non-zero rows
  for (row = 0; row < m_rows; ++row) {
    for (col = 0; col < m_cols; ++col) {
	  if ((*this)(row,col) != 0.0)
	    break;
	}
	if (col < m_cols)      // found a non-zero element in row
      m_rowsel[m_num_rowsel++] = row;
  }

  // Find all non-zero columns
  for (col = 0; col < m_cols; ++col) {
    for (row = 0; row < m_rows; ++row) {
	  if ((*this)(row,col) != 0.0)
	    break;
	}
	if (row < m_rows)      // found a non-zero element in column
      m_colsel[m_num_colsel++] = col;
  }
  
  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                        GMatrixBase private methods                      =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrixBase::init_members(void)
{
  // Initialise GMatrixBase members
  m_rows       = 0;      // Number of rows
  m_cols       = 0;      // Number of columns
  m_elements   = 0;      // Logically used number of elements
  m_alloc      = 0;      // Allocated # of elements (>= m_elements)
  m_num_rowsel = 0;      // Number of selected rows (for comp. decomp.)
  m_num_colsel = 0;      // Number of selected columns (for comp. decomp.)
  m_data       = NULL;   // Matrix data
  m_colstart   = NULL;   // Column start indices (m_cols+1)
  m_rowsel     = NULL;   // Row selection (for compressed decomposition)
  m_colsel     = NULL;   // Column selection (for compressed decomposition)
  
  // Return
  return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrixBase::copy_members(const GMatrixBase& m)
{
  // Copy matrix attributes
  m_rows       = m.m_rows;
  m_cols       = m.m_cols;
  m_elements   = m.m_elements;
  m_alloc      = m.m_alloc;
  m_num_rowsel = m.m_num_rowsel;
  m_num_colsel = m.m_num_colsel;

  // Allocate memory for column start array and copy content
  m_colstart = new int[m_cols+1];
  if (m_colstart == NULL)
    throw GException::mem_alloc(G_COPY_MEMBERS, m_cols+1);
  for (int i = 0; i <= m_cols; ++i)
    m_colstart[i] = m.m_colstart[i];
  
  // Allocate memory for elements and copy them (only if there are elements)
  if (m.m_data != NULL && m_alloc > 0) {
    m_data = new double[m_alloc];
    if (m_data == NULL)
	  throw GException::mem_alloc(G_COPY_MEMBERS, m_alloc);
    for (int i = 0; i < m_elements; ++i)
      m_data[i] = m.m_data[i];
  }
  
  // If there is a row selection then copy it
  if (m.m_rowsel != NULL && m_num_rowsel > 0) {
    m_rowsel = new int[m_num_rowsel];
    if (m_rowsel == NULL)
	  throw GException::mem_alloc(G_COPY_MEMBERS, m_num_rowsel);
    for (int i = 0; i < m_num_rowsel; ++i)
      m_rowsel[i] = m.m_rowsel[i];
  }

  // If there is a column selection then copy it
  if (m.m_colsel != NULL && m_num_colsel > 0) {
    m_colsel = new int[m_num_colsel];
    if (m_colsel == NULL)
	  throw GException::mem_alloc(G_COPY_MEMBERS, m_num_colsel);
    for (int i = 0; i < m_num_colsel; ++i)
      m_colsel[i] = m.m_colsel[i];
  }

  // Return
  return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GMatrixBase::free_members(void)
{
  // De-allocate only if memory has indeed been allocated
  if (m_colsel   != NULL) delete [] m_colsel;
  if (m_rowsel   != NULL) delete [] m_rowsel;
  if (m_data     != NULL) delete [] m_data;
  if (m_colstart != NULL) delete [] m_colstart;

  // Properly mark members as free
  m_colstart = NULL;
  m_data     = NULL;
  m_rowsel   = NULL;
  m_colsel   = NULL;
  
  // Return
  return;
}
