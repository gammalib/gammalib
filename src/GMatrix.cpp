/***************************************************************************
 *                       GMatrix.cpp  -  matrix class                      *
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
#include "GMatrix.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                           GMatrix constructor                           *
 ***************************************************************************/
GMatrix::GMatrix(unsigned rows, unsigned cols)
{
  // Initialise private members for clean destruction
  m_rows     = 0;
  m_cols     = 0;
  m_elements = 0;
  m_data     = NULL;
  m_colstart = NULL;

  // Determine number of physical elements in matrix
  unsigned elements = rows*cols;

  // Throw exception if requested matrix size is zero
  if (elements == 0)
    throw empty("GMatrix constructor");
	
  // Allocate matrix array and column start index array. Throw an exception 
  // if allocation failed
  m_data     = new double[elements];
  m_colstart = new unsigned[cols+1];
  if (m_data == NULL || m_colstart == NULL)
	throw mem_alloc("GMatrix constructor", elements);
	
  // Store matrix size (logical and physical)
  m_rows     = rows;
  m_cols     = cols;
  m_elements = elements;
  
  // Set-up column start indices
  m_colstart[0] = 0;
  for (unsigned col = 1; col <= m_cols; ++col)
    m_colstart[col] = m_colstart[col-1] + m_rows;
	
  // Initialise matrix elements to 0.0
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] = 0.0;
}


/***************************************************************************
 *                         GMatrix copy constructor                        *
 ***************************************************************************/
GMatrix::GMatrix(const GMatrix& m)
{
  m_rows     = 0;
  m_cols     = 0;
  m_elements = 0;
  m_data     = NULL;
  m_colstart = NULL;
  m_data     = new double[m.m_elements];
  m_colstart = new unsigned[m.m_cols+1];
  if (m_data == NULL || m_colstart == NULL)
	throw mem_alloc("GMatrix copy constructor", m.m_elements);
  m_rows     = m.m_rows;
  m_cols     = m.m_cols;
  m_elements = m.m_elements;
  for (unsigned i = 0; i <= m_cols; ++i)
    m_colstart[i] = m.m_colstart[i];
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] = m.m_data[i];
}


/***************************************************************************
 *                           GMatrix destructor                            *
 ***************************************************************************/
GMatrix::~GMatrix()
{
  // De-allocate only if memory has indeed been allocated
  if (m_colstart != NULL) delete[] m_colstart;
  if (m_data     != NULL) delete[] m_data;
}


/***************************************************************************
 *                        GMatrix assignment operator                      *
 ***************************************************************************/
GMatrix& GMatrix::operator= (const GMatrix& m)
{
  m_rows     = 0;
  m_cols     = 0;
  m_elements = 0;
  m_data     = NULL;
  m_colstart = NULL;
  m_data     = new double[m.m_elements];
  m_colstart = new unsigned[m.m_cols+1];
  if (m_data == NULL || m_colstart == NULL)
	throw mem_alloc("GMatrix assignment operator", m.m_elements);
  m_rows     = m.m_rows;
  m_cols     = m.m_cols;
  m_elements = m.m_elements;
  for (unsigned i = 0; i <= m_cols; ++i)
    m_colstart[i] = m.m_colstart[i];
  for (unsigned i = 0; i < m_elements; ++i)
    m_data[i] = m.m_data[i];
  return *this;
}


/***************************************************************************
 *                          GMatrix output operator                        *
 ***************************************************************************/
ostream& operator<< (ostream& os, const GMatrix& m)
{
  for (unsigned row = 0; row < m.m_rows; ++row) {
    for (unsigned col = 0; col < m.m_cols; ++col) {
      os << m(row,col);
	  if (col != m.m_cols-1)
	    os << ", ";
	}
	if (row != m.m_rows-1)
	  os << endl;
  }
  return os;
}


/***************************************************************************
 *                       GMatrix * Gvector multiplication                  *
 ***************************************************************************/
GVector GMatrix::operator* (const GVector& v) const
{
  // Raise an exception if the matrix and vector dimensions are not compatible
  if (m_cols != v.m_num)
    throw vec_mat_mismatch("GMatrix*GVector operator", v.m_num, m_rows, m_cols);

  // Perform vector multiplication
  GVector result(m_rows);
  for (unsigned row = 0; row < m_rows; ++row) {
    double sum = 0.0;
	for (unsigned col = 0; col < m_cols; ++col)
	  sum += (*this)(row,col) * v(col);
	result(row) = sum;
  }
  return result;
}


/***************************************************************************
 *                       GMatrix * GMatrix multiplication                  *
 ***************************************************************************/
GMatrix GMatrix::operator* (const GMatrix& m) const
{
  // Raise an exception if the matrix dimensions are not compatible
  if (m_cols != m.m_rows)
    throw GMatrix::dim_mult_mismatch("GMatrix multiplication", m_cols, m.m_rows);

  // Perform matrix multiplication
  GMatrix result(m_rows,m.m_cols);
  for (unsigned row = 0; row < m_rows; ++row) {
    for (unsigned col = 0; col < m.m_cols; ++col) {
	  double sum = 0.0;
	  for (unsigned i = 0; i < m_cols; ++i)
	    sum += (*this)(row,i)*m(i,col);
      result(row,col) = sum;
    }
  }
  return result;
}


/***************************************************************************
 *                          Class exception handlers                       *
 ***************************************************************************/
// Matrix row or column out of range
GMatrix::out_of_range::out_of_range(string origin, unsigned row, unsigned col, 
                                    unsigned rows, unsigned cols)
{
  m_origin = origin;
  ostringstream s_row;
  ostringstream s_col;
  ostringstream s_rows;
  ostringstream s_cols;
  s_row  << row;
  s_col  << col;
  s_rows << rows-1;
  s_cols << cols-1;
  m_message = "Matrix element (" + s_row.str() + "," + s_col.str() +
              ") out of range ([0," + s_rows.str() + "], [0," +
			  s_cols.str() + "])";
}

// Matrix - Vector mismatch
GMatrix::vec_mat_mismatch::vec_mat_mismatch(string origin, unsigned num, 
                                            unsigned rows, unsigned cols)
{
  m_origin = origin;
  ostringstream s_num;
  ostringstream s_rows;
  ostringstream s_cols;
  s_num  << num;
  s_rows << rows;
  s_cols << cols;
  m_message = "Vector dimension [" + s_num.str() + 
              "] is incompatible with matrix size [" + 
			  s_rows.str() + "," + s_cols.str() + "]";
}

// Matrix mismatch (additional etc.)
GMatrix::dim_add_mismatch::dim_add_mismatch(string origin, unsigned rows1, unsigned rows2,
                                            unsigned cols1, unsigned cols2)
{
  m_origin = origin;
  ostringstream s_rows1;
  ostringstream s_rows2;
  ostringstream s_cols1;
  ostringstream s_cols2;
  s_rows1 << rows1;
  s_rows2 << rows2;
  s_cols1 << cols1;
  s_cols2 << cols2;
  m_message = "Size [" + s_rows1.str() + "," + s_cols1.str() +
              "] of first matrix is incompatible with size [" +
			  s_rows2.str() + "," + s_cols2.str() +
			  "] of second matrix";
}

// Matrix mismatch (multiplication)
GMatrix::dim_mult_mismatch::dim_mult_mismatch(string origin, unsigned cols, unsigned rows)
{
  m_origin = origin;
  ostringstream s_rows;
  ostringstream s_cols;
  s_rows << rows;
  s_cols << cols;
  m_message = "Number of columns [" + s_cols.str() + 
              "] in first matrix is incompatible with number of rows [" + 
			  s_rows.str() + "] in second matrix";
}
