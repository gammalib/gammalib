/***************************************************************************
 *                   GException.cpp  -  exception handler                  *
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
#include "GException.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                  Exception handler base class definition                *
 ***************************************************************************/
const char* GExceptionHandler::what() const throw()
{
  string message = "*** ERROR in " + m_origin + ": " + m_message;
  return message.c_str();
}


/***************************************************************************
 *                        Memory allocation exception                      *
 ***************************************************************************/
GException::mem_alloc::mem_alloc(string origin, unsigned num)
{
  ostringstream s_num;
  s_num << num;
  m_origin  = origin;
  m_message = "Memory allocation error (" + s_num.str() + " elements)";
}


/***************************************************************************
 *                           Empty object exception                        *
 ***************************************************************************/
GException::empty::empty(string origin)
{
  m_origin  = origin;
  m_message = "Zero-size allocation";
}


/***************************************************************************
 *                          Vector index out of range                      *
 ***************************************************************************/
GException::out_of_range::out_of_range(string origin, int inx, int elements)
{
  m_origin = origin;
  if (elements > 0) {
    ostringstream s_inx;
    ostringstream s_elements;
    s_inx      << inx;
    s_elements << elements-1;
    m_message = "Vector index (" + s_inx.str() + ") out of range [0," + 
                s_elements.str() + "]";
  }
  else {
    m_message = "Empty vector";
  }
}



/***************************************************************************
 *                      Matrix row or column out of range                  *
 ***************************************************************************/
GException::out_of_range::out_of_range(string origin, int row, int col, int rows, int cols)
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


/***************************************************************************
 *                          Vector dimensions differ                       *
 ***************************************************************************/
GException::vector_mismatch::vector_mismatch(string origin, int size1, int size2)
{
  m_origin = origin;
  ostringstream s_size1;
  ostringstream s_size2;
  s_size1 << size1;
  s_size2 << size2;
  m_message = "Vector dimensions differ (" + s_size1.str() + " <-> " + 
              s_size2.str() + ")";
}


/***************************************************************************
 *                   Invalid vector dimension for cross product            *
 ***************************************************************************/
GException::vector_bad_cross_dim::vector_bad_cross_dim(int elements)
{
  ostringstream s_elements;
  s_elements << elements;
  m_message = "Vector cross product only defined for 3 dimensions but vector size is " + 
              s_elements.str(); 
}


/***************************************************************************
 *                     Mismatch between matrix and vector                  *
 ***************************************************************************/
GException::matrix_vector_mismatch::matrix_vector_mismatch(string origin, int num, int rows, int cols)
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


/***************************************************************************
 *                       Matrix mismatch in operation                      *
 ***************************************************************************/
GException::matrix_mismatch::matrix_mismatch(string origin, int rows1, int cols1, int rows2, int cols2)
{
  m_origin = origin;
  ostringstream s_rows1;
  ostringstream s_rows2;
  ostringstream s_cols1;
  ostringstream s_cols2;
  s_rows1 << rows1;
  s_cols1 << cols1;
  s_rows2 << rows2;
  s_cols2 << cols2;
  m_message = "Matrix mismatch: M1(" + s_rows1.str() + "," + s_cols1.str() +
              ") incompatible with M2(" + s_rows2.str() + "," + 
			  s_cols2.str() + ")";
}


/***************************************************************************
 *                           Matrix not rectangular                        *
 ***************************************************************************/
GException::matrix_not_rectangular::matrix_not_rectangular(string origin, int rows, int cols)
{
  m_origin = origin;
  ostringstream s_rows;
  ostringstream s_cols;
  s_rows << rows;
  s_cols << cols;
  m_message = "Matrix is not rectangular [" + s_rows.str() + "," + s_cols.str() +
              "]";
}


/***************************************************************************
 *                     Matrix is not positive definite                     *
 ***************************************************************************/
GException::matrix_not_pos_definite::matrix_not_pos_definite(string origin, int row, double sum)
{
  m_origin = origin;
  ostringstream s_rows;
  ostringstream s_sum;
  s_rows << row;
  s_sum << scientific << sum;
  m_message = "Matrix is not positive definite (sum " + s_sum.str() + 
              " occured in row/column " + s_rows.str() + ")";
}


/***************************************************************************
 *                         Matrix is not symmetric                         *
 ***************************************************************************/
GException::matrix_not_symmetric::matrix_not_symmetric(string origin, int cols, int rows)
{
  m_origin = origin;
  ostringstream s_rows;
  ostringstream s_cols;
  s_rows << rows;
  s_cols << cols;
  m_message = "Matrix is not symmetric [" + s_rows.str() + "," +
			  s_cols.str() + "]";
}


/***************************************************************************
 *                       Matrix has not been factorised                    *
 ***************************************************************************/
GException::matrix_not_factorised::matrix_not_factorised(string origin, string type)
{
  m_origin  = origin;
  m_message = "Matrix has not been factorised using " + type;
}


/***************************************************************************
 *                       All matrix elements are zero                      *
 ***************************************************************************/
GException::matrix_zero::matrix_zero(string origin)
{
  m_origin = origin;
  m_message = "All matrix elements are zero";
}


/***************************************************************************
 *                     Invalid ordering scheme requested                   *
 ***************************************************************************/
GException::invalid_order::invalid_order(string origin, int order, int min_order, int max_order)
{
  m_origin = origin;
  ostringstream s_order;
  ostringstream s_min_order;
  ostringstream s_max_order;
  s_order     << order;
  s_min_order << min_order;
  s_max_order << max_order;
  m_message = "Invalid ordering type " + s_order.str() + 
              "requested; must be comprised in [" + s_min_order.str() +
			  "," + s_max_order.str() + "]";
}
