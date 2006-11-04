/***************************************************************************
 *                   GException.hpp  -  exception handler                  *
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

#ifndef GEXCEPTION_HPP
#define GEXCEPTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                  Exception handler base class definition                *
 ***************************************************************************/
class GExceptionHandler : public exception {
public:
  GExceptionHandler() { }
 ~GExceptionHandler() throw() { }
  virtual const char* what() const throw();
protected:
  string m_origin;
  string m_message;
};


/***************************************************************************
 *                         Exception class definition                      *
 ***************************************************************************/
class GException : public GExceptionHandler {
public:

  // Memory allocation exception class
  class mem_alloc : public GExceptionHandler {
  public:
    mem_alloc(string origin, unsigned num);
  };

  // Empty object exception class
  class empty : public GExceptionHandler {
  public:
    empty(string origin);
  };

  // Index out of range
  class out_of_range : public GExceptionHandler {
  public:
    out_of_range(string origin, int row, int col, int rows, int cols);
  };
  
  // Vector - Matrix mismatch
  class matrix_vector_mismatch : public GExceptionHandler {
  public:
    matrix_vector_mismatch(string origin, int num, int rows, int cols);
  };

  // Matrix dimensions mismatch
  class matrix_mismatch : public GExceptionHandler {
  public:
    matrix_mismatch(string origin, int rows1, int cols1, int rows2, int cols2);
  };

  // Matrix not rectangular
  class matrix_not_rectangular : public GExceptionHandler {
  public:
	matrix_not_rectangular(string origin, int rows, int cols);
  };

  // Matrix not positive definite
  class matrix_not_pos_definite : public GExceptionHandler {
  public:
    matrix_not_pos_definite(string origin, int row, double sum);
  };

  // Matrix not symmetric
  class matrix_not_symmetric : public GExceptionHandler {
  public:
    matrix_not_symmetric(string origin, int cols, int rows);
  };

  // Matrix not factorised
  class matrix_not_factorised : public GExceptionHandler {
  public:
    matrix_not_factorised(string origin, string type);
  };

  // All matrix elements are zero
  class matrix_zero : public GExceptionHandler  {
  public:
    matrix_zero(string origin);
  };

  // Invalid ordering scheme
  class invalid_order : public GExceptionHandler {
  public:
    invalid_order(string origin, int order, int min_order, int max_order);
  };

};

#endif /* GEXCEPTION_HPP */
