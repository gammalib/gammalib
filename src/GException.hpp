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
using namespace std;

/* __ Enumerators ________________________________________________________ */

/* __ Structures _________________________________________________________ */

/* __ Prototypes _________________________________________________________ */


/***************************************************************************
 *                  Exception handler base class definition                *
 ***************************************************************************/
class GException : public exception {
public:
  GException() { }
 ~GException() throw() { }
  virtual const char* what() const throw();
protected:
  string m_origin;
  string m_message;
};


/***************************************************************************
 *                    Specific exception class definitions                 *
 ***************************************************************************/
class mem_alloc : public GException {
public:
  mem_alloc(string origin, unsigned num);
};
class empty : public GException {
public:
  empty(string origin);
};

#endif /* GEXCEPTION_HPP */
