/***************************************************************************
 *                GFitPar.i  -  Fit parameter class SWIG file              *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitPar.hpp"
%}


/***************************************************************************
 *                         GFitPar class definition                        *
 ***************************************************************************/
class GFitPar {
public:
  // Constructors and destructors
  GFitPar(void);
  GFitPar(const GFitPar& p);
 ~GFitPar();

  // Fit parameter methods
  double error(void) const;               // Get error value
  double min(void) const;                 // Get minimum value boundary
  double max(void) const;                 // Get maximum value boundary
  void   min_set(double& min);            // Set minimum value boundary
  void   max_set(double& max);            // Set maximum value boundary
  void   min_remove(void);                // Remove minimum value boundary
  void   max_remove(void);                // Remove maximum value boundary
  void   free(void);                      // Free parameter
  void   fix(void);                       // Fix parameter
  
private:
  // Private methods
  void constructor(void);
  void init_members(void);
  void copy_members(const GFitPar& p);
  void free_members(void);

  // Private data members
  double       m_value;                   // Parameter value
  double       m_error;                   // Parameter error
  double       m_value_min;               // Minimum parameter value
  double       m_value_max;               // Maximum parameter value
  int          m_has_error;               // Parameter has error
  int          m_has_min;                 // Parameter has a minimum boundary
  int          m_has_max;                 // Parameter has a maximum boundary
  int          m_free;                    // Parameter is free
  std::string  m_name;                    // Parameter name
};


/***************************************************************************
 *                      GFitPar class SWIG extension                       *
 ***************************************************************************/
%extend GFitPar {
  char *__str__() {
    std::ostringstream buffer;
	buffer << *self;
	std::string str = buffer.str();
	char* ptr = (char*)str.c_str();
	return ptr;
  }
  double value(void) const {
    return double(*self);
  }
  void set(const double val) {
    (*self) = val;
  }
  void min_set(const double val) {
    double arg = val;
    (*self).min_set(arg);
  }
  void max_set(const double val) {
    double arg = val;
    (*self).max_set(arg);
  }
  GFitPar copy() {
    return (*self);
  }
}
