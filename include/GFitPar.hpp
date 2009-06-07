/***************************************************************************
 *                   GFitPar.hpp  -  fit parameter class                   *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2007 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GFITPAR_HPP
#define GFITPAR_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include <iostream>                           // std::ostream


/***************************************************************************
 *                         GFitPar class definition                        *
 ***************************************************************************/
class GFitPar {
  // Friend classes

  // Operator friends

  // I/O friends
  friend std::ostream& operator<< (std::ostream& os, const GFitPar& p);

  // Friend functions

public:
  // Constructors and destructors
  GFitPar(void);
  GFitPar(const GFitPar& p);
 ~GFitPar();
 
  // Conversion operators
  operator double() const;                // Get parameter value

  // Fit parameter operators
  GFitPar& operator= (const GFitPar& p);  // Copy parameter
  GFitPar& operator= (const double& v);   // Assign parameter value

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
 *                               Inline members                            *
 ***************************************************************************/
// Conversion operator: double value = parameter
inline
GFitPar::operator double() const
{
  return m_value;
}

// Return associated parameter error
inline
double GFitPar::error() const
{
  return m_error;
}

// Return minimum parameter boundary
inline
double GFitPar::min() const
{
  return m_value_min;
}

// Return maximum parameter boundary
inline
double GFitPar::max() const
{
  return m_value_max;
}

// Remove minimum parameter boundary
inline
void GFitPar::min_remove(void)
{
  m_value_min = 0.0;
  m_has_min   = 0;
  return;
}

// Remove maximum parameter boundary
inline
void GFitPar::max_remove(void)
{
  m_value_max = 0.0;
  m_has_max   = 0;
  return;
}

// Free parameter
inline
void GFitPar::free(void)
{
  m_free = 1;
  return;
}

// Fix parameter
inline
void GFitPar::fix(void)
{
  m_free = 0;
  return;
}


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/

#endif /* GFITPAR_HPP */
