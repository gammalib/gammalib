/***************************************************************************
 *                     GFitPar.cpp  -  fit parameter class                 *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2007 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GFitPar.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Method name definitions ____________________________________________ */
//#define G_OP_MUL_VEC   "GSparseMatrix::operator* (const GVector&) const"

/* __ Macros _____________________________________________________________ */
//#define G_MIN(a,b) (((a) < (b)) ? (a) : (b))

/* __ Coding definitions _________________________________________________ */
//#define G_USE_MEMCPY

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_SPARSE_PENDING                    // Analyse pending values



/*==========================================================================
 =                                                                         =
 =                      GFitPar constructors/destructors                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                               Constructor                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitPar::GFitPar(void)
{
  // Initialise private members for clean destruction
  init_members();

  // Construct fit parameter
  constructor();
  
  // Return
  return;
}


/***************************************************************************
 *                             Copy constructor                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitPar::GFitPar(const GFitPar& p)
{ 
  // Initialise private members for clean destruction
  init_members();

  // Copy members
  copy_members(p);

  // Return
  return;
}


/***************************************************************************
 *                                Destructor                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitPar::~GFitPar()
{
  // Free members
  free_members();

  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                            GFitPar operators                            =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitPar& GFitPar::operator= (const GFitPar& p)
{ 
  // Execute only if object is not identical
  if (this != &p) {

    // Free members
    free_members();

    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(p);

  } // endif: object was not identical
  
  // Return
  return *this;
}


/*==========================================================================
 =                                                                         =
 =                             GFitPar methods                             =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                          Assign value to parameter                      *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitPar& GFitPar::operator= (const double& v)
{
  // Assign value
  m_value = v;
  
  // If there is a minimum boundary then make sure that the value is above 
  // the boundary
  if (m_has_min && m_value < m_value_min)
    m_value = m_value_min;

  // If there is a maximum boundary then make sure that the value is below 
  // the boundary
  if (m_has_max && m_value > m_value_max)
    m_value = m_value_max;
	
  // Return
  return *this;
}


/***************************************************************************
 *                        Set minimum parameter boundary                   *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitPar::min_set(double& min)
{
  // Assign minimum boundary and set flag
  m_value_min = min;
  m_has_min   = 1;
  
  // If there is a maximum boundary then make sure that the minimum boundary
  // is below the maximum boundary
  if (m_has_max && m_value_min > m_value_max)
    m_value_min = m_value_max;

  // Make sure that the parameter falls within the boundaries
  if (m_has_min && m_value < m_value_min)
    m_value = m_value_min;
  if (m_has_max && m_value > m_value_max)
    m_value = m_value_max;
  
  // Return
  return;
}


/***************************************************************************
 *                        Set maximum parameter boundary                   *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitPar::max_set(double& max)
{
  // Assign maximum boundary and set flag
  m_value_max = max;
  m_has_max   = 1;
  
  // If there is a minimum boundary then make sure that the minimum boundary
  // is below the maximum boundary
  if (m_has_min && m_value_max < m_value_min)
    m_value_max = m_value_min;

  // Make sure that the parameter falls within the boundaries
  if (m_has_min && m_value < m_value_min)
    m_value = m_value_min;
  if (m_has_max && m_value > m_value_max)
    m_value = m_value_max;
  
  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                        GFitPar private functions                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Constructor method                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitPar::constructor(void)
{

  // Return
  return;
}



/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitPar::init_members(void)
{
  // Initialise sparse matrix members
  m_value     = 0.0;
  m_error     = 0.0;
  m_value_min = 0.0;
  m_value_max = 0.0;
  m_has_error = 0;                  // Parameter has no associated error
  m_has_min   = 0;                  // Parameter has no maximum boundary
  m_has_max   = 0;                  // Parameter has no minimum boundary
  m_free      = 1;                  // Parameter is free
  m_name      = "";                 // Parameter has blank name
  
  // Return
  return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitPar::copy_members(const GFitPar& p)
{
  // Copy members
  m_value     = p.m_value;
  m_error     = p.m_error;
  m_value_min = p.m_value_min;
  m_value_max = p.m_value_max;
  m_has_error = p.m_has_error;
  m_has_min   = p.m_has_min;
  m_has_max   = p.m_has_max;
  m_free      = p.m_free;
  m_name      = p.m_name;

  // Return
  return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitPar::free_members(void)
{
  
  // Return
  return;
}


/*==========================================================================
 =                                                                         =
 =                              GFitPar friends                            =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
ostream& operator<< (ostream& os, const GFitPar& p)
{
  // Put header in stream
  os << "=== GFitPar ===" << endl;
  os << " Value .....................: " << p.m_value;
  if (p.m_free == 0)
    os << " (fixed)" << endl;
  else
    os << " (free)" << endl;
  if (p.m_has_error)
    os << " Error .....................: " << p.m_error << endl;
  if (p.m_has_min)
    os << " Minimum boundary ..........: " << p.m_value_min << endl;
  if (p.m_has_max)
    os << " Maximum boundary ..........: " << p.m_value_max << endl;

  // Return output stream
  return os;
}


/*==========================================================================
 =                                                                         =
 =                      Other functions used by GFitPar                    =
 =                                                                         =
 ==========================================================================*/
