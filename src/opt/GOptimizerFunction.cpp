/***************************************************************************
 *  GOptimizerFunction.cpp  -  Abstract base class for optimizer function  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/**
 * @file GOptimizerFunction.cpp
 * @brief GOptimizerFunction abstract base class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GOptimizerFunction.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =               GOptimizerFunction constructors/destructors               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GOptimizerFunction::GOptimizerFunction()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] fct Optimizer function from which the instance should be built.
 ***************************************************************************/
GOptimizerFunction::GOptimizerFunction(const GOptimizerFunction& fct)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(fct);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GOptimizerFunction::~GOptimizerFunction()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GOptimizerFunction operators                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] fct Optimizer function to be assigned.
 ***************************************************************************/
GOptimizerFunction& GOptimizerFunction::operator= (const GOptimizerFunction& fct)
{
    // Execute only if object is not identical
    if (this != &fct) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(fct);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                     GOptimizerFunction item iterator                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Iterator Constructor
 ***************************************************************************/
GOptimizerFunction::iterator::iterator()
{
    // Initialise iterator
    m_base  = NULL;
    m_index = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Iterator Constructor
 ***************************************************************************/
GOptimizerFunction::iterator::iterator(GOptimizerFunction *fct)
{
    // Initialise iterator
    m_base  = fct;
    m_index = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Iterator prefix operator
 ***************************************************************************/
GOptimizerFunction::iterator& GOptimizerFunction::iterator::operator++(void)
{
    // Increment index
    m_index++;
    
    // Move to next item. If end is reached then set index to -1
    if (m_base->next_item())
        m_index = -1;
    
    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Iterator postfix operator
 ***************************************************************************/
GOptimizerFunction::iterator GOptimizerFunction::iterator::operator++(int junk)
{
    // Save actual iterator
    GOptimizerFunction::iterator actual = *this;
    
    // Increment using prefix iterator
    ++(*this);
    
    // Return actual iterator
    return actual;
}


/***********************************************************************//**
 * @brief Set iterator on first elements
 ***************************************************************************/
GOptimizerFunction::iterator GOptimizerFunction::begin(void)
{
    // Allocate iterator object
    GOptimizerFunction::iterator iter(this);
    
    // Initialise index
    iter.m_index = 0;
    
    // Set iterator to first data/model item
    if (iter.m_base != NULL)
        iter.m_base->first_item();

    // Return
    return iter;
}


/***********************************************************************//**
 * @brief Set iterator to end
 ***************************************************************************/
GOptimizerFunction::iterator GOptimizerFunction::end(void)
{
    // Allocate iterator object
    GOptimizerFunction::iterator iter(this);

    // End is flagged by an index of -1
    iter.m_index = -1;

    // Return
    return iter;
}


/*==========================================================================
 =                                                                         =
 =                    GOptimizerFunction public methods                    =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                    GOptimizerFunction private methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GOptimizerFunction::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] fct GOptimizerFunction members which should be copied.
 ***************************************************************************/
void GOptimizerFunction::copy_members(const GOptimizerFunction& fct)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GOptimizerFunction::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GOptimizerFunction friends                        =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                 Other functions used by GOptimizerFunction              =
 =                                                                         =
 ==========================================================================*/
