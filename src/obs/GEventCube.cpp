/***************************************************************************
 *          GEventCube.cpp  -  Abstract event cube container class         *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
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
 * @file GEventCube.cpp
 * @brief GEventCube container class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GEventCube.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COPY_MEMBERS          "GEventCube::copy_members(const GEventCube&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                    GEventCube constructors/destructors                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GEventCube::GEventCube() : GEvents()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Event cube from which the instance should be built.
 ***************************************************************************/
GEventCube::GEventCube(const GEventCube& cube) : GEvents(cube)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GEventCube::~GEventCube()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GEventCube operators                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] cube Event cube to be assigned.
 ***************************************************************************/
GEventCube& GEventCube::operator= (const GEventCube& cube)
{
    // Execute only if object is not identical
    if (this != &cube) {

        // Copy base class members
        this->GEvents::operator=(cube);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(cube);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                       GEventCube public methods                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return number of elements in cube
 ***************************************************************************/
int GEventCube::elements(void) const
{
    // Return
    return m_elements;
}


/***********************************************************************//**
 * @brief Return event cube dimension
 ***************************************************************************/
int GEventCube::dim(void) const
{
    // Return
    return m_dim;
}


/***********************************************************************//**
 * @brief Return event cube axis dimension
 *
 * @param[in] axis Cube axis number (starting from 0).
 ***************************************************************************/
int GEventCube::naxis(int axis) const
{
    // Throw error if index is out of range
    #if defined(G_RANGE_CHECK)
    if (m_naxis == NULL || axis < 0 || axis >= m_dim)
        throw GException::out_of_range(G_NAXIS, axis, 0, m_dim-1);
    #endif
        
    // Return axis dimension
    return m_naxis[axis];
}


/*==========================================================================
 =                                                                         =
 =                       GEventCube private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GEventCube::init_members(void)
{
    // Initialise members
    m_elements = 0;
    m_dim      = 0;
    m_naxis    = NULL;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube GEventCube members to be copied.
 ***************************************************************************/
void GEventCube::copy_members(const GEventCube& cube)
{
    // Copy attributes
    m_elements = cube.m_elements;
    m_dim      = cube.m_dim;

    // If the cube is not empty then copy it
    if (m_dim > 0 && cube.m_naxis != NULL) {
    
        // Allocate memory
        m_naxis = new int[m_dim];
        if (m_naxis == NULL)
            throw GException::mem_alloc(G_COPY_MEMBERS, m_dim);

        // Copy axis dimensions
        for (int i = 0; i < m_dim; ++i)
            m_naxis[i] = cube.m_naxis[i];

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GEventCube::free_members(void)
{
    // Free memory
    if (m_naxis != NULL) delete [] m_naxis;

    // Signal free pointers
    m_naxis = NULL;

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GEventCube friends                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream into which the event cube will be dumped
 * @param[in] cube Event cube to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GEventCube& cube)
{
    // Put event cube in stream
    os << "=== GEventCube ===" << std::endl;
    if (cube.m_naxis != NULL) {
        os << " Cube is undefined" << std::endl;
    }
    else {
        os << " Number of cube elements ...: " << cube.elements() << std::endl;
        os << " Cube dimension ............: " << cube.dim() << std::endl;
        for (int axis = 0; axis < cube.m_dim; ++axis) {
            os << " Axis " << axis+1 << " dimension ..........: " 
               << cube.naxis(axis) << std::endl;
        }
    }
        
    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GEventList                    =
 =                                                                         =
 ==========================================================================*/
