/***************************************************************************
 *                 GNodeArray.cpp  -  Array of nodes class                 *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
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
#include <iostream>
#include <math.h>
#include "GException.hpp"
#include "GNodeArray.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                    GNodeArray constructors/destructors                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GNodeArray::GNodeArray()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param array Array to be copied
 ***************************************************************************/
GNodeArray::GNodeArray(const GNodeArray& array)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(array);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GNodeArray::~GNodeArray()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GNodeArray operators                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param array Array to be assigned
 ***************************************************************************/
GNodeArray& GNodeArray::operator= (const GNodeArray& array)
{
    // Execute only if object is not identical
    if (this != &array) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(array);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                         GNodeArray public methods                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set node array
 *
 * @param[in] num Number of nodes
 * @param[in] array Values of node array.
 ***************************************************************************/
void GNodeArray::nodes(const int& num, const double* array)
{
    // Free any existing memory
    free_members();
    
    // Initialise members
    init_members();
    
    // Set number of nodes
    m_nodes = num;
        
    // Allocate memory
    if (m_nodes > 0) {
        m_node = new double[m_nodes];
        m_step = new double[m_nodes];
    }
    
    // Copy node data
    for (int i = 0; i < m_nodes; ++i)
        m_node[i] = array[i];

    // Setup distance array
    for (int i = 0; i < m_nodes-1; ++i)
        m_step[i] = m_node[i+1] - m_node[i];

    // Evaluate linear slope and offset
    m_linear_slope  = double(m_nodes) / (m_node[m_nodes-1] - m_node[0]);
    m_linear_offset = -m_linear_slope * m_node[0];
    
    // Check if nodes form a linear array
    m_is_linear = 1;
    for (int i = 0; i < m_nodes-1; ++i) {
        double eps = m_linear_slope * m_node[i] + m_linear_offset - double(i);
        if (fabs(eps) > 1.0e-6) {
            cout << "WARNING: Node " << i 
                 << "indexing is invalid (eps=" << eps << ")" << endl;
            m_is_linear = 0;
        }
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set node array from vector
 *
 * @param[in] vector Vector from which node array will be built.
 ***************************************************************************/
void GNodeArray::nodes(const GVector& vector)
{
    // Free any existing memory
    free_members();
    
    // Initialise members
    init_members();
    
    // Set number of nodes
    m_nodes = vector.size();
        
    // Allocate memory
    if (m_nodes > 0) {
        m_node = new double[m_nodes];
        m_step = new double[m_nodes];
    }
    
    // Copy node data
    for (int i = 0; i < m_nodes; ++i)
        m_node[i] = vector(i);

    // Setup distance array
    for (int i = 0; i < m_nodes-1; ++i)
        m_step[i] = m_node[i+1] - m_node[i];

    // Evaluate linear slope and offset
    m_linear_slope  = double(m_nodes) / (m_node[m_nodes-1] - m_node[0]);
    m_linear_offset = -m_linear_slope * m_node[0];
    
    // Check if nodes form a linear array
    m_is_linear = 1;
    for (int i = 0; i < m_nodes-1; ++i) {
        double eps = m_linear_slope * m_node[i] + m_linear_offset - double(i);
        if (fabs(eps) > 1.0e-6) {
            cout << "WARNING: Node " << i 
                 << "indexing is invalid (eps=" << eps << ")" << endl;
            m_is_linear = 0;
        }
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set indices and weighting factors for interpolation
 *
 * @param[in] value Value for which the interpolation should be done.
 ***************************************************************************/
void GNodeArray::set_value(const double& value)
{
    // If array is linear
    if (m_is_linear) {
    
        // Set left index
        m_inx_left = m_linear_slope * value + m_linear_offset;
        
        // Keep index in valid range
        if (m_inx_left < 0)               m_inx_left = 0;
        else if (m_inx_left >= m_nodes-1) m_inx_left = m_nodes - 2;
        
        // Set right index
        m_inx_right = m_inx_left + 1;

        // Set weighting factors
        m_wgt_right = (value - m_node[m_inx_left]) / m_step[m_inx_left];
        m_wgt_left  = 1.0 - m_wgt_right;

    } // endif: array is linear
    
    // ... otherwise search the relevant bin
    {
        // TBD
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GNodeArray private methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GNodeArray::init_members(void)
{
    // Initialise members
    m_nodes         = 0;
    m_node          = NULL;
    m_step          = NULL;
    m_is_linear     = 0;
    m_linear_slope  = 0.0;
    m_linear_offset = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] array Array to be copied
 ***************************************************************************/
void GNodeArray::copy_members(const GNodeArray& array)
{
    // Copy number of bins
    m_nodes         = array.m_nodes;
    m_is_linear     = array.m_is_linear;
    m_linear_slope  = array.m_linear_slope;
    m_linear_offset = array.m_linear_offset;
    
    // Copy nodes
    if (m_nodes > 0) {
        m_node = new double[m_nodes];
        m_step = new double[m_nodes];
        memcpy(m_node, array.m_node, m_nodes*sizeof(double));
        memcpy(m_step, array.m_step, m_nodes*sizeof(double));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GNodeArray::free_members(void)
{
    // Free memory
    if (m_node != NULL) delete [] m_node;
    if (m_step != NULL) delete [] m_step;

    // Signal that memory is free
    m_node = NULL;
    m_step = NULL;

    // Return
    return;
}
