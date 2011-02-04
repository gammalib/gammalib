/***************************************************************************
 *                 GNodeArray.cpp  -  Array of nodes class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GNodeArray.cpp
 * @brief GNodeArray class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <iostream>
#include "GException.hpp"
#include "GNodeArray.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                               "GNodeArray::operator() (int)"
#define G_INTERPOLATE "GNodeArray::interpolate(std::vector<double>&, double&"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GNodeArray::GNodeArray(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] array Instance from which object should be constructed.
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
GNodeArray::~GNodeArray(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] array Instance that should be assigned.
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


/***********************************************************************//**
 * @brief Node access operator
 *
 * @param[in] index Node index (0,1,...).
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 ***************************************************************************/
double& GNodeArray::operator() (int index)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_ACCESS, index, size()-1);
    #endif

    // Return node
    return m_node[index];
}


/***********************************************************************//**
 * @brief Node access operator (const version)
 *
 * @param[in] index Node index (0,1,...).
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 ***************************************************************************/
const double& GNodeArray::operator() (int index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size())
        throw GException::out_of_range(G_ACCESS, index, size()-1);
    #endif

    // Return node
    return m_node[index];
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GNodeArray::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GNodeArray* GNodeArray::clone(void) const
{
    return new GNodeArray(*this);
}


/***********************************************************************//**
 * @brief Set node array
 *
 * @param[in] num Number of nodes
 * @param[in] array Node values \f$x_i\f$.
 ***************************************************************************/
void GNodeArray::nodes(const int& num, const double* array)
{
    // Free any existing memory
    free_members();
    
    // Initialise members
    init_members();
    
    // Set node values
    for (int i = 0; i < num; ++i)
        m_node.push_back(array[i]);

    // Setup node distances and linear array handling
    setup();

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
    
    // Set node values
    for (int i = 0; i < vector.size(); ++i)
        m_node.push_back(vector(i));

    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set node array from vector
 *
 * @param[in] vector Vector from which node array will be built.
 ***************************************************************************/
void GNodeArray::nodes(const std::vector<double>& vector)
{
    // Free any existing memory
    free_members();
    
    // Initialise members
    init_members();
    
    // Set node values
    m_node = vector;

    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append one node to array.
 *
 * @param[in] node Node to be appended to array.
 ***************************************************************************/
void GNodeArray::append(const double& node)
{
    // Add node
    m_node.push_back(node);
    
    // Setup node distances and linear array handling
    setup();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Interpolate value.
 *
 * @param[in] value Value \f$x\f$ at which interpolation should be done.
 * @param[in] vector Vector \f$y_i\f$ that should be interpolated.
 *
 * @exception GException::not_enough_nodes
 *            Not enough nodes for interpolation in node array.
 * @exception GException::vector_mismatch
 *            Size of node vector does not match the size of vector argument.
 *
 * This method performs a linear interpolation of values \f$y_i\f$. The
 * corresponding values \f$x_i\f$ are stored in the node array.
 ***************************************************************************/
double GNodeArray::interpolate(const double& value,
                               const std::vector<double>& vector)
{
    // Throw exception if there are not enough nodes
    if (m_node.size() < 2)
        throw GException::not_enough_nodes(G_INTERPOLATE, m_node.size());

    // Throw exception if vectors have not the same size
    if (m_node.size() != vector.size())
        throw GException::vector_mismatch(G_INTERPOLATE, m_node.size(),
                                          vector.size());
    
    // Set interpolation value
    set_value(value);

    // Interpolate
    double y = vector[inx_left()]  * wgt_left() +
               vector[inx_right()] * wgt_right();

    // Return
    return y;
}


/***********************************************************************//**
 * @brief Set indices and weighting factors for interpolation
 *
 * @param[in] value Value for which the interpolation should be done.
 *
 * Set the indices that bound the specified value and the corresponding
 * weighting factors for linear interpolation. If the array has a linear
 * form (i.e. the nodes are equidistant), an analytic formula is used to
 * determine the boundary indices. If the nodes are not equidistant the
 * boundary indices are searched by bisection.
 ***************************************************************************/
void GNodeArray::set_value(const double& value)
{
    // Get number of nodes
    int nodes = m_node.size();
    
    // Continue only if we have at least 2 nodes
    if (nodes > 1) {
    
        // If array is linear then get left index from analytic formula
        if (m_is_linear) {

            // Set left index
            m_inx_left = int(m_linear_slope * value + m_linear_offset);

            // Keep index in valid range
            if (m_inx_left < 0)             m_inx_left = 0;
            else if (m_inx_left >= nodes-1) m_inx_left = nodes - 2;

        } // endif: array is linear

        // ... otherwise search the relevant indices by bisection
        else {
            // Set left index if value is before first node
            if (value < m_node[0])
                m_inx_left = 0;

            // Set left index if value is after last node
            else if (value >  m_node[nodes-1])
                m_inx_left = nodes - 2;

            // Set left index by bisection
            else {
                int low  = 0;
                int high = nodes - 1;
                while ((high - low) > 1) {
                    int mid = (low+high) / 2;
                    if (m_node[mid] > value)
                        high = mid;
                    else if (m_node[mid] <= value)
                        low = mid;
                }
                m_inx_left = low;
            } // endelse: did bisection
        }

        // Set right index
        m_inx_right = m_inx_left + 1;

        // Set weighting factors
        m_wgt_right = (value - m_node[m_inx_left]) / m_step[m_inx_left];
        m_wgt_left  = 1.0 - m_wgt_right;
    
    } // endif: there were at least 2 nodes

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
    m_node.clear();
    m_step.clear();
    m_is_linear     = 0;
    m_linear_slope  = 0.0;
    m_linear_offset = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] array Node array to be copied.
 ***************************************************************************/
void GNodeArray::copy_members(const GNodeArray& array)
{
    // Copy number of bins
    m_node          = array.m_node;
    m_step          = array.m_step;
    m_is_linear     = array.m_is_linear;
    m_linear_slope  = array.m_linear_slope;
    m_linear_offset = array.m_linear_offset;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GNodeArray::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute distance array and linear slope/offset
 ***************************************************************************/
void GNodeArray::setup(void)
{
    // Reset distance vector
    m_step.clear();
    
    // Get number of nodes
    int nodes = m_node.size();
    
    // Continue only if we have nodes at least 2 nodes
    if (nodes > 1) {
    
        // Setup distance array between subsequent nodes
        for (int i = 0; i < nodes-1; ++i)
            m_step.push_back(m_node[i+1] - m_node[i]);

        // Evaluate linear slope and offset
        m_linear_slope  = double(nodes-1) / (m_node[nodes-1] - m_node[0]);
        m_linear_offset = -m_linear_slope * m_node[0];
    
        // Check if nodes form a linear array
        m_is_linear = true;
        for (int i = 0; i < nodes-1; ++i) {
            double eps = m_linear_slope * m_node[i] + m_linear_offset - double(i);
            if (std::abs(eps) > 1.0e-6) {
                m_is_linear = false;
                break;
            }
        }
        
    } // endif: there were at least 2 nodes

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put object in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] array Node array to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os,const GNodeArray& array)
{
    // Put object in stream
    os << "=== GNodeArray ===" << std::endl;
    os << " Number of nodes in array ..: " << array.m_node.size() << std::endl;
    if (array.m_is_linear) {
        os << " Array type ................: linear" << std::endl;
        os << " Linear slope ..............: " << array.m_linear_slope
           << std::endl;
        os << " Linear offset .............: " << array.m_linear_offset
           << std::endl;
    }
    else
        os << " Array type ................: nonlinear" << std::endl;
    os << " Indices and weights .......: (" 
       << array.m_inx_left << "," 
       << array.m_inx_right << ")=("
       << array.m_wgt_left << ","
       << array.m_wgt_right << ")" << std::endl;
    for (int i = 0; i < array.m_node.size(); ++i)
        os << " " << array.m_node.at(i);

    // Return output stream
    return os;
}
