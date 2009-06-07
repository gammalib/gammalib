/***************************************************************************
 *                 GNodeArray.hpp  -  Array of nodes class                 *
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
/**
 * @file GNodeArray.hpp
 * @brief GNodeArray class interface definition.
 * @author J. Knodlseder
 */

#ifndef GNODEARRAY_HPP
#define GNODEARRAY_HPP

/* __ Includes ___________________________________________________________ */
#include "GVector.hpp"


/***********************************************************************//**
 * @class GNodeArray
 *
 * @brief Interface for the node array class.
 *
 * The node array class collects a number of nodes that may be used to 
 * describe a functional relation. This class may be used to perform a linear
 * interpolation between these nodes. Using the GNodeArray::set_value method, 
 * the indices of the nearest nodes are determined and the weighting factors
 * for linear interpolation are computed. These indices and weighting factors
 * can then be accessed via the methods 
 * GNodeArray::inx_left,
 * GNodeArray::inx_right,
 * GNodeArray::wgt_left, and
 * GNodeArray::wgt_right.
 *
 * If the nodes are equally spaced, interpolation is more rapid.
 ***************************************************************************/
class GNodeArray {

public:
    // Constructors and destructors
    GNodeArray();
    GNodeArray(const GNodeArray& array);
    ~GNodeArray();

    // Operators
    GNodeArray& operator= (const GNodeArray & array);

    // Methods
    void   nodes(const int& num, const double* array);
    void   nodes(const GVector& vector);
    void   set_value(const double& value);
    int    inx_left(void) { return m_inx_left; }
    int    inx_right(void) { return m_inx_right; }
    double wgt_left(void) { return m_wgt_left; }
    double wgt_right(void) { return m_wgt_right; }

private:
    // Methods
    void init_members(void);
    void copy_members(const GNodeArray& array);
    void free_members(void);
    void setup(void);
    
    // Data
    int     m_nodes;          //!< Number of nodes
    double* m_node;           //!< Node values
    double* m_step;           //!< Distance to next node
    int     m_is_linear;      //!< Nodes form a linear array
    double  m_linear_slope;   //!< Slope for linear array
    double  m_linear_offset;  //!< Offset for linear array
    int     m_inx_left;       //!< Index of left node for linear interpolation
    int     m_inx_right;      //!< Index of right node for linear interpolation
    double  m_wgt_left;       //!< Weight for left node for linear interpolation
    double  m_wgt_right;      //!< Weight for right node for linear interpolation
};

#endif /* GNODEARRAY_HPP */
