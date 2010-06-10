/***************************************************************************
 *                 GNodeArray.hpp  -  Array of nodes class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
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
#include <vector>
#include <iostream>
#include "GVector.hpp"


/***********************************************************************//**
 * @class GNodeArray
 *
 * @brief Interface for the node array class.
 *
 * The node array class collects nodes \f$x_i\f$ that may be used to describe
 * a functional relation \f$y_i=f(x_i)\f$. This class may be used to perform
 * a linear interpolation between the nodes to determine any value of
 * \f$y=f(x)\f$.
 * Nodes are allocated either from a double precision array, a GVector object
 * or a std::vector using the nodes() method. Alternatively, the node array
 * may be built on the fly using the append() method.
 * Interpolation can be either performed using the interpolate() method
 * or using the set_value(). In the latter case, the node indices and
 * weighting factors can be recovered using inx_left(), inx_right(),
 * wgt_left() and wgt_right().
 * If the nodes are equally spaced, interpolation is more rapid.
 ***************************************************************************/
class GNodeArray {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GNodeArray& array);

public:
    // Constructors and destructors
    GNodeArray(void);
    GNodeArray(const GNodeArray& array);
    ~GNodeArray(void);

    // Operators
    GNodeArray& operator= (const GNodeArray & array);

    // Methods
    void   nodes(const int& num, const double* array);
    void   nodes(const GVector& vector);
    void   nodes(const std::vector<double>& vector);
    void   append(const double& node);
    double interpolate(const double& value, const std::vector<double>& vector);
    void   set_value(const double& value);
    int    size(void) { return m_node.size(); }
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
    std::vector<double> m_node;          //!< Number of nodes
    std::vector<double> m_step;          //!< Distance to next node
    bool                m_is_linear;     //!< Nodes form a linear array
    double              m_linear_slope;  //!< Slope for linear array
    double              m_linear_offset; //!< Offset for linear array
    int                 m_inx_left;      //!< Index of left node for linear interpolation
    int                 m_inx_right;     //!< Index of right node for linear interpolation
    double              m_wgt_left;      //!< Weight for left node for linear interpolation
    double              m_wgt_right;     //!< Weight for right node for linear interpolation
};

#endif /* GNODEARRAY_HPP */
