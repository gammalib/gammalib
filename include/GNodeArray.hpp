/***************************************************************************
 *                  GNodeArray.hpp - Array of nodes class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2014 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GNodeArray.hpp
 * @brief Node array class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GNODEARRAY_HPP
#define GNODEARRAY_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GContainer.hpp"
#include "GVector.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GNodeArray
 *
 * @brief Node array class
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
class GNodeArray : public GContainer {

public:
    // Constructors and destructors
    GNodeArray(void);
    explicit GNodeArray(const int& num, const double* array);
    explicit GNodeArray(const GVector& vector);
    explicit GNodeArray(const std::vector<double>& vector);
    GNodeArray(const GNodeArray& array);
    virtual ~GNodeArray(void);

    // Operators
    GNodeArray&   operator= (const GNodeArray & array);
    double&       operator[](const int& index);
    const double& operator[](const int& index) const;

    // Methods
    void          clear(void);
    GNodeArray*   clone(void) const;
    std::string   classname(void) const;
    double&       at(const int& index);
    const double& at(const int& index) const;
    int           size(void) const;
    bool          is_empty(void) const;
    void          append(const double& node);
    void          insert(const int& index, const double& node);
    void          remove(const int& index);
    void          reserve(const int& num);
    void          extend(const GNodeArray& nodes);
    void          nodes(const int& num, const double* array);
    void          nodes(const GVector& vector);
    void          nodes(const std::vector<double>& vector);
    double        interpolate(const double& value,
                              const std::vector<double>& vector) const;
    void          set_value(const double& value) const;
    const int&    inx_left(void) const;
    const int&    inx_right(void) const;
    const double& wgt_left(void) const;
    const double& wgt_right(void) const;
    void          load(const std::string& filename,
                       const std::string& extname = "NODES");
    void          save(const std::string& filename, const bool& clobber = false,
                       const std::string& extname = "NODES") const;
    void          read(const GFitsTable& table);
    void          write(GFits& file,
                        const std::string& extname = "NODES") const;
    std::string   print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GNodeArray& array);
    void free_members(void);
    void setup(void) const;
    
    // Node values
    std::vector<double> m_node;                   //!< Array of nodes

    // Evaluation cache
    mutable bool                m_need_setup;     //!< Call of setup is required
    mutable bool                m_is_linear;      //!< Nodes form a linear array
    mutable bool                m_has_last_value; //!< Last value is valid
    mutable std::vector<double> m_step;           //!< Distance to next node
    mutable double              m_last_value;     //!< Last requested value
    mutable double              m_linear_slope;   //!< Slope for linear array
    mutable double              m_linear_offset;  //!< Offset for linear array
    mutable int                 m_inx_left;       //!< Index of left node for linear interpolation
    mutable int                 m_inx_right;      //!< Index of right node for linear interpolation
    mutable double              m_wgt_left;       //!< Weight for left node for linear interpolation
    mutable double              m_wgt_right;      //!< Weight for right node for linear interpolation
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GNodeArray").
 ***************************************************************************/
inline
std::string GNodeArray::classname(void) const
{
    return ("GNodeArray");
}


/***********************************************************************//**
 * @brief Node access operator
 *
 * @param[in] index Node index [0,...,size()-1].
 * @return Node value.
 *
 * Returns a reference to the node with the specified @p index. No range
 * checking is performed on @p index. As this operator may change the
 * values of the node array, the setup method needs to be called before
 * doing the interpolation.
 ***************************************************************************/
inline
double& GNodeArray::operator[](const int& index)
{
    m_need_setup = true;
    return (m_node[index]);
}


/***********************************************************************//**
 * @brief Node access operator (const version)
 *
 * @param[in] index Node index [0,...,size()-1].
 * @return Node value.
 *
 * Returns a reference to the node with the specified @p index. No range
 * checking is performed on @p index.
 ***************************************************************************/
inline
const double& GNodeArray::operator[](const int& index) const
{
    return (m_node[index]);
}


/***********************************************************************//**
 * @brief Return number of nodes in node array
 *
 * @return Number of nodes in node array.
 *
 * Returns the number of nodes in the node array.
 ***************************************************************************/
inline
int GNodeArray::size(void) const
{
    return (m_node.size());
}


/***********************************************************************//**
 * @brief Signals if there are no nodes in node array
 *
 * @return True if node array is empty, false otherwise.
 *
 * Signals if the node array does not contain any node.
 ***************************************************************************/
inline
bool GNodeArray::is_empty(void) const
{
    return (m_node.empty());
}


/***********************************************************************//**
 * @brief Reserves space for nodes in node array
 *
 * @param[in] num Number of nodes.
 *
 * Reserves space for @p num nodes in the node array.
 ***************************************************************************/
inline
void GNodeArray::reserve(const int& num)
{
    m_node.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Returns left node index
 *
 * @return Left node index.
 *
 * Returns the left node index to be used for interpolation.
 ***************************************************************************/
inline
const int& GNodeArray::inx_left(void) const
{
    return m_inx_left;
}


/***********************************************************************//**
 * @brief Returns right node index
 *
 * @return Right node index.
 *
 * Returns the right node index to be used for interpolation.
 ***************************************************************************/
inline
const int& GNodeArray::inx_right(void) const
{
    return m_inx_right;
}


/***********************************************************************//**
 * @brief Returns left node weight
 *
 * @return Left node weight.
 *
 * Returns the weighting factor for the left node to be used for
 * interpolation.
 ***************************************************************************/
inline
const double& GNodeArray::wgt_left(void) const
{
    return m_wgt_left;
}


/***********************************************************************//**
 * @brief Returns right node weight
 *
 * @return Right node weight.
 *
 * Returns the weighting factor for the right node to be used for
 * interpolation.
 ***************************************************************************/
inline
const double& GNodeArray::wgt_right(void) const
{
    return m_wgt_right;
}

#endif /* GNODEARRAY_HPP */
