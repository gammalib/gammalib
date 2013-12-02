/***************************************************************************
 *                     GSkyPixel.hpp - Sky map pixel class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GSkyPixel.hpp
 * @brief Sky map pixel class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSKYPIXEL_HPP
#define GSKYPIXEL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GSkyPixel
 *
 * @brief Sky map pixel class
 *
 * This class implements a sky map pixel. A sky map pixel may be either
 * 1-dimensional (1D) or 2-dimensional (2D). 1D pixels are set and retrieved
 * using the index() method. The X and Y components of 2D pixels are set and
 * retrieved using the x() and y() methods. The xy() method allows setting
 * of both the X and Y components simultaneously.
 *
 * There are integer and double precision variants for 1D and 2D pixel
 * constructors:
 *
 *    GSkyPixel int_1D_pixel(1);          // Set index to 1
 *    GSkyPixel double_1D_pixel(1.0);     // Set index to 1.0
 *    GSkyPixel int_2D_pixel(1,2);        // Set (x,y)=(1,2)
 *    GSkyPixel double_2D_pixel(1.0,2.0); // Set (x,y)=(1.0,2.0)
 *
 * The 1D constructors allow for automatic type conversion:
 *
 *    GSkyPixel int_1D_pixel = 1;         // Set index to 1
 *    GSkyPixel double_1D_pixel = 1.0;    // Set index to 1.0
 *
 * Note that the class also provides conversion operators to int and double
 * allowing to write:
 *
 *    int    index  = GSkyPixel(1);       // Retrieve integer index
 *    double dindex = GSkyPixel(1.0);     // Retrieve double precision index
 *
 * The dimensionality of a pixel is checked using the is1D() and is2D()
 * methods. The size() method returns 1 for 1D pixels and 2 for 2D pixels.
 * If a pixel is not initialised, e.g.
 *
 *    GSkyPixel pixel;                    // Not initialised
 *
 * the is1D() and is2D() methods return false and the size() method returns
 * 0.
 ***************************************************************************/
class GSkyPixel : public GBase {

public:
    // Constructors and destructors
    GSkyPixel(void);
    GSkyPixel(const int& index);
    GSkyPixel(const double& index);
    GSkyPixel(const int& x, const int& y);
    GSkyPixel(const double& x, const double& y);
    GSkyPixel(const GSkyPixel& pixel);
    virtual ~GSkyPixel(void);

    // Operators
    GSkyPixel& operator=(const GSkyPixel& pixel);
    operator int() const;
    operator double() const;
    
    // Methods
    void          clear(void);
    GSkyPixel*    clone(void) const;
    int           size(void) const;
    bool          is1D(void) const;
    bool          is2D(void) const;
    void          index(const double& index);
    void          x(const double& x);
    void          y(const double& y);
    void          xy(const double& x, const double& y);
    const double& index(void) const;
    const double& x(void) const;
    const double& y(void) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GSkyPixel& pixel);
    void free_members(void);

    // Private data area
    int    m_size;       //!< Pixel dimension (0=undefined, 1=1D, 2=2D)
    double m_x;          //!< X index
    double m_y;          //!< Y index
};


/***********************************************************************//**
 * @brief Return pixel dimension
 *
 * @return Pixel dimension [0,1,2].
 *
 * Returns the pixel dimension. If the pixel is not defined the method
 * returns 0.
 ***************************************************************************/
inline
int GSkyPixel::size(void) const
{
    return m_size;
}


/***********************************************************************//**
 * @brief Check if pixel is 1D
 *
 * @return True if pixel is 1D.
 ***************************************************************************/
inline
bool GSkyPixel::is1D(void) const
{
    return (m_size == 1);
}


/***********************************************************************//**
 * @brief Check if pixel is 2D
 *
 * @return True if pixel is 2D.
 ***************************************************************************/
inline
bool GSkyPixel::is2D(void) const
{
    return (m_size == 2);
}


/***********************************************************************//**
 * @brief Set sky map pixel index
 *
 * @param[in] index Pixel index.
 *
 * Sets the 1D pixel index.
 ***************************************************************************/
inline
void GSkyPixel::index(const double& index)
{
    // Set pixel index
    m_size = 1;
    m_x    = index;
    m_y    = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set x value of sky map pixel
 *
 * @param[in] x X pixel index.
 *
 * Sets the X value of a 2D pixel index. The method does not alter the Y
 * pixel index.
 ***************************************************************************/
inline
void GSkyPixel::x(const double& x)
{
    m_size = 2;
    m_x    = x;
    return;
}


/***********************************************************************//**
 * @brief Set y value of sky map pixel
 *
 * @param[in] y Y pixel index.
 *
 * Sets the Y value of a 2D pixel index. The method does not alter the X
 * pixel index.
 ***************************************************************************/
inline
void GSkyPixel::y(const double& y)
{
    m_size = 2;
    m_y    = y;
    return;
}


/***********************************************************************//**
 * @brief Set x and y value of sky map pixel
 *
 * @param[in] x X pixel index.
 * @param[in] y Y pixel index.
 *
 * Sets the X and Y value of a 2D pixel index.
 ***************************************************************************/
inline
void GSkyPixel::xy(const double& x, const double& y)
{
    m_size = 2;
    m_x    = x;
    m_y    = y;
    return;
}


/***********************************************************************//**
 * @brief Return sky map pixel index
 *
 * @return Sky map pixel index
 *
 * Returns the 1D index of a sky map pixel. The method does not check whether
 * the pixel is indeed a 1D pixel.
 ***************************************************************************/
inline
const double& GSkyPixel::index(void) const
{
    return m_x;
}


/***********************************************************************//**
 * @brief Return x value of sky map pixel
 *
 * @return X value of sky map pixel
 *
 * Returns the X value of a 2D sky map pixel. The method does not check
 * whether the pixel is indeed a 2D pixel.
 ***************************************************************************/
inline
const double& GSkyPixel::x(void) const
{
    return m_x;
}


/***********************************************************************//**
 * @brief Return y value of sky pixel
 *
 * @return Y value of sky map pixel
 *
 * Returns the Y value of a 2D sky map pixel. The method does not check
 * whether the pixel is indeed a 2D pixel.
 ***************************************************************************/
inline
const double& GSkyPixel::y(void) const
{
    return m_y;
}

#endif /* GSKYPIXEL_HPP */
