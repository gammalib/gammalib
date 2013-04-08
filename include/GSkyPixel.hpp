/***************************************************************************
 *        GSkyPixel.cpp - Class that implements a 2D sky pixel index       *
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
 * @brief Sky pixel class definition
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
 * @brief GSkyPixel class interface defintion
 ***************************************************************************/
class GSkyPixel : public GBase {

public:
    // Constructors and destructors
    GSkyPixel(void);
    explicit GSkyPixel(double x, double y);
    GSkyPixel(const GSkyPixel& pixel);
    virtual ~GSkyPixel(void);

    // Operators
    GSkyPixel& operator= (const GSkyPixel& pixel);
    
    // Methods
    void        clear(void);
    GSkyPixel*  clone(void) const;
    void        x(const double& x);
    void        y(const double& y);
    double      x(void) const;
    double      y(void) const;
    std::string print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GSkyPixel& pixel);
    void free_members(void);

    // Private data area
    double m_x;          //!< WCS x index
    double m_y;          //!< WCS y index
};

#endif /* GSKYPIXEL_HPP */
