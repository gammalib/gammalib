/***************************************************************************
 *       GSkyPixel.cpp  -  Class that implements a 2D sky pixel index      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GSkyPixel.hpp
 * @brief GSkyPixel class definition.
 * @author J. Knodlseder
 */

#ifndef GSKYPIXEL_HPP
#define GSKYPIXEL_HPP

/* __ Includes ___________________________________________________________ */
#include <sstream>


/***********************************************************************//**
 * @class GSkyPixel
 *
 * @brief GSkyPixel class interface defintion
 ***************************************************************************/
class GSkyPixel {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GSkyPixel& pixel);

public:
    // Constructors and destructors
    GSkyPixel(void);
    GSkyPixel(const GSkyPixel& pixel);
    virtual ~GSkyPixel(void);

    // Operators
    GSkyPixel& operator= (const GSkyPixel& pixel);
    
    // Methods
    void   x(const double& x);
    void   y(const double& y);
    double x(void) const;
    double y(void) const;

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
