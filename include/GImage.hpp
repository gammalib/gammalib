/***************************************************************************
 *                 GImage.hpp  -  Image abstract base class                *
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
 * @file GImage.hpp
 * @brief GImage abstract base class definition.
 * @author J. Knodlseder
 */

#ifndef GIMAGE_HPP
#define GIMAGE_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GImage
 *
 * @brief Abstract GImage class interface defintion
 ***************************************************************************/
class GImage {

public:
    // Constructors and destructors
    GImage();
    GImage(const GImage& image);
    virtual ~GImage();

    // Operators
    GImage& operator= (const GImage& image);

    // Methods
    
protected:
    // Protected methods
    void            init_members(void);
    void            copy_members(const GImage& image);
    void            free_members(void);
    virtual GImage* clone(void) const = 0;

    // Protected data area
private:
};

#endif /* GIMAGE_HPP */
