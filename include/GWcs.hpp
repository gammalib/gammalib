/***************************************************************************
 *          GWcs.hpp  -  World Coordinate System virtual base class        *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2010 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GWcs.hpp
 * @brief GWcs virtual base class definition.
 * @author J. Knodlseder
 */

#ifndef GWCS_HPP
#define GWCS_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsHDU.hpp"
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GWcs
 *
 * @brief GWcs virtual base class interface defintion
 ***************************************************************************/
class GWcs {

  // Friend classes
  friend class GSkymap;

public:
    // Constructors and destructors
    GWcs();
    GWcs(const GWcs& wcs);
    virtual ~GWcs();

    // Operators
    virtual GWcs& operator= (const GWcs& wcs);

    // Virtual Methods
    virtual std::string type(void) const = 0;
    virtual void        read(const GFitsHDU* hdu) = 0;
    virtual void        write(GFitsHDU* hdu) = 0;
    virtual GSkyDir     pix2dir(const int& pix) = 0;
    virtual int         dir2pix(GSkyDir dir) const = 0;
    virtual double      omega(const int& pix) const = 0;
    virtual int         npix(void) const = 0;
    virtual int         naxes(void) const = 0;
    virtual int         naxis(const int& axis) const = 0;

    // Implemented methods
    std::string coordsys(void) const;
    void        coordsys(const std::string& coordsys);

private:
    // Private methods
    void          init_members(void);
    void          copy_members(const GWcs& wcs);
    void          free_members(void);
    virtual GWcs* clone(void) const = 0;

protected:
    // Protected data area
    int      m_coordsys;     //!< 0=celestial, 1=galactic
};

#endif /* GWCS_HPP */
