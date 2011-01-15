/***************************************************************************
 *          GWcs.i  -  World Coordinate System base class SWIG file        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GWcs.i
 * @brief GWcs class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcs.hpp"
#include "GTools.hpp"
%}


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
    GWcs(void);
    explicit GWcs(const std::string& coords,
                  const double& crval1, const double& crval2,
                  const double& crpix1, const double& crpix2,
                  const double& cdelt1, const double& cdelt2);
    GWcs(const GWcs& wcs);
    virtual ~GWcs(void);

    // Pure virtual methods (not implemented)
    virtual void  clear(void) = 0;
    virtual GWcs* clone(void) const = 0;
    virtual void  read(const GFitsHDU* hdu) = 0;
    virtual void  write(GFitsHDU* hdu) const = 0;

    // Virtual methods
    virtual std::string type(void) const;
    virtual std::string coordsys(void) const;
    virtual void        coordsys(const std::string& coordsys);
    virtual double      omega(const int& pix) const;
    virtual double      omega(const GSkyPixel& pix) const;
    virtual GSkyDir     pix2dir(const int& pix) const;
    virtual int         dir2pix(GSkyDir dir) const;
    virtual GSkyDir     xy2dir(const GSkyPixel& pix) const;
    virtual GSkyPixel   dir2xy(GSkyDir dir) const;
};




/***********************************************************************//**
 * @brief GWcs class extension
 ***************************************************************************/
%extend GWcs {
    char *__str__() {
        return tochar(self->print());
    }
    bool __is__(const GWcs &a) {
            return (*self) == a;
    }
};
