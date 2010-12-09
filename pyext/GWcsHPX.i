/***************************************************************************
 *              GWcsHPX.i  -  Healpix projection class SWIG file           *
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
 * @file GWcsHPX.i
 * @brief GWcsHPX class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GWcsHPX.hpp"
%}


/***********************************************************************//**
 * @class GWcsHPX
 *
 * @brief GWcsHPX class interface defintion
 ***************************************************************************/
class GWcsHPX : public GWcs {
public:
    // Constructors and destructors
    GWcsHPX(void);
    explicit GWcsHPX(const int& nside, const std::string& ordering = "NESTED",
                     const std::string& coordsys = "GAL");
    explicit GWcsHPX(const GFitsHDU* hdu);
    GWcsHPX(const GWcsHPX& wcs);
    virtual ~GWcsHPX(void);

    // Implemented pure virtual methods
    void     clear(void);
    GWcsHPX* clone(void) const;
    void     read(const GFitsHDU* hdu);
    void     write(GFitsHDU* hdu) const;

    // Overloaded base class methods
    double      omega(const int& pix) const;
    GSkyDir     pix2dir(const int& pix);
    int         dir2pix(GSkyDir dir) const;

    // Class specific methods
    int         npix(void) const;
    int         nside(void) const;
    std::string ordering(void) const;
    void        ordering(const std::string& ordering);
};


/***********************************************************************//**
 * @brief GWcsHPX class extension
 ***************************************************************************/
%extend GWcsHPX {
    GWcsHPX copy() {
        return (*self);
    }
};
