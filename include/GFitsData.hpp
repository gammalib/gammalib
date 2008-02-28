/***************************************************************************
 *         GFitsData.hpp  - FITS data handling abstract base class         *
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

#ifndef GFITSDATA_HPP
#define GFITSDATA_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                        GFitsData class definition                       *
 ***************************************************************************/
class GFitsData {

public:
    // Constructors and destructors
    GFitsData();
    GFitsData(const GFitsData& data);
    virtual ~GFitsData();

    // Operators
    virtual GFitsData& operator= (const GFitsData& data);

    // Methods
    virtual void       open(__fitsfile*  fptr) = 0;
    virtual void       close(void) = 0;
    virtual GFitsData* clone(void) const = 0;
    
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GFitsData& data);
    void free_members(void);

    // Protected data area
};

#endif /* GFITSDATA_HPP */
