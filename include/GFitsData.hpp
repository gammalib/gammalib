/***************************************************************************
 *         GFitsData.hpp  - FITS data handling abstract base class         *
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
 * @file GFitsData.hpp
 * @brief GFitsData class definition.
 * @author J. Knodlseder
 */

#ifndef GFITSDATA_HPP
#define GFITSDATA_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GFitsData
 *
 * @brief Interface for abstract FITS data class
 ***************************************************************************/
class GFitsData {

    // Friend classes
    friend class GFitsHDU;

public:
    // Constructors and destructors
    GFitsData(void);
    GFitsData(const GFitsData& data);
    virtual ~GFitsData(void);

    // Operators
    virtual GFitsData& operator= (const GFitsData& data);

protected:
    // Protected methods
    virtual void       open(void* vptr) = 0;
    virtual void       save(void) = 0;
    virtual void       close(void) = 0;
    virtual GFitsData* clone(void) const = 0;
    virtual void       connect(void* vptr);

    // Private methods
    void init_members(void);
    void copy_members(const GFitsData& data);
    void free_members(void);

    // Protected data area
    void* m_fitsfile;         //!< Pointer on FITS file data
};

#endif /* GFITSDATA_HPP */
