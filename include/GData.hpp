/***************************************************************************
 *                  GData.hpp  -  Data abstract base class                 *
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
 * @file GData.hpp
 * @brief GData abstract base class definition.
 * @author J. Knodlseder
 */

#ifndef GDATA_HPP
#define GDATA_HPP

/* __ Includes ___________________________________________________________ */

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GData
 *
 * @brief Abstract GData class interface defintion
 ***************************************************************************/
class GData {

public:
    // Constructors and destructors
    GData();
    GData(const GData& d);
    virtual ~GData();

    // Operators
    GData& operator= (const GData& data);

    // Methods
    
protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GData& data);
    void           free_members(void);
    virtual GData* clone(void) const = 0;

    // Protected data area
private:
};

#endif /* GDATA_HPP */
