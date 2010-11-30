/***************************************************************************
 *          GEventCube.hpp  -  Abstract event cube container class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEventCube.hpp
 * @brief GEventCube container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTCUBE_HPP
#define GEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvents.hpp"
#include "GFits.hpp"


/***********************************************************************//**
 * @class GEventCube
 *
 * @brief GEventCube container class interface defintion.
 ***************************************************************************/
class GEventCube : public GEvents {

    // Friend classes
    friend class GData;
    friend class GObservation;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GEventCube& cube);

public:
    // Constructors and destructors
    GEventCube(void);
    GEventCube(const GEventCube& cube);
    virtual ~GEventCube(void);

    // Operators
    virtual GEventCube& operator= (const GEventCube& cube);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEventCube* clone(void) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual GEventBin*  pointer(int index) = 0;
    virtual int         number(void) const = 0;

    // Implemented pure virtual methods
    int  size(void) const { return m_elements; }
    int  dim(void) const { return m_dim; }
    int  naxis(int axis) const;
    bool islist(void) const { return false; }
    bool iscube(void) const { return true; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEventCube& cube);
    void free_members(void);

    // Protected data area
    int  m_elements;         //!< Number of cube elements
    int  m_dim;              //!< Cube dimension
    int* m_naxis;            //!< Number of bins in each axis

private:
};

#endif /* GEVENTCUBE_HPP */
