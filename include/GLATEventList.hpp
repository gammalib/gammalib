/***************************************************************************
 *                GLATEventList.hpp  -  LAT Event list class               *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATEventList.hpp
 * @brief GLATEventList class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATEVENTLIST_HPP
#define GLATEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventList.hpp"
#include "GLATEventAtom.hpp"
#include "GFits.hpp"


/***********************************************************************//**
 * @class GLATEventList
 *
 * @brief GLATEventList class interface defintion.
 ***************************************************************************/
class GLATEventList : public GEventList {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATEventList& list);

public:
    // Constructors and destructors
    GLATEventList();
    GLATEventList(const GLATEventList& list);
    virtual ~GLATEventList();

    // Operators
    GLATEventList& operator= (const GLATEventList& list);

    // Methods
	void           load(const std::string& filename);
    void           load(GFitsHDU* hdu);
    GLATEventAtom* pointer(int index) const;
    
protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GLATEventList& list);
    void           free_members(void);
    GLATEventList* clone(void) const;

    // Protected data area

private:
};

#endif /* GLATEVENTLIST_HPP */
