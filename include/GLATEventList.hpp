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
    void           load_events(GFitsHDU* hdu);
    void           load_ds_keys(GFitsHDU* hdu);

    // Diffuse response information
    int          m_num_difrsp;       //!< Number of diffuse response models
    std::string* m_difrsp_label;     //!< Diffuse response model labels
    
    // Data selection information
    int          m_ds_num;           //!< Number of data selection keys
    std::string* m_ds_type;          //!< Data selection types
    std::string* m_ds_unit;          //!< Data selection units
    std::string* m_ds_value;         //!< Data selection values
    std::string* m_ds_reference;     //!< Data selection references

private:
};

#endif /* GLATEVENTLIST_HPP */
