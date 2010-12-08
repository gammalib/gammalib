/***************************************************************************
 *                GLATEventList.hpp  -  LAT Event list class               *
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
 * @file GLATEventList.hpp
 * @brief GLATEventList class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATEVENTLIST_HPP
#define GLATEVENTLIST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventList.hpp"
#include "GLATEventAtom.hpp"
#include "GLATObservation.hpp"
#include "GLATResponse.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GLATEventList
 *
 * @brief GLATEventList class interface defintion.
 ***************************************************************************/
class GLATEventList : public GEventList {

public:
    // Constructors and destructors
    GLATEventList(void);
    GLATEventList(const GLATEventList& list);
    virtual ~GLATEventList(void);

    // Operators
    GLATEventList& operator= (const GLATEventList& list);

    // Implemented pure virtual base class methods
    void           clear(void);
    GLATEventList* clone(void) const;
    int            size(void) const { return m_num; }
    void           load(const std::string& filename);
    GLATEventAtom* pointer(int index);
    int            number(void) const { return m_num; }
    std::string    print(void) const;

    // Other methods

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATEventList& list);
    void free_members(void);
    void load_ft1(GFitsTable* hdu);
    void load_events(GFitsTable* hdu);
    void load_ds_keys(GFitsTable* hdu);

    // Protected data area
    int            m_num;            //!< Number of events
    GLATEventAtom* m_events;         //!< Pointer to events

    // Diffuse response information
    int            m_num_difrsp;     //!< Number of diffuse response models
    std::string*   m_difrsp_label;   //!< Diffuse response model labels

    // Data selection information
    int            m_ds_num;         //!< Number of data selection keys
    std::string*   m_ds_type;        //!< Data selection types
    std::string*   m_ds_unit;        //!< Data selection units
    std::string*   m_ds_value;       //!< Data selection values
    std::string*   m_ds_reference;   //!< Data selection references

private:
};

#endif /* GLATEVENTLIST_HPP */
