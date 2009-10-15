/***************************************************************************
 *                   GData.hpp  -  Data container class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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
 * @brief GData container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GDATA_HPP
#define GDATA_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GEvent.hpp"


/***********************************************************************//**
 * @class GData
 *
 * @brief GData container class interface defintion.
 ***************************************************************************/
class GData {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GData& data);

public:
    // Constructors and destructors
    GData();
    GData(const GData& data);
    virtual ~GData();

    // Operators
    GData& operator= (const GData& data);

    // Methods
	void          add(GObservation &obs);
    void          models(const GModels& models) { m_models=models; return; }
    int           elements(void) const { return m_num; }
    GObservation* observation(int index) const;
    GModels*      models(void) { return &m_models; }

    // Event iterator
    class iterator {
    friend class GData;
    public:
        iterator();
        iterator(GData *data);
        ~iterator() { return; }
        iterator& operator++(void);                // Prefix
        iterator  operator++(int junk);            // Postfix
        bool      operator==(const iterator& it) const 
                  { return ((m_index == it.m_index) && (m_event == it.m_event)); }
        bool      operator!=(const iterator& it) const
                  { return ((m_index != it.m_index) || (m_event != it.m_event)); }
        GEvent&   operator*(void) { return *m_event; }
        GEvent*   operator->(void) { return &(*m_event); }
    protected:
        int               m_index;   //!< Actual observation index [0,m_num-1]
        GEvents::iterator m_event;   //!< Iterator on actual event
        GEvents::iterator m_end;     //!< Iterator on observation end
        GObservation*     m_obs;     //!< Pointer to actual observation
        GData*            m_data;    //!< Pointer to GData object
    };
    iterator begin(void);
    iterator end(void);

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GData& data);
    void           free_members(void);

    // Protected data area
	int            m_num;            //!< Number of observations
	GObservation** m_obs;            //!< Pointers to observations
    GModels        m_models;         //!< Models
	
};

#endif /* GDATA_HPP */
