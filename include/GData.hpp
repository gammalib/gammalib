/***************************************************************************
 *                   GData.hpp  -  Data container class                    *
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
	void add(GObservation &obs);

    // Event iterator
    class iterator {
    friend class GData;
    public:
        iterator(GData *data);
        ~iterator();
        iterator& operator++(void);                // Prefix
        iterator  operator++(int junk);            // Postfix
        bool      operator==(const iterator& it) const;
        bool      operator!=(const iterator& it) const;
        GEvent&   operator*(void);
        GEvent*   operator->(void);
    protected:
        int           m_obs_index;    //!< Actual observation index [0,m_num-1]
        int           m_event_index;  //!< Actual event index
        int           m_num_events;   //!< Total number of events in actual observation
        GObservation *m_obs;          //!< Pointer to actual observation
        GData        *m_data;         //!< Pointer to GData object
    };
    iterator begin(void);
    iterator end(void);

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GData& data);
    void           free_members(void);

    // Protected data area
	int            m_num;       //!< Number of observations
	GObservation **m_obs;       //!< Pointers to observations
	
private:
};

#endif /* GDATA_HPP */
