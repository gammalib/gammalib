/***************************************************************************
 *                        GData.hpp  -  Data class                         *
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
 * @brief GData class definition.
 * @author J. Knodlseder
 */

#ifndef GDATA_HPP
#define GDATA_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"


/***********************************************************************//**
 * @class GData
 *
 * @brief GData class interface defintion
 ***************************************************************************/
class GData {

public:
    // Constructors and destructors
    GData();
    GData(const GData& data);
    virtual ~GData();

    // Operators
    GData& operator= (const GData& data);

    // Methods
    
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
