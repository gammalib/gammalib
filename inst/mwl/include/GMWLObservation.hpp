/***************************************************************************
 *       GMWLObservation.hpp  -  Multi-wavelength observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GMWLObservation.hpp
 * @brief GMWLObservation class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMWLOBSERVATION_HPP
#define GMWLOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"
#include "GObservation.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GMWLResponse.hpp"
#include "GMWLPointing.hpp"


/***********************************************************************//**
 * @class GMWLObservation
 *
 * @brief Interface class for multi-wavelength observations
 *
 * This class implements a multi-wavelength observation. A multi-wavelength
 * observation contains spectral points obtained with an unspecified
 * instrument. The spectral points are given in physical units.
 ***************************************************************************/
class GMWLObservation : public GObservation {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GMWLObservation& obs);
    friend GLog&         operator<< (GLog& log, const GMWLObservation& obs);

public:
    // Constructors and destructors
    GMWLObservation(void);
    GMWLObservation(const std::string& filename);
    GMWLObservation(const GMWLObservation& obs);
    virtual ~GMWLObservation(void);

    // Operators
    GMWLObservation& operator= (const GMWLObservation& obs);

    // Implement pure virtual methods
    void             clear(void);
    GMWLObservation* clone(void) const;
    void             response(const std::string& rspname, std::string caldb = "");
    GMWLResponse*    response(const GTime& time) const;
    GMWLPointing*    pointing(const GTime& time) const;
    std::string      instrument(void) const { return m_instrument; }

    // Other methods
    void        load(const std::string& filename);
    void        instrument(const std::string& instrument) { m_instrument=instrument; }
    std::string print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLObservation& obs);
    void free_members(void);

    // Protected members
    std::string   m_instrument;   //!< Instrument name
    GMWLResponse* m_response;     //!< Pointer to response functions
    GMWLPointing* m_pointing;     //!< Pointer to pointing direction
};

#endif /* GMWLOBSERVATION_HPP */
