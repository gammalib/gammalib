/***************************************************************************
 *               GCTAObservation.hpp  -  CTA Observation class             *
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
 * @file GCTAObservation.hpp
 * @brief GCTAObservation class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAOBSERVATION_HPP
#define GCTAOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GCTAPointing.hpp"
#include "GCTAResponse.hpp"
#include "GTime.hpp"
#include "GModel.hpp"


/***********************************************************************//**
 * @class GCTAObservation
 *
 * @brief Interface class for CTA observations.
 ***************************************************************************/
class GCTAObservation : public GObservation {

public:
    // Constructors and destructors
    GCTAObservation(void);
    GCTAObservation(const GCTAObservation& obs);
    virtual ~GCTAObservation(void);

    // Operators
    GCTAObservation& operator= (const GCTAObservation& obs);

    // Implemented pure virtual base class methods
    void             clear(void);
    GCTAObservation* clone(void) const;
    GCTAResponse*    response(void) const;
    GCTAPointing*    pointing(const GTime& time) const;
    std::string      instrument(void) const;
    std::string      print(void) const;

    // Other methods
    void load_unbinned(const std::string& filename);
    void load_binned(const std::string& filename);
    void response(const std::string& irfname, std::string caldb = "");
    void pointing(const GCTAPointing& pointing);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAObservation& obs);
    void free_members(void);

    // Npred integration methods
    double npred_temp(const GModel& model) const;
    double npred_grad_temp(const GModel& model, int ipar) const;

    // Protected members
    GCTAResponse* m_response;   //!< Pointer to instrument response functions
    GCTAPointing* m_pointing;   //!< Pointer to pointing direction
};

#endif /* GCTAOBSERVATION_HPP */
