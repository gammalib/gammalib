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
#include "GCTAResponse.hpp"
#include "GCTAPointing.hpp"
#include "GTime.hpp"
#include "GModel.hpp"


/***********************************************************************//**
 * @class GCTAObservation
 *
 * @brief Interface class for CTA observations.
 ***************************************************************************/
class GCTAObservation : public GObservation {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAObservation& obs);

public:
    // Constructors and destructors
    GCTAObservation(void);
    GCTAObservation(const GCTAObservation& obs);
    virtual ~GCTAObservation(void);

    // Operators
    GCTAObservation& operator= (const GCTAObservation& obs);

    // Implement pure virtual methods
    GCTAObservation* clone(void) const;
    void             response(const std::string& irfname, std::string caldb = "");
    GResponse*       response(const GTime& time) const;
    GPointing*       pointing(const GTime& time) const;
    std::string      instrument(void) const;

    // Other methods
    void load_unbinned(const std::string& filename);
    void load_binned(const std::string& filename);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAObservation& obs);
    void free_members(void);

    // Npred integration methods
    double npred_temp(const GModel& model) const;
    double npred_grad_temp(const GModel& model, int ipar) const;

    // Protected data area
    GCTAResponse* m_response;     //!< Pointer to instrument response functions
    GCTAPointing* m_pointing;     //!< Pointer to pointing direction

};

#endif /* GCTAOBSERVATION_HPP */
