/***************************************************************************
 *     GCTAModelRadialAcceptance.hpp  -  Radial acceptance model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAModelRadialAcceptance.hpp
 * @brief GCTAModelRadialAcceptance class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAMODELRADIALACCEPTANCE_HPP
#define GCTAMODELRADIALACCEPTANCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelData.hpp"
#include "GModelPar.hpp"
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GXmlElement.hpp"
#include "GCTAEventList.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialAcceptance
 *
 * @brief Radial acceptance model interface definition.
 *
 * This class implements a radial acceptance model for CTA.
 ***************************************************************************/
class GCTAModelRadialAcceptance : public GModelData {

public:
    // Constructors and destructors
    GCTAModelRadialAcceptance(void);
    explicit GCTAModelRadialAcceptance(const GXmlElement& xml);
    GCTAModelRadialAcceptance(const GCTAModelRadialAcceptance& model);
    virtual ~GCTAModelRadialAcceptance(void);

    // Operators
    GCTAModelRadialAcceptance& operator= (const GCTAModelRadialAcceptance& model);

    // Implemented pure virtual methods
    void                       clear(void);
    GCTAModelRadialAcceptance* clone(void) const;
    std::string                type(void) const { return "RadialAcceptance"; }
    double                     eval(const GEvent& event, const GObservation& obs);
    double                     eval_gradients(const GEvent& event, const GObservation& obs);
    GCTAEventList*             mc(const GObservation& obs, GRan& ran);
    void                       read(const GXmlElement& xml);
    void                       write(GXmlElement& xml) const;
    std::string                print(void) const;

    // Other methods
    GCTAModelRadial* radial(void)   const { return m_radial; }
    GModelSpectral*  spectral(void) const { return m_spectral; }
    GModelTemporal*  temporal(void) const { return m_temporal; }

protected:
    // Protected methods
    void             init_members(void);
    void             copy_members(const GCTAModelRadialAcceptance& model);
    void             free_members(void);
    void             set_pointers(void);
    bool             valid_model(void) const;
    GCTAModelRadial* xml_radial(const GXmlElement& radial) const;
    GModelSpectral*  xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal*  xml_temporal(const GXmlElement& temporal) const;

    // Proteced data members
    GCTAModelRadial* m_radial;       //!< Radial model
    GModelSpectral*  m_spectral;     //!< Spectral model
    GModelTemporal*  m_temporal;     //!< Temporal model
};

#endif /* GCTAMODELRADIALACCEPTANCE_HPP */
