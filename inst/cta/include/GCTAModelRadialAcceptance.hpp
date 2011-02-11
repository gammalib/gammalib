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
#include <cmath>
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
    double                     eval(const GEvent& event, const GObservation& obs) const;
    double                     eval_gradients(const GEvent& event, const GObservation& obs) const;
    double                     npred(const GEnergy& obsEng, const GTime& obsTime,
                                     const GObservation& obs) const;
    GCTAEventList*             mc(const GObservation& obs, GRan& ran) const;
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

    // ROI integration kernel
    class roi_kern : public GIntegrand {
    public:
        roi_kern(const GCTAModelRadial* parent, const double& roi, const double& dist) :
                 m_parent(parent),
                 m_roi(roi),   m_cosroi(cos(roi)),
                 m_dist(dist), m_cosdist(cos(dist)), m_sindist(sin(dist)) { return; }
        double eval(double r);
    protected:
        const GCTAModelRadial* m_parent;   //!< Pointer to radial model
        double                 m_roi;      //!< ROI radius in radians
        double                 m_cosroi;   //!< Cosine of ROI radius
        double                 m_dist;     //!< Distance between pointing and ROI centre in radians
        double                 m_cosdist;  //!< Cosine of distance
        double                 m_sindist;  //!< Sinus of distance
    };

    // Proteced data members
    GCTAModelRadial* m_radial;       //!< Radial model
    GModelSpectral*  m_spectral;     //!< Spectral model
    GModelTemporal*  m_temporal;     //!< Temporal model
};

#endif /* GCTAMODELRADIALACCEPTANCE_HPP */
