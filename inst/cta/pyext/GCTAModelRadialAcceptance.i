/***************************************************************************
 * GCTAModelRadialAcceptance.i  -  Radial acceptance model class python I/F*
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
 * @file GCTAModelRadialAcceptance.i
 * @brief GCTAModelRadialAcceptance class python interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialAcceptance.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialAcceptance
 *
 * @brief Radial acceptance model interface definition.
 ***************************************************************************/
class GCTAModelRadialAcceptance : public GModelData {
public:
    // Constructors and destructors
    GCTAModelRadialAcceptance(void);
    explicit GCTAModelRadialAcceptance(const GXmlElement& xml);
    GCTAModelRadialAcceptance(const GCTAModelRadialAcceptance& model);
    virtual ~GCTAModelRadialAcceptance(void);

    // Implemented pure virtual methods
    void                       clear(void);
    GCTAModelRadialAcceptance* clone(void) const;
    std::string                type(void) const;
    double                     eval(const GEvent& event, const GObservation& obs);
    double                     eval_gradients(const GEvent& event, const GObservation& obs);
    double                     npred(const GEnergy& obsEng, const GTime& obsTime,
                                     const GObservation& obs) const;
    GCTAEventList*             mc(const GObservation& obs, GRan& ran) const;
    void                       read(const GXmlElement& xml);
    void                       write(GXmlElement& xml) const;

    // Other methods
    GCTAModelRadial* radial(void)   const;
    GModelSpectral*  spectral(void) const;
    GModelTemporal*  temporal(void) const;
};


/***********************************************************************//**
 * @brief GCTAModelRadial class extension
 ***************************************************************************/
%extend GCTAModelRadial {
};
