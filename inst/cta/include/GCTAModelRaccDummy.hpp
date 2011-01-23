/***************************************************************************
 *     GCTAModelRaccDummy.hpp  -  Radial acceptance dummy model class      *
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
 * @file GCTAModelRaccDummy.hpp
 * @brief GCTAModelRaccDummy class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAMODELRACCDUMMY_HPP
#define GCTAMODELRACCDUMMY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelData.hpp"
#include "GModelPar.hpp"
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GCTAModelRaccDummy
 *
 * @brief Radial acceptance dummy model interface definition.
 *
 * This class implements a radial acceptance dummy model for CTA.
 ***************************************************************************/
class GCTAModelRaccDummy : public GModelData {

public:
    // Constructors and destructors
    GCTAModelRaccDummy(void);
    explicit GCTAModelRaccDummy(const GXmlElement& xml);
    GCTAModelRaccDummy(const GCTAModelRaccDummy& model);
    virtual ~GCTAModelRaccDummy(void);

    // Operators
    GCTAModelRaccDummy& operator= (const GCTAModelRaccDummy& model);

    // Implemented pure virtual methods
    void                clear(void);
    GCTAModelRaccDummy* clone(void) const;
    std::string         type(void) const { return "RadialAcceptanceDummy"; }
    double              eval(const GEvent& event, const GObservation& obs);
    double              eval_gradients(const GEvent& event, const GObservation& obs);
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
    std::string         print(void) const;

    // Other methods

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRaccDummy& model);
    void free_members(void);
    void set_pointers(void);

    // Protected members
};

#endif /* GCTAMODELRACCDUMMY_HPP */
