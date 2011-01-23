/***************************************************************************
 * GCTAModelRaccDummy.i  -  Radial acceptance dummy model class python I/F *
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
 * @file GCTAModelRaccDummy.i
 * @brief GCTAModelRaccDummy class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRaccDummy.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRaccDummy
 *
 * @brief Radial acceptance dummy model interface definition
 ***************************************************************************/
class GCTAModelRaccDummy : public GModelData, public GModelFactorized {

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
    std::string         type(void) const;
    double              eval(const GEvent& event, const GObservation& obs);
    double              eval_gradients(const GEvent& event, const GObservation& obs);
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GCTAModelRaccDummy class extension
 ***************************************************************************/
%extend GCTAModelRaccDummy {
    GCTAModelRaccDummy copy() {
        return (*self);
    }
};
