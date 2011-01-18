/***************************************************************************
 *   GModelTemporalConst.i  -  Temporal constant model class python I/F    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelTemporalConst.i
 * @brief GModelTemporalConst class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelTemporalConst.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelTemporalConst
 *
 * @brief Python interface definition for the constant model class.
 ***************************************************************************/
class GModelTemporalConst  : public GModelTemporal {
public:
    // Constructors and destructors
    GModelTemporalConst(void);
    GModelTemporalConst(const GModelTemporalConst& model);
    virtual ~GModelTemporalConst(void);

    // Implemented virtual methods
    void                 clear(void);
    GModelTemporalConst* clone(void) const;
    int                  size(void) const { return m_npars; }
    std::string          type(void) const { return "Constant"; }
    double               eval(const GTime& srcTime);
    double               eval_gradients(const GTime& srcTime);
    GTimes               mc(const double& rate,
                            const GTime& tmin, const GTime& tmax,
                            GRan& ran);
    void                 read(const GXmlElement& xml);
    void                 write(GXmlElement& xml) const;

    // Other methods
    double norm(void) const { return m_norm.real_value(); }
};


/***********************************************************************//**
 * @brief GModelTemporalConst class extension
 ***************************************************************************/
%extend GModelTemporalConst {
    char *__str__() {
        return tochar(self->print());
    }
    GModelTemporalConst copy() {
        return (*self);
    }
};
