/***************************************************************************
 * GModelTemporalConst.i  -  Temporal constant model class SWIG interface  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
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
 * @brief GModelTemporalConst class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelTemporalConst.hpp"
%}


/***********************************************************************//**
 * @class GModelTemporalConst
 *
 * @brief SWIG interface definition for the constant model class.
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
    GModelPar*           par(int index) const;
    double               eval(const GTime& srcTime);
    double               eval_gradients(const GTime& srcTime);
};


/***********************************************************************//**
 * @brief GModelTemporalConst class extension
 ***************************************************************************/
%extend GModelTemporalConst {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    GModelTemporalConst copy() {
        return (*self);
    }
};
