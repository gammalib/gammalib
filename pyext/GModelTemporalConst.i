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

    // Methods
    int        npars(void) const { return m_npars; }
    GModelPar* par(int index) const;
    double     eval(const GTime& srcTime);
    double     eval_gradients(const GTime& srcTime);
};


/***********************************************************************//**
 * @brief GModelTemporalConst class extension
 ***************************************************************************/
%extend GModelTemporalConst {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
        str_buffer[1000] = '\0';
        return str_buffer;
    }
    GModelTemporalConst copy() {
        return (*self);
    }
};
