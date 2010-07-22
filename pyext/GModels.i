/***************************************************************************
 *            GModels.i  -  Model container class SWIG interface           *
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
 * @file GModels.i
 * @brief GModels class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModels.hpp"
%}


/***********************************************************************//**
 * @class GModels
 *
 * @brief GModels class SWIG interface defintion.
 ***************************************************************************/
class GModels : public GOptimizerPars {
public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    ~GModels(void);
 
    // Methods
    void   append(const GModel& model);
    double eval(const GInstDir& obsDir, const GEnergy& obsEng,
                const GTime& obsTime, const GResponse& rsp,
                const GPointing& pnt);
    double eval_gradients(const GInstDir& obsDir, const GEnergy& obsEng,
                          const GTime& obsTime, const GResponse& rsp,
                          const GPointing& pnt);
    int    size(void) const { return m_elements; }
};


/***********************************************************************//**
 * @brief GModels class extension
 ***************************************************************************/
%extend GModels {
    char *__str__() {
        static char str_buffer[10001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 10001);
        str_buffer[10000] = '\0';
        return str_buffer;
    }
    GModel* __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return (*self)(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, 
                                       self->size());
    }
    GModels copy() {
        return (*self);
    }
};
