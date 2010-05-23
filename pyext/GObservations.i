/***************************************************************************
 *     GObservations.i  -  Observations container class SWIG interface     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GObservations.i
 * @brief GObservations class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GObservations.hpp"
%}


/***********************************************************************//**
 * @class GObservations
 *
 * @brief Interface for the observations container class.
 ***************************************************************************/
class GObservations {
public:
    // Constructors and destructors
    GObservations(void);
    GObservations(const GObservations& obs);
    virtual ~GObservations(void);

    // Methods
    void     append(GObservation &obs);
    int      size(void) const { return m_num; }
    void     models(const GModels& models) { m_models=models; return; }
    GModels* models(void) { return &m_models; }
    void     optimize(GOptimizer& opt);

};


/***********************************************************************//**
 * @brief GObservations class extension
 ***************************************************************************/
%extend GObservations {
    char *__str__() {
        static char str_buffer[1001];
        std::ostringstream buffer;
        buffer << *self;
        std::string str = buffer.str();
        strncpy(str_buffer, (char*)str.c_str(), 1001);
        str_buffer[1000] = '\0';
        return str_buffer;
    }
    GObservation& __getitem__(int index) {
        if (index >= 0 && index < self->size())
            return (*self)(index);
        else
            throw GException::out_of_range("__getitem__(int)", index, 
                                           self->size());
    }
};
