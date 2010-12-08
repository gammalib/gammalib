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
    void     clear(void);
    int      size(void) const { return m_num; }
    void     append(GObservation &obs);
    void     models(const GModels& models) { m_models=models; return; }
    void     models(const std::string& filename);
    GModels* models(void) { return &m_models; }
    void     optimize(GOptimizer& opt);
    double   npred(void) const { return m_npred; }
};


/***********************************************************************//**
 * @brief GObservations class extension
 ***************************************************************************/
%extend GObservations {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    GObservation& __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return (*self)(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GObservation& val) {
        if (index>=0 && index < self->size())
            (*self)(index) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    GObservations copy() {
        return (*self);
    }
};
