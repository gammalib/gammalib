/***************************************************************************
 *            GModelPar.i  -  Model parameter class SWIG interface         *
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
 * @file GModelPar.i
 * @brief GModelPar class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelPar.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelPar
 *
 * @brief GModelPar class SWIG interface defintion.
 ***************************************************************************/
class GModelPar {
public:
    // Constructors and destructors
    GModelPar(void);
    GModelPar(const GModelPar& par);
    ~GModelPar();

    // Methods
    std::string name(void) const { return m_name; }
    std::string unit(void) const { return m_unit; }
    double      real_value(void) const { return m_value*m_scale; }
    double      real_error(void) const { return m_error*m_scale; }
    double      real_gradient(void) const { return m_gradient*m_scale; }
    double      real_min(void) const { return m_min*m_scale; }
    double      real_max(void) const { return m_max*m_scale; }
    double      value(void) const { return m_value; }
    double      error(void) const { return m_error; }
    double      gradient(void) const { return m_gradient; }
    double      min(void) const { return m_min; }
    double      max(void) const { return m_max; }
    double      scale(void) const { return m_scale; }
    bool        isfree(void) const { return m_free; }
    bool        hasmin(void) const { return m_hasmin; }
    bool        hasmax(void) const { return m_hasmax; }
    void        name(const std::string& name) { m_name=name; return; }
    void        unit(const std::string& unit) { m_unit=unit; return; }
    void        value(const double& value);
    void        error(const double& error)  { m_error=error; return; }
    void        gradient(const double& gradient) { m_gradient=gradient; return; }
    void        min(const double& min);
    void        max(const double& max);
    void        scale(const double& scale) { m_scale=scale; return; }
    void        range(const double& min, const double& max);
    void        remove_min(void) { m_hasmin=false; return; }
    void        remove_max(void) { m_hasmax=false; return; }
    void        remove_range(void) { m_hasmin=false; m_hasmax=false; return; }
    void        free(void) { m_free=true; return; }
    void        fix(void) { m_free=false; return; }
    void        read(const GXmlElement& xml);
    void        write(GXmlElement& xml) const;
};


/***********************************************************************//**
 * @brief GModelPar class extension
 ***************************************************************************/
%extend GModelPar {
    char *__str__() {
        return tochar(self->print());
    }
    GModelPar copy() {
        return (*self);
    }
};
