/***************************************************************************
 *                  GModelPar.hpp  -  Model parameter class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelPar.hpp
 * @brief GModelPar class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELPAR_HPP
#define GMODELPAR_HPP

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include <iostream>


/***********************************************************************//**
 * @class GModelPar
 *
 * @brief GModelPar class interface defintion.
 ***************************************************************************/
class GModelPar {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModelPar& par);

public:
    // Constructors and destructors
    GModelPar(void);
    GModelPar(const GModelPar& par);
    ~GModelPar();
 
    // Operators
    GModelPar& operator= (const GModelPar& par);

    // Methods
    std::string name(void) const { return m_name; }
    std::string unit(void) const { return m_unit; }
    double      real_value(void) const { return m_value*m_scale; }
    double      real_error(void) const { return m_error*m_scale; }
    double      real_min(void) const { return m_min*m_scale; }
    double      real_max(void) const { return m_max*m_scale; }
    double      value(void) const { return m_value; }
    double      error(void) const { return m_error; }
    double      min(void) const { return m_min; }
    double      max(void) const { return m_max; }
    double      scale(void) const { return m_scale; }
    bool        free(void) const { return m_free; }
    bool        hasmin(void) const { return m_hasmin; }
    bool        hasmax(void) const { return m_hasmax; }
    void        name(const std::string& name) { m_name=name; return; }
    void        unit(const std::string& unit) { m_unit=unit; return; }
    void        value(const double& value);
    void        error(const double& error);
    void        min(const double& min);
    void        max(const double& max);
    void        scale(const double& scale) { m_scale=scale; return; }
    void        range(const double& min, const double& max);
    void        remove_min(void) { m_hasmin=false; return; }
    void        remove_max(void) { m_hasmax=false; return; }
    void        remove_range(void) { m_hasmin=false; m_hasmax=false; return; }
    void        free(void) { m_free=true; return; }
    void        fix(void) { m_free=false; return; }
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelPar& par);
    void free_members(void);

    // Proteced data members
    std::string  m_name;         //!< Parameter name
    std::string  m_unit;         //!< Parameter unit
    double       m_value;        //!< Parameter value
    double       m_error;        //!< Uncertainty in parameter value
    double       m_min;          //!< Parameter minimum
    double       m_max;          //!< Parameter maximum
    double       m_scale;        //!< Parameter scale (real = m_value * m_scale)
    bool         m_free;         //!< Parameter is free
    bool         m_hasmin;       //!< Parameter has minimum boundary
    bool         m_hasmax;       //!< Parameter has maximum boundary
};

#endif /* GMODELPAR_HPP */
