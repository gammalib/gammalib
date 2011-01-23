/***************************************************************************
 *        GModelSpectralFunc.hpp  -  Spectral function model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralFunc.hpp
 * @brief Spectral function model class definition
 * @author J. Knodlseder
 */

#ifndef GMODELSPECTRALFUNC_HPP
#define GMODELSPECTRALFUNC_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"
#include "GNodeArray.hpp"


/***********************************************************************//**
 * @class GModelSpectralFunc
 *
 * @brief Spectral function interface definition
 *
 * This class implements an arbitrary function as the spectral component of
 * the model. The function is defined as
 * \f[I(E)=norm f(E)\f]
 * where
 * \f$norm\f$ is the normalization of the function.
 ***************************************************************************/
class GModelSpectralFunc  : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralFunc(void);
    explicit GModelSpectralFunc(const std::string& filename);
    explicit GModelSpectralFunc(const GXmlElement& xml);
    GModelSpectralFunc(const GModelSpectralFunc& model);
    virtual ~GModelSpectralFunc(void);

    // Operators
    GModelSpectralFunc& operator= (const GModelSpectralFunc& model);

    // Implemented pure virtual methods
    void                clear(void);
    GModelSpectralFunc* clone(void) const;
    int                 size(void) const { return m_npars; }
    std::string         type(void) const { return "FileFunction"; }
    double              eval(const GEnergy& srcEng);
    double              eval_gradients(const GEnergy& srcEng);
    double              flux(const GEnergy& emin, const GEnergy& emax) const;
    GEnergy             mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
    std::string         print(void) const;

    // Other methods
    double norm(void) const { return m_norm.real_value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralFunc& model);
    void free_members(void);
    void load_nodes(const std::string& filename);

    // Implemented pure virtual methods
    GModelPar** par(void) { return m_par; }

    // Protected members
    int                 m_npars;     //!< Number of parameters
    GModelPar*          m_par[1];    //!< Pointers to parameters
    GModelPar           m_norm;      //!< Normalization factor
    GNodeArray          m_nodes;     //!< Nodes of function
    std::vector<double> m_values;    //!< Function values at nodes
    std::string         m_filename;  //!< Name of file function
};

#endif /* GMODELSPECTRALFUNC_HPP */
