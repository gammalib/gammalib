/***************************************************************************
 *        GModelSpectralConst.hpp  -  Spectral constant model class        *
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
 * @file GModelSpectralConst.hpp
 * @brief GModelSpectralConst class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELSPECTRALCONST_HPP
#define GMODELSPECTRALCONST_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectralConst
 *
 * @brief Spectral constant interface definition
 *
 * This class implements a constant as the spectral component of the
 * gamma-ray data model. The function is defined as
 * \f[I(E)=norm\f]
 * where
 * \f$norm\f$ is the normalization constant.
 ***************************************************************************/
class GModelSpectralConst  : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralConst(void);
    explicit GModelSpectralConst(const GXmlElement& xml);
    GModelSpectralConst(const GModelSpectralConst& model);
    virtual ~GModelSpectralConst(void);

    // Operators
    GModelSpectralConst& operator= (const GModelSpectralConst& model);

    // Implemented pure virtual methods
    void                 clear(void);
    GModelSpectralConst* clone(void) const;
    int                  size(void) const { return m_npars; }
    std::string          type(void) const { return "ConstantValue"; }
    double               eval(const GEnergy& srcEng);
    double               eval_gradients(const GEnergy& srcEng);
    double               flux(const GEnergy& emin, const GEnergy& emax) const;
    GEnergy              mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    void                 read(const GXmlElement& xml);
    void                 write(GXmlElement& xml) const;
    std::string          print(void) const;

    // Other methods
    double norm(void) const { return m_norm.real_value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralConst& model);
    void free_members(void);

    // Implemented pure virtual methods
    GModelPar** par(void) { return m_par; }

    // Protected members
    int        m_npars;     //!< Number of parameters
    GModelPar* m_par[1];    //!< Pointers to parameters
    GModelPar  m_norm;      //!< Normalization factor
};

#endif /* GMODELSPECTRALCONST_HPP */
