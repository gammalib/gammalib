/***************************************************************************
 *        GModelSpectralPlaw.hpp  -  Spectral power law model class        *
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
 * @file GModelSpectralPlaw.hpp
 * @brief GModelSpectralPlaw class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELSPECTRALPLAW_HPP
#define GMODELSPECTRALPLAW_HPP

/* __ Includes ___________________________________________________________ */
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectralPlaw
 *
 * @brief Spectral power law interface definition.
 *
 * This class implements a power law as the spectral component of the
 * gamma-ray data model. The power law is defined as
 * \f[I(E)=norm (E/pivot)^{index}\f]
 * where
 * \f$norm\f$ is the normalization or prefactor,
 * \f$pivot\f$ is the pivot energy, and
 * \f$index\f$ is the spectral index.
 ***************************************************************************/
class GModelSpectralPlaw  : public GModelSpectral {

public:
    // Constructors and destructors
    GModelSpectralPlaw(void);
    explicit GModelSpectralPlaw(const double& norm, const double& index);
    explicit GModelSpectralPlaw(const GXmlElement& xml);
    GModelSpectralPlaw(const GModelSpectralPlaw& model);
    virtual ~GModelSpectralPlaw(void);

    // Operators
    GModelSpectralPlaw& operator= (const GModelSpectralPlaw& model);

    // Implemented pure virtual methods
    void                clear(void);
    GModelSpectralPlaw* clone(void) const;
    int                 size(void) const { return m_npars; }
    std::string         type(void) const { return "PowerLaw"; }
    double              eval(const GEnergy& srcEng);
    double              eval_gradients(const GEnergy& srcEng);
    double              flux(const GEnergy& emin, const GEnergy& emax) const;
    GEnergy             mc(const GEnergy& emin, const GEnergy& emax, GRan& ran) const;
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
    std::string         print(void) const;

    // Other methods
    void   autoscale(void);
    double norm(void) const { return m_norm.real_value(); }
    double index(void) const { return m_index.real_value(); }
    double pivot(void) const { return m_pivot.real_value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralPlaw& model);
    void free_members(void);

    // Implemented pure virtual methods
    GModelPar** par(void) { return m_par; }

    // Protected members
    int        m_npars;           //!< Number of parameters
    GModelPar* m_par[3];          //!< Pointers to parameters
    GModelPar  m_norm;            //!< Normalization factor
    GModelPar  m_index;           //!< Spectral index
    GModelPar  m_pivot;           //!< Pivot energy
};

#endif /* GMODELSPECTRALPLAW_HPP */
