/***************************************************************************
 *    GModelSpectralExpPlaw.hpp  -  Exponential cut off power law model    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralExpPlaw.hpp
 * @brief GModelSpectralExpPlaw class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELSPECTRALEXPPLAW_HPP
#define GMODELSPECTRALEXPPLAW_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GEnergy.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpectralExpPlaw
 *
 * @brief Spectral power law interface definition.
 *
 * This class implements a power law as the spectral component of the
 * gamma-ray data model. The power law is defined as
 * \f[I(E)=norm (E/pivot)^{index} \exp(-E/ecut)\f]
 * where
 * \f$norm\f$ is the normalization or prefactor,
 * \f$pivot\f$ is the pivot energy,
 * \f$index\f$ is the spectral index, and
 * \f$ecut\f$ is the cut off energy.
 ***************************************************************************/
class GModelSpectralExpPlaw  : public GModelSpectral {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModelSpectralExpPlaw& model);
    friend GLog&         operator<< (GLog& log, const GModelSpectralExpPlaw& model);

public:
    // Constructors and destructors
    explicit GModelSpectralExpPlaw(void);
    explicit GModelSpectralExpPlaw(const double& norm, const double& index, const double& ecut);
    explicit GModelSpectralExpPlaw(const GXmlElement& xml);
    GModelSpectralExpPlaw(const GModelSpectralExpPlaw& model);
    virtual ~GModelSpectralExpPlaw(void);

    // Operators
    GModelSpectralExpPlaw& operator= (const GModelSpectralExpPlaw& model);

    // Implemented pure virtual methods
    void                   clear(void);
    GModelSpectralExpPlaw* clone(void) const;
    int                    size(void) const { return m_npars; }
    std::string            name(void) const { return "ExpCutoff"; }
    GModelPar*             par(int index) const;
    double                 eval(const GEnergy& srcEng);
    double                 eval_gradients(const GEnergy& srcEng);
    void                   read(const GXmlElement& xml);
    void                   write(GXmlElement& xml) const;
    std::string            print(void) const;

    // Other methods
    GModelPar* par_norm(void) { return &m_norm; }
    GModelPar* par_index(void) { return &m_index; }
    GModelPar* par_ecut(void) { return &m_ecut; }
    GModelPar* par_pivot(void) { return &m_pivot; }
    void       autoscale(void);
    double     norm(void) const { return m_norm.real_value(); }
    double     index(void) const { return m_index.real_value(); }
    double     ecut(void) const { return m_ecut.real_value(); }
    double     pivot(void) const { return m_pivot.real_value(); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpectralExpPlaw& model);
    void free_members(void);

    // Data area
    int        m_npars;           //!< Number of parameters
    GModelPar* m_par[4];          //!< Pointers to parameters
    GModelPar  m_norm;            //!< Normalization factor
    GModelPar  m_index;           //!< Spectral index
    GModelPar  m_ecut;            //!< Exponential cut off energy
    GModelPar  m_pivot;           //!< Pivot energy
};

#endif /* GMODELSPECTRALEXPPLAW_HPP */
