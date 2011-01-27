/***************************************************************************
 *      GCTAModelRadialGauss.hpp  -  Radial Gaussian CTA model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAModelRadialGauss.hpp
 * @brief GCTAModelRadialGauss class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAMODELRADIALGAUSS_HPP
#define GCTAMODELRADIALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include "GModelPar.hpp"
#include "GXmlElement.hpp"
#include "GCTAModelRadial.hpp"
#include "GCTAInstDir.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialGauss
 *
 * @brief Radial Gaussian CTA model interface definition
 ***************************************************************************/
class GCTAModelRadialGauss  : public GCTAModelRadial {

public:
    // Constructors and destructors
    GCTAModelRadialGauss(void);
    explicit GCTAModelRadialGauss(const double& sigma);
    explicit GCTAModelRadialGauss(const GXmlElement& xml);
    GCTAModelRadialGauss(const GCTAModelRadialGauss& model);
    virtual ~GCTAModelRadialGauss(void);

    // Operators
    GModelPar&            operator() (int index);
    const GModelPar&      operator() (int index) const;
    GCTAModelRadialGauss& operator= (const GCTAModelRadialGauss& model);

    // Implemented pure virtual methods
    void                  clear(void);
    GCTAModelRadialGauss* clone(void) const;
    int                   size(void) const { return m_npars; }
    std::string           type(void) const { return "GaussFunction"; }
    double                eval(const double& offset);
    double                eval_gradients(const double& offset);
    GCTAInstDir           mc(const GCTAInstDir& dir, GRan& ran) const;
    double                omega(void) const;
    void                  read(const GXmlElement& xml);
    void                  write(GXmlElement& xml) const;
    std::string           print(void) const;

    // Other methods
    double sigma(void) const { return m_sigma.real_value(); }
    void   sigma(const double& sigma) { m_sigma.real_value(sigma); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRadialGauss& model);
    void free_members(void);

    // Protected members
    int        m_npars;        //!< Number of parameters
    GModelPar* m_par[1];       //!< Pointers to parameters
    GModelPar  m_sigma;        //!< Gaussian width (deg)
};

#endif /* GCTAMODELRADIALGAUSS_HPP */
