/***************************************************************************
 *      GModelRadialGauss.hpp  -  Radial Gaussian source model class       *
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
 * @file GModelRadialGauss.hpp
 * @brief Radial Gaussian model class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELRADIALGAUSS_HPP
#define GMODELRADIALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelRadial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelRadialGauss
 *
 * @brief Radial Gaussian model class
 *
 * This class implements the spatial component of the factorised source
 * model for a Gaussian source.
 ***************************************************************************/
class GModelRadialGauss : public GModelRadial {

public:
    // Constructors and destructors
    GModelRadialGauss(void);
    explicit GModelRadialGauss(const GSkyDir& dir, const double& sigma);
    explicit GModelRadialGauss(const GXmlElement& xml);
    GModelRadialGauss(const GModelRadialGauss& model);
    virtual ~GModelRadialGauss(void);

    // Operators
    virtual GModelRadialGauss& operator=(const GModelRadialGauss& model);

    // Implemented pure virtual methods
    virtual void               clear(void);
    virtual GModelRadialGauss* clone(void) const;
    virtual std::string        type(void) const { return "GaussFunction"; }
    virtual double             eval(const double& theta) const;
    virtual double             eval_gradients(const double& theta) const;
    virtual GSkyDir            mc(GRan& ran) const;
    virtual double             theta_max(void) const;
    virtual void               read(const GXmlElement& xml);
    virtual void               write(GXmlElement& xml) const;
    virtual std::string        print(void) const;

    // Other methods
    double  sigma(void) const { return m_sigma.real_value(); }
    void    sigma(const double& sigma) { m_sigma.real_value(sigma); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelRadialGauss& model);
    void free_members(void);

    // Protected members
    GModelPar m_sigma;      //!< Gaussian width (deg)
};

#endif /* GMODELRADIALGAUSS_HPP */
