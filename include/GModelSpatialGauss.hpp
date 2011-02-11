/***************************************************************************
 *     GModelSpatialGauss.hpp  -  Spatial Gaussian source model class      *
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
 * @file GModelSpatialGauss.hpp
 * @brief Gaussian spatial model class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIALGAUSS_HPP
#define GMODELSPATIALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialGauss
 *
 * @brief Gaussian spatial model class
 *
 * This class implements the spatial component of the factorised source
 * model for a Gaussian source.
 ***************************************************************************/
class GModelSpatialGauss : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialGauss(void);
    explicit GModelSpatialGauss(const GSkyDir& dir, const double& sigma);
    explicit GModelSpatialGauss(const GXmlElement& xml);
    GModelSpatialGauss(const GModelSpatialGauss& model);
    virtual ~GModelSpatialGauss(void);

    // Operators
    virtual GModelSpatialGauss& operator=(const GModelSpatialGauss& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpatialGauss* clone(void) const;
    virtual std::string         type(void) const { return "GaussFunction"; }
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(void) const;

    // Other methods
    double  ra(void) const { return m_ra.real_value(); }
    double  dec(void) const { return m_dec.real_value(); }
    double  sigma(void) const { return m_sigma.real_value(); }
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
    void    sigma(const double& sigma) { m_sigma.real_value(sigma); }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelSpatialGauss& model);
    void free_members(void);

    // Protected members
    GModelPar m_ra;         //!< Right Ascension (deg)
    GModelPar m_dec;        //!< Declination (deg)
    GModelPar m_sigma;      //!< Gaussian width (deg)
};

#endif /* GMODELSPATIALGAUSS_HPP */
