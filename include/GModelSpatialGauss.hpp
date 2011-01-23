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
 * @brief GModelSpatialGauss class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELSPATIALGAUSS_HPP
#define GMODELSPATIALGAUSS_HPP

/* __ Includes ___________________________________________________________ */
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/***********************************************************************//**
 * @class GModelSpatialGauss
 *
 * @brief Gaussian source model interface definition.
 *
 * This class implements the spatial component of the factorised source
 * model for a Gaussian source.
 ***************************************************************************/
class GModelSpatialGauss  : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialGauss(void);
    explicit GModelSpatialGauss(const GSkyDir& dir, const double& sigma);
    explicit GModelSpatialGauss(const GXmlElement& xml);
    GModelSpatialGauss(const GModelSpatialGauss& model);
    virtual ~GModelSpatialGauss(void);

    // Operators
    GModelSpatialGauss& operator= (const GModelSpatialGauss& model);

    // Implemented pure virtual methods
    void                clear(void);
    GModelSpatialGauss* clone(void) const;
    int                 size(void) const { return m_npars; }
    std::string         type(void) const { return "GaussFunction"; }
    double              eval(const GSkyDir& srcDir);
    double              eval_gradients(const GSkyDir& srcDir);
    GSkyDir             mc(GRan& ran) const;
    void                read(const GXmlElement& xml);
    void                write(GXmlElement& xml) const;
    std::string         print(void) const;
    bool                isptsource(void) const { return false; }

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

    // Implemented pure virtual methods
    GModelPar** par(void) { return m_par; }

    // Protected members
    int        m_npars;           //!< Number of parameters
    GModelPar* m_par[3];          //!< Pointers to parameters
    GModelPar  m_ra;              //!< Right Ascension (deg)
    GModelPar  m_dec;             //!< Declination (deg)
    GModelPar  m_sigma;           //!< Gaussian width (deg)
};

#endif /* GMODELSPATIALGAUSS_HPP */
