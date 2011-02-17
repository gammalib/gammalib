/***************************************************************************
 *      GModelRadial.hpp  -  Abstract radial spatial model base class      *
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
 * @file GModelRadial.hpp
 * @brief Abstract radial spatial model base class interface definition
 * @author J. Knodlseder
 */

#ifndef GMODELRADIAL_HPP
#define GMODELRADIAL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"
#include "GRan.hpp"


/***********************************************************************//**
 * @class GModelRadial
 *
 * @brief Abstract radial model base class
 *
 * This class implements the radial spatial component of the factorized
 * gamma-ray source model. Typical examples of radial components are a
 * Gaussian or a shell source.
 ***************************************************************************/
class GModelRadial : public GModelSpatial {

public:
    // Constructors and destructors
    GModelRadial(void);
    GModelRadial(const GModelRadial& model);
    explicit GModelRadial(const GXmlElement& xml);
    virtual ~GModelRadial(void);

    // Operators
    virtual GModelRadial& operator=(const GModelRadial& model);

    // Pure virtual methods
    virtual void          clear(void) = 0;
    virtual GModelRadial* clone(void) const = 0;
    virtual std::string   type(void) const = 0;
    virtual double        eval(const double& theta) const = 0;
    virtual double        eval_gradients(const double& theta) const = 0;
    virtual GSkyDir       mc(GRan& ran) const = 0;
    virtual std::string   print(void) const = 0;

    // Implemented virtual methods
    virtual double        eval(const GSkyDir& srcDir) const;
    virtual double        eval_gradients(const GSkyDir& srcDir) const;
    virtual void          read(const GXmlElement& xml);
    virtual void          write(GXmlElement& xml) const;

    // Other methods
    double                ra(void) const { return m_ra.real_value(); }
    double                dec(void) const { return m_dec.real_value(); }
    GSkyDir               dir(void) const;
    void                  dir(const GSkyDir& dir);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GModelRadial& model);
    void free_members(void);

    // Proteced members
    GModelPar  m_ra;    //!< Right Ascension (deg)
    GModelPar  m_dec;   //!< Declination (deg)
};

#endif /* GMODELRADIAL_HPP */
