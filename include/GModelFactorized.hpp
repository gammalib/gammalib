/***************************************************************************
 *            GModelFactorized.hpp  -  Model factorization class           *
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
 * @file GModelFactorized.hpp
 * @brief GModelFactorized class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELFACTORIZED_HPP
#define GMODELFACTORIZED_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GModelSpatial.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GXmlElement.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GModelFactorized
 *
 * @brief Model factorization class interface defintion
 *
 * A factorized model is composed of a spatial, a spectral and a temporal
 * model component. The spatial component describes the model prediction
 * as function of the position on the sky (GSkyDir), the spectral component
 * as function of energy (GEnergy), and the temporal component as function
 * of time (GTime). This class provides support methods to any derived class
 * that requires such a factorization.
 ***************************************************************************/
class GModelFactorized {

public:
    // Constructors and destructors
    GModelFactorized(void);
    GModelFactorized(const GModelFactorized& model);
    virtual ~GModelFactorized(void);

    // Operators
    GModelFactorized& operator= (const GModelFactorized& model);

    // Methods
    GModelSpatial*  spatial(void)  const { return m_spatial; }
    GModelSpectral* spectral(void) const { return m_spectral; }
    GModelTemporal* temporal(void) const { return m_temporal; }
    
protected:
    // Protected methods
    void                    init_members(void);
    void                    copy_members(const GModelFactorized& model);
    void                    free_members(void);
    std::vector<GModelPar*> set_par_pointers(void);
    void                    xml_read(const GXmlElement& xml);
    void                    xml_write(GXmlElement& xml, const std::string& name,
                                      const std::string& instruments) const;
    GModelSpatial*          xml_spatial(const GXmlElement& spatial) const;
    GModelSpectral*         xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal*         xml_temporal(const GXmlElement& temporal) const;
    bool                    valid_model(void) const;
    std::string             print_model(void) const;

    // Proteced data members
    GModelSpatial*  m_spatial;       //!< Spatial model
    GModelSpectral* m_spectral;      //!< Spectral model
    GModelTemporal* m_temporal;      //!< Temporal model
};

#endif /* GMODELFACTORIZED_HPP */
