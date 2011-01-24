/***************************************************************************
 *             GModelFactorized.i  -  Model factorization class            *
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
 * @file GModelFactorized.i
 * @brief GModelFactorized class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelFactorized.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModelFactorized
 *
 * @brief Model factorization class interface defintion
 ***************************************************************************/
class GModelFactorized {
public:
    // Constructors and destructors
    GModelFactorized(void);
    GModelFactorized(const GModelFactorized& model);
    virtual ~GModelFactorized(void);

    // Methods
    GModelSpatial*  spatial(void)  const { return m_spatial; }
    GModelSpectral* spectral(void) const { return m_spectral; }
    GModelTemporal* temporal(void) const { return m_temporal; }
};


/***********************************************************************//**
 * @brief GModelFactorized class extension
 ***************************************************************************/
%extend GModelFactorized {
};
