/***************************************************************************
 *        GModelSpectral.hpp  -  Abstract spectral model base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectral.hpp
 * @brief GModelSpectral abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMODELSPECTRAL_HPP
#define GMODELSPECTRAL_HPP

/* __ Includes ___________________________________________________________ */
#include "GModelPar.hpp"
#include "GEnergy.hpp"


/***********************************************************************//**
 * @class GModelSpectral
 *
 * @brief Abstract interface definition for the spectral model class.
 *
 * This class implements the spectral component of the factorised source
 * model.
 ***************************************************************************/
class GModelSpectral {

    // Friend classes
    friend class GModel;

public:
    // Constructors and destructors
    GModelSpectral(void);
    GModelSpectral(const GModelSpectral& model);
    virtual ~GModelSpectral();

    // Operators
    virtual GModelSpectral& operator= (const GModelSpectral& model);

    // Virtual methods
    virtual int        npars(void) const = 0;
    virtual GModelPar* par(int index) const = 0;
    virtual double     eval(const GEnergy& srcEng) = 0;
    virtual double     eval_gradients(const GEnergy& srcEng) = 0;

protected:
    // Protected methods
    void                    init_members(void);
    void                    copy_members(const GModelSpectral& model);
    void                    free_members(void);
    virtual GModelSpectral* clone(void) const = 0;
};

#endif /* GMODELSPECTRAL_HPP */
