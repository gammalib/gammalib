/***************************************************************************
 *            GOptimizerLM.hpp  -  Levenberg Marquardt optimizer           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GOptimizerLM.hpp
 * @brief GOptimizerLM base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GOPTIMIZERLM_HPP
#define GOPTIMIZERLM_HPP

/* __ Includes ___________________________________________________________ */
#include "GOptimizer.hpp"


/***********************************************************************//**
 * @class GOptimizerLM
 *
 * @brief GOptimizerLM class interface defintion.
 ***************************************************************************/
class GOptimizerLM : public GOptimizer {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GOptimizerLM& opt);

public:

    // Constructors and destructors
    GOptimizerLM();
    GOptimizerLM(const GOptimizerFunction& fct, const GOptimizerPars &pars);
    GOptimizerLM(const GOptimizerFunction& fct, const GModels &models);
    GOptimizerLM(const GOptimizerLM& opt);
    virtual ~GOptimizerLM();

    // Operators
    GOptimizerLM& operator= (const GOptimizerLM& opt);

    // Methods
    void optimize(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerLM& opt);
    void free_members(void);
    
    // Protected data area

};

#endif /* GOPTIMIZERLM_HPP */
