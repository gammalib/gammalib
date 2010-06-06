/***************************************************************************
 *       GIntegrand.hpp  -  Integrand function abstract base class         *
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
 * @file GIntegrand.hpp
 * @brief GIntegrand abstract virtual base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GINTEGRAND_HPP
#define GINTEGRAND_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>


/***********************************************************************//**
 * @class GIntegrand
 *
 * @brief GIntegrand class interface defintion.
 ***************************************************************************/
class GIntegrand {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GIntegrand& arg);

public:

    // Constructors and destructors
    GIntegrand(void);
    GIntegrand(const GIntegrand& arg);
    virtual ~GIntegrand(void);

    // Operators
    GIntegrand& operator= (const GIntegrand& arg);

    // Methods
    virtual double eval(double x) = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GIntegrand& arg);
    void free_members(void);
    
    // Protected data area

};

#endif /* GINTEGRAND_HPP */
