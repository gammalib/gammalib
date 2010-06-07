/***************************************************************************
 *                test_GNumerics.hpp  -  test numerics modules             *
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

#ifndef TEST_NUMERICS_HPP
#define TEST_NUMERICS_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @class Gauss
 *
 * @brief Gaussian function.
 ***************************************************************************/
class Gauss : public GIntegrand {
public:
    Gauss(const double& sigma) : m_sigma(sigma) { return; }
    double eval(double x) {
        double arg = -0.5*x*x/m_sigma/m_sigma;
        double val = 1.0/sqrt(twopi)/m_sigma * exp(arg);
        return val;
    }
protected:
    double m_sigma;
};

#endif /* TEST_NUMERICS_HPP */
