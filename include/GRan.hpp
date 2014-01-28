/***************************************************************************
 *                 GRan.hpp - Random number generator class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2014 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GRan.hpp
 * @brief Random number generator class definition
 * @author Juergen Knoedlseder
 */

#ifndef GRAN_HPP
#define GRAN_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GRan
 *
 * @brief Random number generator class
 *
 * This class implements a random number generator.
 ***************************************************************************/
class GRan : public GBase {

public:
    // Constructors and destructors
    GRan(void);
    GRan(unsigned long long int seed);
    GRan(const GRan& ran);
    virtual ~GRan(void);
 
    // Operators
    GRan& operator= (const GRan& ran);

    // Methods
    void                   clear(void);
    GRan*                  clone(void) const;
    void                   seed(unsigned long long int seed);
    unsigned long long int seed(void) const;
    unsigned long int      int32(void);
    unsigned long long int int64(void);
    double                 uniform(void);
    double                 exp(const double& lambda);
    double                 poisson(const double& lambda);
    double                 chisq2(void);
    std::string            print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void                   init_members(unsigned long long int seed = 41L);
    void                   copy_members(const GRan& ran);
    void                   free_members(void);

    // Protected data members
    unsigned long long int m_seed;    //!< Random number generator seed
    unsigned long long int m_value1;  //!< Value 1
    unsigned long long int m_value2;  //!< Value 2
    unsigned long long int m_value3;  //!< Value 3

    // Poisson cache
    double m_old_lambda;    //!< Old lambda value
    double m_sqrt_lambda;   //!< sqrt(2*lambda)
    double m_log_lambda;    //!< log(lambda)
    double m_exp_lambda;    //!< exp(-lambda)
};


/***********************************************************************//**
 * @brief Return seed value
 *
 * @return Seed value.
 ***************************************************************************/
inline
unsigned long long int GRan::seed(void) const
{
    return m_seed;
}

#endif /* GRAN_HPP */
