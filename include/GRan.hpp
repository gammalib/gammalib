/***************************************************************************
 *                 GRan.hpp - Randon number generator class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2012 by Juergen Knoedlseder                         *
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
 * @brief Randon number generator class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GRAN_HPP
#define GRAN_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GBase.hpp"
#include "GLog.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GRan
 *
 * @brief Random number generator class
 *
 * This class implements a random number generator.
 ***************************************************************************/
class GRan : GBase {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GRan& ran);
    friend GLog&         operator<< (GLog& log, const GRan& ran);

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
    unsigned long int      int32(void);
    unsigned long long int int64(void);
    double                 uniform(void);
    double                 exp(const double& lambda);
    double                 chisq2(void);
    std::string            print(void) const;
  
protected:
    // Protected methods
    void                   init_members(unsigned long long int seed = 41L);
    void                   copy_members(const GRan& ran);
    void                   free_members(void);

    // Protected data members
    unsigned long long int m_seed;  //!< Random number generator seed
    unsigned long long int m_u;     //!< u
    unsigned long long int m_v;     //!< v
    unsigned long long int m_w;     //!< w
};

#endif /* GRAN_HPP */
