/***************************************************************************
 *       GLATEfficiency.hpp - Fermi/LAT IRF efficiency factor functor      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2018 by Juergen Knoedlseder                         *
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
 * @file GLATEfficiency.hpp
 * @brief Fermi/LAT IRF efficiency factor functor class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATEFFICIENCY_HPP
#define GLATEFFICIENCY_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GLATEfficiency
 *
 * @brief Interface for the Fermi/LAT efficiency factor functor
 ***************************************************************************/
class GLATEfficiency : public GBase {

public:
    // Constructors and destructors
    GLATEfficiency(void);
    explicit GLATEfficiency(const std::vector<double>& pars);
    GLATEfficiency(const GLATEfficiency& eff);
    virtual ~GLATEfficiency(void);

    // Operators
    GLATEfficiency& operator=(const GLATEfficiency& eff);
    double          operator()(const double& logE) const;

    // Methods
    void                clear(void);
    GLATEfficiency*     clone(void) const;
    std::string         classname(void) const;
    std::vector<double> pars(void) const;
    std::string         print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATEfficiency& eff);
    void free_members(void);
    
    // Protected members
    double m_a0;        //!< Energy domain 1 scale
    double m_a1;        //!< Energy domain 2 scale
    double m_a2;        //!< Energy domain 3 scale
    double m_b0;        //!< Energy domain 1 offset
    double m_b1;        //!< Energy domain 2 offset
    double m_b2;        //!< Energy domain 3 offset
    double m_logEb1;    //!< Separation of energy domains 1/2
    double m_logEb2;    //!< Separation of energy domains 2/3
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATEfficiency").
 ***************************************************************************/
inline
std::string GLATEfficiency::classname(void) const
{
    return ("GLATEfficiency");
}

#endif /* GLATEFFICIENCY_HPP */
