/***************************************************************************
 *           GCOMD1Response.hpp - COMPTEL D1 module response class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMD1Response.hpp
 * @brief COMPTEL D1 module response class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMD1RESPONSE_HPP
#define GCOMD1RESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GCaldb.hpp"
#include "GNodeArray.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GFitsTable;


/***********************************************************************//**
 * @class GCOMD1Response
 *
 * @brief Interface for the COMPTEL D1 module response class
 ***************************************************************************/
class GCOMD1Response : public GBase {

public:
    // Constructors and destructors
    GCOMD1Response(void);
    GCOMD1Response(const GCOMD1Response& rsp);
    GCOMD1Response(const GCaldb& caldb, const std::string& sdaname);
    ~GCOMD1Response(void);

    // Operators
    GCOMD1Response& operator=(const GCOMD1Response& rsp);
    double          operator()(const double& etrue, const double& ereco) const;

    // Methods
    void            clear(void);
    GCOMD1Response* clone(void) const;
    std::string     classname(void) const;
    void            caldb(const GCaldb& caldb);
    const GCaldb&   caldb(void) const;
    void            load(const std::string& sdaname);
    void            read(const GFitsTable& table);
    double          position(const double& etrue) const;
    double          sigma(const double& etrue) const;
    double          amplitude(const double& etrue) const;
    double          emin(const double& etrue) const;
    double          ewidth(const double& etrue) const;
    double          emax(const double& etrue) const;
    std::string     print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GCOMD1Response& rsp);
    void free_members(void);
    void update_cache(const double& etrue) const;

    // Private data members
    GCaldb              m_caldb;      //!< Calibration database
    GNodeArray          m_energies;   //!< Input energies
    std::vector<double> m_positions;  //!< Photo peak position in MeV
    std::vector<double> m_sigmas;     //!< Photo peak width in MeV
    std::vector<double> m_amplitudes; //!< Photo peak amplitude
    std::vector<double> m_emins;      //!< Lower energy threshold of D1
    std::vector<double> m_ewidths;    //!< Lower energy threshold width of D1
    std::vector<double> m_emaxs;      //!< Upper energy limit of D1

    // Pre-computation cache
    mutable double m_energy;
    mutable double m_position;
    mutable double m_sigma;
    mutable double m_amplitude;
    mutable double m_emin;
    mutable double m_ewidth;
    mutable double m_emax;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMD1Response").
 ***************************************************************************/
inline
std::string GCOMD1Response::classname(void) const
{
    return ("GCOMD1Response");
}


/***********************************************************************//**
 * @brief Return calibration database
 *
 * @return Calibration database.
 ***************************************************************************/
inline
const GCaldb& GCOMD1Response::caldb(void) const
{
    return (m_caldb);
}


/***********************************************************************//**
 * @brief Set calibration database
 *
 * @param[in] caldb Calibration database.
 *
 * Sets the calibration database containing the COMPTEL D1 module response.
 ***************************************************************************/
inline
void GCOMD1Response::caldb(const GCaldb& caldb)
{
    m_caldb = caldb;
    return;
}


/***********************************************************************//**
 * @brief Return photo peak position
 *
 * @param[in] etrue True energy (MeV).
 * @return Photo peak position (MeV).
 ***************************************************************************/
inline
double GCOMD1Response::position(const double& etrue) const
{
    update_cache(etrue);
    return (m_position);
}


/***********************************************************************//**
 * @brief Return photo peak standard deviation
 *
 * @param[in] etrue True energy (MeV).
 * @return Photo peak standard deviation (MeV).
 ***************************************************************************/
inline
double GCOMD1Response::sigma(const double& etrue) const
{
    update_cache(etrue);
    return (m_sigma);
}


/***********************************************************************//**
 * @brief Return photo peak amplitude
 *
 * @param[in] etrue True energy (MeV).
 * @return Photo peak amplitude.
 ***************************************************************************/
inline
double GCOMD1Response::amplitude(const double& etrue) const
{
    update_cache(etrue);
    return (m_amplitude);
}


/***********************************************************************//**
 * @brief Return minimum energy
 *
 * @param[in] etrue True energy (MeV).
 * @return Minimum energy (MeV).
 ***************************************************************************/
inline
double GCOMD1Response::emin(const double& etrue) const
{
    update_cache(etrue);
    return (m_emin);
}


/***********************************************************************//**
 * @brief Return energy threshold width
 *
 * @param[in] etrue True energy (MeV).
 * @return Energy threshold width (MeV).
 ***************************************************************************/
inline
double GCOMD1Response::ewidth(const double& etrue) const
{
    update_cache(etrue);
    return (m_ewidth);
}


/***********************************************************************//**
 * @brief Return maximum energy
 *
 * @param[in] etrue True energy (MeV).
 * @return Maximum energy (MeV).
 ***************************************************************************/
inline
double GCOMD1Response::emax(const double& etrue) const
{
    update_cache(etrue);
    return (m_emax);
}

#endif /* GCOMD1RESPONSE_HPP */
