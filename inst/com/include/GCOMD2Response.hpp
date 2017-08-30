/***************************************************************************
 *           GCOMD2Response.hpp - COMPTEL D2 module response class         *
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
 * @file GCOMD2Response.hpp
 * @brief COMPTEL D2 module response class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMD2RESPONSE_HPP
#define GCOMD2RESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GCaldb.hpp"
#include "GNodeArray.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GFitsTable;


/***********************************************************************//**
 * @class GCOMD2Response
 *
 * @brief Interface for the COMPTEL D2 module response class
 ***************************************************************************/
class GCOMD2Response : public GBase {

public:
    // Constructors and destructors
    GCOMD2Response(void);
    GCOMD2Response(const GCOMD2Response& rsp);
    GCOMD2Response(const GCaldb& caldb, const std::string& sdbname);
    ~GCOMD2Response(void);

    // Operators
    GCOMD2Response& operator=(const GCOMD2Response& rsp);
    double          operator()(const double& etrue, const double& ereco) const;

    // Methods
    void            clear(void);
    GCOMD2Response* clone(void) const;
    std::string     classname(void) const;
    void            caldb(const GCaldb& caldb);
    const GCaldb&   caldb(void) const;
    void            load(const std::string& sdbname);
    void            read(const GFitsTable& table);
    std::string     print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GCOMD2Response& rsp);
    void free_members(void);
    void update_cache(const double& etrue) const;

    // Private data members
    GCaldb              m_caldb;       //!< Calibration database
    GNodeArray          m_energies;    //!< Input energies
    std::vector<double> m_positions;   //!< Photo peak position in MeV
    std::vector<double> m_sigmas;      //!< Photo peak width in MeV
    std::vector<double> m_amplitudes;  //!< Photo peak amplitude
    std::vector<double> m_escapes1;    //!< Amplitude of first escape peak
    std::vector<double> m_escapes2;    //!< Amplitude of second escape peak
    std::vector<double> m_tails;       //!< Amplitude of Compton tail
    std::vector<double> m_backgrounds; //!< Amplitude of Compton background
    std::vector<double> m_emins;       //!< Lower energy threshold of D2
    std::vector<double> m_ewidths;     //!< Lower energy threshold width of D2
    std::vector<double> m_emaxs;       //!< Upper energy limit of D2

    // Pre-computation cache
    mutable double m_energy;
    mutable double m_position;
    mutable double m_sigma;
    mutable double m_amplitude;
    mutable double m_escape1;
    mutable double m_escape2;
    mutable double m_tail;
    mutable double m_background;
    mutable double m_emin;
    mutable double m_ewidth;
    mutable double m_emax;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMD2Response").
 ***************************************************************************/
inline
std::string GCOMD2Response::classname(void) const
{
    return ("GCOMD2Response");
}


/***********************************************************************//**
 * @brief Return calibration database
 *
 * @return Calibration database.
 ***************************************************************************/
inline
const GCaldb& GCOMD2Response::caldb(void) const
{
    return m_caldb;
}


/***********************************************************************//**
 * @brief Set calibration database
 *
 * @param[in] caldb Calibration database.
 *
 * Sets the calibration database containing the COMPTEL D2 module response.
 ***************************************************************************/
inline
void GCOMD2Response::caldb(const GCaldb& caldb)
{
    m_caldb = caldb;
    return;
}

#endif /* GCOMD2RESPONSE_HPP */
