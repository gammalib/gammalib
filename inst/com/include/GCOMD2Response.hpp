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
    double          position(const double& etrue) const;
    double          sigma(const double& etrue) const;
    double          amplitude(const double& etrue) const;
    double          escape1(const double& etrue) const;
    double          escape2(const double& etrue) const;
    double          comptontail(const double& etrue) const;
    double          background(const double& etrue) const;
    double          emin(const double& etrue) const;
    double          ewidth(const double& etrue) const;
    double          emax(const double& etrue) const;
    double          emin(void) const;
    double          emax(void) const;
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
    mutable double m_energy;       //!< Incident total energy (MeV)
    mutable double m_position;     //!< Position of photo peak (MeV)
    mutable double m_sigma;        //!< Width of photo peak (MeV)
    mutable double m_amplitude;    //!< Amplitude of photo peak
    mutable double m_escape1;      //!< Amplitude of first escape peak
    mutable double m_escape2;      //!< Amplitude of second escape peak
    mutable double m_tail;         //!< Amplitude of Compton tail
    mutable double m_background;   //!> Amplitude of Compton background
    mutable double m_emin;         //!< Lower energy threshold of D2 (MeV)
    mutable double m_ewidth;       //!< Lower energy threshold width of D2 (MeV)
    mutable double m_emax;         //!< Upper energy limit of D2 (MeV)
    mutable double m_pos_escape1;  //!< Position of first escape peak (MeV)
    mutable double m_pos_escape2;  //!< Position of second escape peak (MeV)
    mutable double m_compton_edge; //!< Position of Compton edge (MeV)
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


/***********************************************************************//**
 * @brief Return photo peak position
 *
 * @param[in] etrue True energy (MeV).
 * @return Photo peak position (MeV).
 ***************************************************************************/
inline
double GCOMD2Response::position(const double& etrue) const
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
double GCOMD2Response::sigma(const double& etrue) const
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
double GCOMD2Response::amplitude(const double& etrue) const
{
    update_cache(etrue);
    return (m_amplitude);
}


/***********************************************************************//**
 * @brief Return first escape peak amplitude
 *
 * @param[in] etrue True energy (MeV).
 * @return First escape peak amplitude.
 ***************************************************************************/
inline
double GCOMD2Response::escape1(const double& etrue) const
{
    update_cache(etrue);
    return (m_escape1);
}


/***********************************************************************//**
 * @brief Return second escape peak amplitude
 *
 * @param[in] etrue True energy (MeV).
 * @return Second escape peak amplitude.
 ***************************************************************************/
inline
double GCOMD2Response::escape2(const double& etrue) const
{
    update_cache(etrue);
    return (m_escape2);
}


/***********************************************************************//**
 * @brief Return Compton tail amplitude
 *
 * @param[in] etrue True energy (MeV).
 * @return Compton tail amplitude.
 ***************************************************************************/
inline
double GCOMD2Response::comptontail(const double& etrue) const
{
    update_cache(etrue);
    return (m_tail);
}


/***********************************************************************//**
 * @brief Return background amplitude
 *
 * @param[in] etrue True energy (MeV).
 * @return Background amplitude.
 ***************************************************************************/
inline
double GCOMD2Response::background(const double& etrue) const
{
    update_cache(etrue);
    return (m_background);
}


/***********************************************************************//**
 * @brief Return minimum energy
 *
 * @param[in] etrue True energy (MeV).
 * @return Minimum energy (MeV).
 ***************************************************************************/
inline
double GCOMD2Response::emin(const double& etrue) const
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
double GCOMD2Response::ewidth(const double& etrue) const
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
double GCOMD2Response::emax(const double& etrue) const
{
    update_cache(etrue);
    return (m_emax);
}


/***********************************************************************//**
 * @brief Return minimum D2 input energy (MeV)
 *
 * @return Minimum energy D2 input energy (MeV).
 *
 * Returns the minimum D2 input energy (MeV). In case that no information
 * has been read from a SDB file so far, the method returns 0.
 ***************************************************************************/
inline
double GCOMD2Response::emin(void) const
{
    double emin = (m_energies.size() > 0) ? m_energies[0] : 0.0;
    return (emin);
}


/***********************************************************************//**
 * @brief Return maximum D2 input energy (MeV)
 *
 * @return Maximum energy D2 input energy (MeV).
 *
 * Returns the maximum D2 input energy (MeV). In case that no information
 * has been read from a SDB file so far, the method returns 0.
 ***************************************************************************/
inline
double GCOMD2Response::emax(void) const
{
    double emax = (m_energies.size() > 0) ? m_energies[m_energies.size()-1] : 0.0;
    return (emax);
}

#endif /* GCOMD2RESPONSE_HPP */
