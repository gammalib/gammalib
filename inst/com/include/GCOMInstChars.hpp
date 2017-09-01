/***************************************************************************
 *       GCOMInstChars.hpp - COMPTEL Instrument Characteristics class      *
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
 * @file GCOMInstChars.hpp
 * @brief COMPTEL Instrument Characteristics interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMINSTCHARS_HPP
#define GCOMINSTCHARS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GCaldb.hpp"
#include "GNodeArray.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GFitsTable;


/***********************************************************************//**
 * @class GCOMInstChars
 *
 * @brief Interface for the COMPTEL Instrument Characteristics class
 ***************************************************************************/
class GCOMInstChars : public GBase {

public:
    // Constructors and destructors
    GCOMInstChars(void);
    GCOMInstChars(const GCOMInstChars& ict);
    GCOMInstChars(const GCaldb& caldb, const std::string& ictname);
    ~GCOMInstChars(void);

    // Operators
    GCOMInstChars& operator=(const GCOMInstChars & ict);

    // Methods
    void           clear(void);
    GCOMInstChars* clone(void) const;
    std::string    classname(void) const;
    void           caldb(const GCaldb& caldb);
    const GCaldb&  caldb(void) const;
    void           load(const std::string& ictname);
    double         prob_D1inter(const double& energy) const;
    double         prob_no_multihit(const double& energy) const;
    double         trans_D1(const double& energy) const;
    double         trans_V1(const double& energy) const;
    double         prob_D2inter(const double& energy, const double& phigeo) const;
    double         trans_D2(const double& energy, const double& phigeo) const;
    double         trans_V23(const double& energy, const double& phigeo) const;
    double         prob_no_selfveto(const double& energy, const double& zenith) const;
    double         multi_scatter(const double& energy, const double& phigeo) const;
    double         psd_correction(const double& energy, const double& phigeo) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void   init_members(void);
    void   copy_members(const GCOMInstChars& ict);
    void   free_members(void);
    void   read_coeffs(const GFitsTable& table, GNodeArray& energies,
                       std::vector<double>& coeffs);
    void   read_pos(const GFitsTable& table, std::vector<double>& x,
                    std::vector<double>& y);
    void   read_selfveto(const GFitsTable& table);
    double ne213a_mfpath(const double& energy) const;
    double kn_cross_section(const double& k) const;
    double min_coeff(const std::vector<double>& coeffs) const;
    double max_coeff(const std::vector<double>& coeffs) const;

    // Private data members
    GCaldb              m_caldb;             //!< Calibration database
    GNodeArray          m_d1inter_energies;  //!< D1 interaction coefficient energies (MeV)
    std::vector<double> m_d1inter_coeffs;    //!< D1 interaction coefficients
    GNodeArray          m_d2inter_energies;  //!< D2 interaction coefficient energies (MeV)
    std::vector<double> m_d2inter_coeffs;    //!< D2 interaction coefficients
    GNodeArray          m_alu_energies;      //!< Al interaction coefficient energies (MeV)
    std::vector<double> m_alu_coeffs;        //!< Al interaction coefficients
    GNodeArray          m_aboved1_energies;  //!< Above D1 attenuation coefficient energies (MeV)
    std::vector<double> m_aboved1_coeffs;    //!< Above D1 attenuation coefficients
    GNodeArray          m_veto_energies;     //!< Veto dome attenuation coefficient energies (MeV)
    std::vector<double> m_veto_coeffs;       //!< Veto dome attenuation coefficients
    GNodeArray          m_selfveto_energies; //!< Selfveto energies (MeV)
    GNodeArray          m_selfveto_zeniths;  //!< Selfveto zenith angle (deg)
    std::vector<double> m_selfveto_coeffs;   //!< Selfveto coefficients (probability)
    GNodeArray          m_d1multi_energies;  //!< D1 multihit attenuation coefficient energies (MeV)
    std::vector<double> m_d1multi_coeffs;    //!< D1 multihit attenuation coefficients (probability)
    GNodeArray          m_d2multi_energies;  //!< D2 multihit attenuation coefficient energies (MeV)
    std::vector<double> m_d2multi_coeffs;    //!< D2 multihit attenuation coefficients (probability)
    std::vector<double> m_d1pos_x;           //!< D1 x-position (cm)
    std::vector<double> m_d1pos_y;           //!< D1 y-position (cm)
    std::vector<double> m_d2pos_x;           //!< D2 x-position (cm)
    std::vector<double> m_d2pos_y;           //!< D2 y-position (cm)
    GNodeArray          m_mfpath_energies;   //!< NE213A mean free path energies (MeV)
    std::vector<double> m_mfpath_coeffs;     //!< NE213A mean free path values
    double              m_d1dens;            //!< D1 density (g/cm^-3)
    double              m_d1rad;             //!< D1 radius (cm)
    double              m_d1thick;           //!< D1 thickness (cm)
    double              m_d2dens;            //!< D2 density (g/cm^-3)
    double              m_d2rad;             //!< D2 radius (cm)
    double              m_d2thick;           //!< D2 thickness (cm)
    double              m_thbar;             //!< Average D2 incident angle (deg)
    double              m_delz;              //!< Distance between D1 and D2 levels (cm)
    double              m_aldens;            //!< Density of aluminium plate above D2 (g/cm^-3)
    double              m_althick;           //!< Thickness of aluminium plate above D2 (cm)
    double              m_abdens;            //!< Density above D1 (g/cm^-3)
    double              m_abthick;           //!< Thickness above D1 (cm)
    double              m_vetodens;          //!< Density of veto domes (g/cm^-3)
    double              m_v1thick;           //!< Thickness of V1 veto dome (cm)
    double              m_vthick;            //!< Thickness of V2 and V3 veto domes together (cm)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMInstChars").
 ***************************************************************************/
inline
std::string GCOMInstChars::classname(void) const
{
    return ("GCOMInstChars");
}


/***********************************************************************//**
 * @brief Return calibration database
 *
 * @return Calibration database.
 ***************************************************************************/
inline
const GCaldb& GCOMInstChars::caldb(void) const
{
    return m_caldb;
}


/***********************************************************************//**
 * @brief Set calibration database
 *
 * @param[in] caldb Calibration database.
 *
 * Sets the calibration database containing the COMPTEL D1 module response.
 ***************************************************************************/
inline
void GCOMInstChars::caldb(const GCaldb& caldb)
{
    m_caldb = caldb;
    return;
}

#endif /* GCOMINSTCHARS_HPP */
