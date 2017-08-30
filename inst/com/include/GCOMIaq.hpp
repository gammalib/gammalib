/***************************************************************************
 *         GCOMIaq.hpp - COMPTEL instrument response representation        *
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
 * @file GCOMIaq.hpp
 * @brief COMPTEL instrument response representation class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMIAQ_HPP
#define GCOMIAQ_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GMath.hpp"
#include "GEbounds.hpp"
#include "GFunction.hpp"
#include "GFitsImageFloat.hpp"
#include "GCOMD1Response.hpp"
#include "GCOMD2Response.hpp"
#include "GCOMInstChars.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GFilename;


/***********************************************************************//**
 * @class GCOMIaq
 *
 * @brief Interface for the COMPTEL instrument response representation class
 ***************************************************************************/
class GCOMIaq : public GBase {

public:
    // Constructors and destructors
    GCOMIaq(void);
    GCOMIaq(const GCOMIaq& iaq);
    GCOMIaq(const double&   phigeo_max, const double& phigeo_bin_size,
            const double&   phibar_max, const double& phibar_bin_size,
            const GEbounds& ebounds);
    ~GCOMIaq(void);

    // Operators
    GCOMIaq& operator=(const GCOMIaq & iaq);

    // Methods
    void        clear(void);
    GCOMIaq*    clone(void) const;
    std::string classname(void) const;
    void        save(const GFilename& filename, const bool& clobber) const;
    std::string print(const GChatter& chatter = NORMAL) const;
    void        iaqwei(const GEnergy& energy, const double& weight);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GCOMIaq& iaq);
    void free_members(void);

    // RESPSI methods
    //void   iaqwei(const GEnergy& energy, const double& weight);
    double respsc(const double& etrue1, const double& etrue2, const double& phibar);

    // Reponse integration kernel
    class response_kernel : public GFunction {
    public:
        response_kernel(const GCOMD1Response& response_d1,
                        const GCOMD2Response& response_d2,
                        const double&         etrue1,
                        const double&         etrue2,
                        const double&         phibar,
                        const double&         etmin,
                        const double&         etmax,
                        const double&         e1min,
                        const double&         e1max,
                        const double&         e2min,
                        const double&         e2max) :
                        m_response_d1(response_d1),
                        m_response_d2(response_d2),
                        m_etrue1(etrue1),
                        m_etrue2(etrue2),
                        m_cos_phibar(std::cos(phibar*gammalib::deg2rad)),
                        m_sin_phibar(std::sin(phibar*gammalib::deg2rad)),
                        m_etmin(etmin),
                        m_etmax(etmax),
                        m_e1min(e1min),
                        m_e1max(e1max),
                        m_e2min(e2min),
                        m_e2max(e2max) {}
        double eval(const double& energy1);
    protected:
        const GCOMD1Response& m_response_d1; //!< Reference to D1 module response
        const GCOMD2Response& m_response_d2; //!< Reference to D2 module response
        double                m_etrue1;      //!< True D1 energy (MeV)
        double                m_etrue2;      //!< True D2 energy (MeV)
        double                m_cos_phibar;  //!< cos(phibar)
        double                m_sin_phibar;  //!< sin(phibar)
        double                m_etmin;       //!< Minimum total energy (MeV)
        double                m_etmax;       //!< Maximum total energy (MeV)
        double                m_e1min;       //!< Minimum D1 energy (MeV)
        double                m_e1max;       //!< Maximum D1 energy (MeV)
        double                m_e2min;       //!< Minimum D2 energy (MeV)
        double                m_e2max;       //!< Maximum D2 energy (MeV)
    };

    // Private data members
    GFitsImageFloat m_iaq;             //!< Response
    GEbounds        m_ebounds;         //!< Energy boundaries
    GCOMD1Response  m_response_d1;     //!< D1 module response
    GCOMD2Response  m_response_d2;     //!< D2 module response
    GCOMInstChars   m_ict;             //!< Instrument characteristics
    double          m_phigeo_max;      //!< Maximum geometrical scatter angle (deg)
    double          m_phibar_max;      //!< Maximum Compton scatter angle (deg)
    double          m_phigeo_bin_size; //!< Bin size in geometrical scatter angle (deg)
    double          m_phibar_bin_size; //!< Bin size in Compton scatter angle (deg)
    double          m_e1min;           //!< Minimum D1 energy (MeV)
    double          m_e1max;           //!< Maximum D1 energy (MeV)
    double          m_e2min;           //!< Minimum D2 energy (MeV)
    double          m_e2max;           //!< Maximum D2 energy (MeV)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMIaq").
 ***************************************************************************/
inline
std::string GCOMIaq::classname(void) const
{
    return ("GCOMIaq");
}

#endif /* GCOMIAQ_HPP */
