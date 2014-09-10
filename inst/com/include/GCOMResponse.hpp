/***************************************************************************
 *                 GCOMResponse.hpp - COMPTEL Response class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @file GCOMResponse.hpp
 * @brief COMPTEL instrument response function class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMRESPONSE_HPP
#define GCOMRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvent.hpp"
#include "GPhoton.hpp"
#include "GObservation.hpp"
#include "GResponse.hpp"
#include "GFitsImage.hpp"
#include "GCaldb.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declaration ________________________________________________ */


/***********************************************************************//**
 * @class GCOMResponse
 *
 * @brief Interface for the COMPTEL instrument response function
 ***************************************************************************/
class GCOMResponse : public GResponse {

public:
    // Constructors and destructors
    GCOMResponse(void);
    GCOMResponse(const GCOMResponse& rsp);
    GCOMResponse(const std::string& rspname, const GCaldb& caldb);
    virtual ~GCOMResponse(void);

    // Operators
    virtual GCOMResponse& operator=(const GCOMResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GCOMResponse* clone(void) const;
    virtual std::string   classname(void) const;
    virtual bool          use_edisp(void) const;
    virtual bool          use_tdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

    // Other Methods
    void               caldb(const GCaldb& caldb);
    const GCaldb&      caldb(void) const;
    const std::string& rspname(void) const;
    void               load(const std::string& rspname);
    void               read(const GFitsImage& hdu);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GCOMResponse& rsp);
    void free_members(void);

    // Private data members
    GCaldb              m_caldb;             //!< Name of or path to the calibration database
    std::string         m_rspname;           //!< Response name
    std::vector<double> m_iaq;               //!< IAQ array
    int                 m_phigeo_bins;       //!< Number of Phigeo bins
    int                 m_phibar_bins;       //!< Number of Phibar bins
    double              m_phigeo_ref_value;  //!< Phigeo reference value (deg)
    double              m_phigeo_ref_pixel;  //!< Phigeo reference pixel (starting from 1)
    double              m_phigeo_bin_size;   //!< Phigeo binsize (deg)
    double              m_phigeo_min;        //!< Phigeo value of first bin (deg)
    double              m_phibar_ref_value;  //!< Phigeo reference value (deg)
    double              m_phibar_ref_pixel;  //!< Phigeo reference pixel (starting from 1)
    double              m_phibar_bin_size;   //!< Phigeo binsize (deg)
    double              m_phibar_min;        //!< Phigeo value of first bin (deg)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMResponse").
 ***************************************************************************/
inline
std::string GCOMResponse::classname(void) const
{
    return ("GCOMResponse");
}


/***********************************************************************//**
 * @brief Signal if energy dispersion will be used
 *
 * @return False.
 ***************************************************************************/
inline
bool GCOMResponse::use_edisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal if time dispersion will be used
 *
 * @return False.
 ***************************************************************************/
inline
bool GCOMResponse::use_tdisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Return calibration database
 *
 * @return Calibration database.
 ***************************************************************************/
inline
const GCaldb& GCOMResponse::caldb(void) const
{
    return m_caldb;
}


/***********************************************************************//**
 * @brief Set path to the calibration database
 *
 * @param[in] caldb Calibration database.
 *
 * Sets the calibration database for the COMPTEL response.
 ***************************************************************************/
inline
void GCOMResponse::caldb(const GCaldb& caldb)
{
    m_caldb = caldb;
    return;
}


/***********************************************************************//**
 * @brief Return response name
 *
 * @return Response name.
 ***************************************************************************/
inline
const std::string& GCOMResponse::rspname(void) const
{
    // Return response name
    return m_rspname;
}

#endif /* GCOMRESPONSE_HPP */
