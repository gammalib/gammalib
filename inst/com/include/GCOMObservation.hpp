/***************************************************************************
 *           GCOMObservation.hpp  -  COMPTEL observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCOMObservation.hpp
 * @brief COMPTEL observation class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMOBSERVATION_HPP
#define GCOMOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GSkymap.hpp"
#include "GCOMResponse.hpp"
#include "GCOMPointing.hpp"


/***********************************************************************//**
 * @class GCOMObservation
 *
 * @brief Interface class for COMPTEL observations
 *
 * This class implements a multi-wavelength observation. A multi-wavelength
 * observation contains spectral points obtained with an unspecified
 * instrument. The spectral points are given in physical units.
 ***************************************************************************/
class GCOMObservation : public GObservation {

public:
    // Constructors and destructors
    GCOMObservation(void);
    explicit GCOMObservation(const std::string& drename,
                             const std::string& drbname,
                             const std::string& drgname,
                             const std::string& drxname);
    GCOMObservation(const GCOMObservation& obs);
    virtual ~GCOMObservation(void);

    // Operators
    virtual GCOMObservation& operator= (const GCOMObservation& obs);

    // Implement pure virtual methods
    virtual void             clear(void);
    virtual GCOMObservation* clone(void) const;
    virtual void             response(const GResponse& rsp);
    virtual GCOMResponse*    response(void) const;
    virtual GCOMPointing*    pointing(void) const;
    virtual std::string      instrument(void) const { return m_instrument; }
    virtual double           ontime(void) const { return m_ontime; }
    virtual double           livetime(void) const { return m_livetime; }
    virtual double           deadc(const GTime& time) const { return m_deadc; }
    virtual void             read(const GXmlElement& xml);
    virtual void             write(GXmlElement& xml) const;
    virtual std::string      print(void) const;

    // Other methods
    void           load(const std::string& drename,
                        const std::string& drbname,
                        const std::string& drgname,
                        const std::string& drxname);
    void           response(const std::string& iaqname,
                            const std::string& caldb = "");
    void           obs_id(const double& id) { m_obs_id=id; }
    void           ontime(const double& ontime) { m_ontime=ontime; }
    void           livetime(const double& livetime) { m_livetime=livetime; }
    void           deadc(const double& deadc) { m_deadc=deadc; }
    double         obs_id(void) const { return m_obs_id; }
    const GSkymap& drb(void) const { return m_drb; }
    const GSkymap& drg(void) const { return m_drg; }
    const GSkymap& drx(void) const { return m_drx; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMObservation& obs);
    void free_members(void);
    void load_dre(const std::string& drename);
    void load_drb(const std::string& drbname);
    void load_drg(const std::string& drgname);
    void load_drx(const std::string& drxname);
    void read_attributes(const GFitsHDU* hdu);
    void write_attributes(GFitsHDU* hdu) const;

    // Protected members
    std::string   m_instrument;  //!< Instrument name
    std::string   m_drename;     //!< DRE filename
    std::string   m_drbname;     //!< DRB filename
    std::string   m_drgname;     //!< DRG filename
    std::string   m_drxname;     //!< DRX filename
    GSkymap       m_drb;         //!< Background model
    GSkymap       m_drg;         //!< Geometry factors
    GSkymap       m_drx;         //!< Exposure map
    GCOMPointing* m_pointing;    //!< Pointer to pointing direction
    GCOMResponse* m_response;    //!< Pointer to response functions
    double        m_obs_id;      //!< Observation ID
    double        m_ontime;      //!< Ontime
    double        m_livetime;    //!< Livetime
    double        m_deadc;       //!< Deadtime correction
};

#endif /* GCOMOBSERVATION_HPP */
