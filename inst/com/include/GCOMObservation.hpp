/***************************************************************************
 *            GCOMObservation.hpp - COMPTEL observation class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
#include <string>
#include "GObservation.hpp"
#include "GTime.hpp"
#include "GFilename.hpp"
#include "GSource.hpp"
#include "GSkyDir.hpp"
#include "GSkyMap.hpp"
#include "GCOMResponse.hpp"
#include "GCOMTim.hpp"
#include "GCOMOads.hpp"
#include "GCOMDri.hpp"
#include "GCOMEventCube.hpp"
#include "GCOMEventList.hpp"

/* __ Forward declarations _______________________________________________ */
class GCaldb;
class GResponse;
class GXmlElement;
class GFitsHDU;
class GCOMStatus;


/***********************************************************************//**
 * @class GCOMObservation
 *
 * @brief Interface class for COMPTEL observations
 *
 * This class implements a COMPTEL observation. Each COMPTEL observation is
 * defined for a given energy range, and is composed of a DRE, DRB, DRG and
 * DRX file. The DRE file contains the event data, the DRB file contains a
 * background model, the DRG file contains geometry factors, and the DRX file
 * contains the exposure.
 ***************************************************************************/
class GCOMObservation : public GObservation {

public:
    // Constructors and destructors
    GCOMObservation(void);
    GCOMObservation(const GFilename& drename,
                    const GFilename& drbname,
                    const GFilename& drgname,
                    const GFilename& drxname);
    GCOMObservation(const GFilename&              evpname,
                    const GFilename&              timname,
                    const std::vector<GFilename>& oadnames);
    GCOMObservation(const GCOMObservation& obs);
    virtual ~GCOMObservation(void);

    // Operators
    virtual GCOMObservation& operator=(const GCOMObservation& obs);

    // Implement pure virtual methods
    virtual void                clear(void);
    virtual GCOMObservation*    clone(void) const;
    virtual std::string         classname(void) const;
    virtual void                response(const GResponse& rsp);
    virtual const GCOMResponse* response(void) const;
    virtual std::string         instrument(void) const;
    virtual double              ontime(void) const;
    virtual double              livetime(void) const;
    virtual double              deadc(const GTime& time = GTime()) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    bool            is_unbinned(void) const;
    bool            is_binned(void) const;
    void            load(const GFilename& drename,
                         const GFilename& drbname,
                         const GFilename& drgname,
                         const GFilename& drxname);
    void            load(const GFilename&              evpname,
                         const GFilename&              timname,
                         const std::vector<GFilename>& oadnames);
    void            response(const GCaldb& caldb, const std::string& rspname);
    void            obs_id(const double& id);
    void            ontime(const double& ontime);
    void            livetime(const double& livetime);
    void            deadc(const double& deadc);
    void            ewidth(const double& ewidth);
    const double&   obs_id(void) const;
    const double&   ewidth(void) const;
    const GSkyMap&  drb(void) const;
    const GSkyMap&  drg(void) const;
    const GSkyMap&  drx(void) const;
    const GCOMDri&  drm(const GSource& source) const;
    const GCOMTim&  tim(void) const;
    const GCOMOads& oads(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMObservation& obs);
    void free_members(void);
    void load_dre(const GFilename& drename);
    void load_drb(const GFilename& drbname);
    void load_drg(const GFilename& drgname);
    void load_drx(const GFilename& drxname);
    void add_drm(const GSource& source);
    bool check_map(const GSkyMap& map) const;
    void read_attributes(const GFitsHDU* hdu);
    void write_attributes(GFitsHDU* hdu) const;

    // Protected members
    std::string            m_instrument; //!< Instrument name
    GCOMResponse           m_response;   //!< Response functions
    double                 m_obs_id;     //!< Observation ID
    double                 m_ontime;     //!< Ontime (sec)
    double                 m_livetime;   //!< Livetime (sec)
    double                 m_deadc;      //!< Deadtime correction

    // Protected members for binned observation
    GFilename              m_drename;    //!< DRE filename
    GFilename              m_drbname;    //!< DRB filename
    GFilename              m_drgname;    //!< DRG filename
    GFilename              m_drxname;    //!< DRX filename
    GSkyMap                m_drb;        //!< Background model
    GSkyMap                m_drg;        //!< Geometry factors
    GSkyMap                m_drx;        //!< Exposure map
    double                 m_ewidth;     //!< Energy width (MeV)

    // Protected members for unbinned observation
    GFilename              m_evpname;    //!< EVP filename
    GFilename              m_timname;    //!< TIM filename
    std::vector<GFilename> m_oadnames;   //!< OAD filenames
    GCOMTim                m_tim;        //!< COMPTEL Good Time Intervals
    GCOMOads               m_oads;       //!< Orbit Aspect Data

    // Protected members for response cache
    std::vector<std::vector<double> > m_pars; //!< Parameters for model cubes
    std::vector<GCOMDri>              m_drms; //!< Convolved model cubes
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMObservation").
 ***************************************************************************/
inline
std::string GCOMObservation::classname(void) const
{
    return ("GCOMObservation");
}


/***********************************************************************//**
 * @brief Return response function
 *
 * @return Response function.
 ***************************************************************************/
inline
const GCOMResponse* GCOMObservation::response(void) const
{
    // Return response pointer
    return &m_response;
}


/***********************************************************************//**
 * @brief Return instrument
 *
 * @return Instrument name.
 ***************************************************************************/
inline
std::string GCOMObservation::instrument(void) const
{
    // Return instrument
    return (m_instrument);
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime (seconds).
 ***************************************************************************/
inline
double GCOMObservation::ontime(void) const
{
    // Return ontime
    return (m_ontime);
}


/***********************************************************************//**
 * @brief Return livetime
 *
 * @return Livetime (seconds).
 ***************************************************************************/
inline
double GCOMObservation::livetime(void) const
{
    // Return livetime
    return (m_livetime);
}


/***********************************************************************//**
 * @brief Return deadtime correction factor
 *
 * @param[in] time Time.
 *
 * @return Deadtime correction factor.
 ***************************************************************************/
inline
double GCOMObservation::deadc(const GTime& time) const
{
    // Return livetime
    return (m_deadc);
}


/***********************************************************************//**
 * @brief Set observation ID
 *
 * @param[in] id Observation ID.
 ***************************************************************************/
inline
void GCOMObservation::obs_id(const double& id)
{
    m_obs_id = id;
    return;
}


/***********************************************************************//**
 * @brief Set ontime
 *
 * @param[in] ontime Ontime.
 ***************************************************************************/
inline
void GCOMObservation::ontime(const double& ontime)
{
    m_ontime = ontime;
    return;
}


/***********************************************************************//**
 * @brief Set livetime
 *
 * @param[in] livetime Livetime.
 ***************************************************************************/
inline
void GCOMObservation::livetime(const double& livetime)
{
    m_livetime = livetime;
    return;
}


/***********************************************************************//**
 * @brief Set deadtime correction factor
 *
 * @param[in] deadc Deadtime correction factor.
 ***************************************************************************/
inline
void GCOMObservation::deadc(const double& deadc)
{
    m_deadc = deadc;
    return;
}


/***********************************************************************//**
 * @brief Set energy width
 *
 * @param[in] ewidth Energy width (MeV).
 ***************************************************************************/
inline
void GCOMObservation::ewidth(const double& ewidth)
{
    m_ewidth = ewidth;
    return;
}


/***********************************************************************//**
 * @brief Return observation ID
 *
 * @return Observation ID.
 ***************************************************************************/
inline
const double& GCOMObservation::obs_id(void) const
{
    // Return observation ID
    return (m_obs_id);
}


/***********************************************************************//**
 * @brief Return energy width
 *
 * @return Energy width (MeV).
 ***************************************************************************/
inline
const double& GCOMObservation::ewidth(void) const
{
    // Return energy width
    return (m_ewidth);
}


/***********************************************************************//**
 * @brief Return background model
 *
 * @return Background model.
 ***************************************************************************/
inline
const GSkyMap& GCOMObservation::drb(void) const
{
    // Return background model
    return (m_drb);
}


/***********************************************************************//**
 * @brief Return geometry factors
 *
 * @return Geometry factors.
 ***************************************************************************/
inline
const GSkyMap& GCOMObservation::drg(void) const
{
    // Return geometry factors
    return (m_drg);
}


/***********************************************************************//**
 * @brief Return exposure
 *
 * @return Exposure.
 ***************************************************************************/
inline
const GSkyMap& GCOMObservation::drx(void) const
{
    // Return exposure
    return (m_drx);
}


/***********************************************************************//**
 * @brief Return COMPTEL Good Time Intervals
 *
 * @return COMPTEL Good Time Intervals.
 ***************************************************************************/
inline
const GCOMTim& GCOMObservation::tim(void) const
{
    // Return COMPTEL Good Time Intervals
    return (m_tim);
}


/***********************************************************************//**
 * @brief Return Orbit Aspect Data
 *
 * @return Orbit Aspect Data
 ***************************************************************************/
inline
const GCOMOads& GCOMObservation::oads(void) const
{
    // Return Orbit Aspect Data
    return (m_oads);
}


/***********************************************************************//**
 * @brief Check whether observation is unbinned
 *
 * @return True if observation is unbinned.
 ***************************************************************************/
inline
bool GCOMObservation::is_unbinned(void) const
{
    return (dynamic_cast<const GCOMEventList*>(m_events) != NULL);
}


/***********************************************************************//**
 * @brief Check whether observation is binned
 *
 * @return True if observation is unbinned.
 ***************************************************************************/
inline
bool GCOMObservation::is_binned(void) const
{
    return (dynamic_cast<const GCOMEventCube*>(m_events) != NULL);
}

#endif /* GCOMOBSERVATION_HPP */
