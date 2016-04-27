/***************************************************************************
 *          GCTAOnOffObservation.hpp - CTA On/Off observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Chia-Chun Lu & Christoph Deil               *
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
 * @file GCTAOnOffObservation.hpp
 * @brief CTA On/Off observation class definition
 * @author Chia-Chun Lu & Christoph Deil
 */

#ifndef GCTAONOFFOBSERVATION_HPP
#define GCTAONOFFOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GObservation.hpp"
#include "GPha.hpp"
#include "GArf.hpp"
#include "GRmf.hpp"
#include "GSkyRegions.hpp"

/* __ Forward declarations _______________________________________________ */
class GModels;
class GOptimizerPars;
class GCTAObservation;


/***********************************************************************//**
 * @class GCTAOnOffObservation
 *
 * @brief CTA On/Off observation class
 ***************************************************************************/
class GCTAOnOffObservation : public GObservation {

public:
    // Constructors and destructors
    GCTAOnOffObservation(void);
    GCTAOnOffObservation(const GEbounds& ereco, const GSkyRegions& on,
                         const GSkyRegions& off);
    GCTAOnOffObservation(const GCTAOnOffObservation& obs);
    virtual ~GCTAOnOffObservation(void);
 
    // Operators
    GCTAOnOffObservation& operator=(const GCTAOnOffObservation& obs);

    // Implemented pure virtual methods
    virtual void                  clear(void);
    virtual GCTAOnOffObservation* clone(void) const;
	virtual std::string           classname(void) const;
    virtual void                  response(const GResponse& rsp);
    virtual const GCTAResponse*   response(void) const;
    virtual std::string           instrument(void) const;
	virtual double                ontime(void) const;
	virtual double                livetime(void) const;
	virtual double                deadc(const GTime& time = GTime()) const;
    virtual void                  read(const GXmlElement& xml);
    virtual void                  write(GXmlElement& xml) const;
	virtual std::string           print(const GChatter& chatter = NORMAL) const;

    // Overloaded virtual methods
	virtual double likelihood(const GModels& models,
                              GVector*       gradient,
                              GMatrixSparse* curvature,
                              double*        npred) const;

    // Other methods
    void        instrument(const std::string& instrument);
    void        on_regions(const GSkyRegions& regions);
    void        off_regions(const GSkyRegions& regions);
    const GPha& on_spec(void) const;
    const GPha& off_spec(void) const;
    const GArf& arf(void) const;
	const GArf& bgd(void) const;
    const GRmf& rmf(void) const;
    void        fill(const GCTAObservation& obs);
    void        compute_response(const GCTAObservation& obs,
	                             const GModels& models,
                                 const GEbounds& etrue);

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GCTAOnOffObservation& obs);
    void   free_members(void);
	void   compute_alpha(const GCTAObservation& obs);
    void   compute_arf(const GCTAObservation& obs);
	void   compute_bgd(const GCTAObservation& obs, const GModels& models);
    void   compute_rmf(const GCTAObservation& obs, const GEbounds& etrue);
    double model_on(const GModels& models, int ibin, GVector* mod_grad) const;
	double model_off(const GModels& models, int ibin, GVector* mod_grad) const;

    // Protected data members
    std::string         m_instrument;   //!< Instrument name
    GCTAResponse*       m_response;     //!< Pointer to instrument response functions
	double              m_ontime;       //!< Ontime (seconds)
    double              m_livetime;     //!< Livetime (seconds)
    double              m_deadc;        //!< Deadtime correction (livetime/ontime)
    GPha 		        m_on_spec;      //!< On counts spectrum
    GPha 		        m_off_spec;     //!< Off counts spectrum
    GArf                m_arf;          //!< Effective area vector
	GArf                m_bgd;          //!< Background rate vector
    GRmf                m_rmf;          //!< Energy dispersion matrix
    GSkyRegions         m_on_regions;   //!< Container of On region
    GSkyRegions         m_off_regions;  //!< Container of Off regions
    std::vector<double> m_alpha;        //!< Ratios of On/Off exposure
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAOnOffObservation").
 ***************************************************************************/
inline
std::string GCTAOnOffObservation::classname(void) const
{
    return ("GCTAOnOffObservation");
}


/***********************************************************************//**
 * @brief Return instrument
 *
 * @return Instrument.
 ***************************************************************************/
inline
std::string GCTAOnOffObservation::instrument(void) const
{
    return (m_instrument);
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime in seconds.
 ***************************************************************************/
inline
double GCTAOnOffObservation::ontime(void) const
{
    return (m_ontime);
}


/***********************************************************************//**
 * @brief Return livetime
 *
 * @return Livetime in seconds.
 ***************************************************************************/
inline
double GCTAOnOffObservation::livetime(void) const
{
    return (m_livetime);
}


/***********************************************************************//**
 * @brief Return deadtime correction factor
 *
 * @param[in] time Time.
 * @return Deadtime correction factor.
 *
 * Returns the deadtime correction factor. Optionally, this method takes a
 * @p time argument that takes provision for returning the deadtime
 * correction factor as function of time.
 *
 * The deadtime correction factor is defined as the livetime divided by the
 * ontime.
 ***************************************************************************/
inline
double GCTAOnOffObservation::deadc(const GTime& time) const
{
    return (m_deadc);
}


/***********************************************************************//**
 * @brief Set instrument
 *
 * @param[in] instrument Instrument.
 ***************************************************************************/
inline
void GCTAOnOffObservation::instrument(const std::string& instrument)
{
    m_instrument = instrument;
    return;
}


/***********************************************************************//**
 * @brief Set On regions
 *
 * @param[in] regions On regions.
 ***************************************************************************/
inline
void GCTAOnOffObservation::on_regions(const GSkyRegions& regions)
{
    m_on_regions = regions;
    return;
}


/***********************************************************************//**
 * @brief Set Off regions
 *
 * @param[in] regions Off regions.
 ***************************************************************************/
inline
void GCTAOnOffObservation::off_regions(const GSkyRegions& regions)
{
    m_off_regions = regions;
    return;
}


/***********************************************************************//**
 * @brief Return ON spectrum
 *
 * @return ON spectrum.
 ***************************************************************************/
inline
const GPha& GCTAOnOffObservation::on_spec(void) const
{
    return (m_on_spec);
}


/***********************************************************************//**
 * @brief Return OFF spectrum
 *
 * @return OFF spectrum.
 ***************************************************************************/
inline
const GPha& GCTAOnOffObservation::off_spec(void) const
{
    return (m_off_spec);
}


/***********************************************************************//**
 * @brief Return Auxiliary Response File
 *
 * @return Auxiliary Response File.
 ***************************************************************************/
inline
const GArf& GCTAOnOffObservation::arf(void) const
{
    return (m_arf);
}


/***********************************************************************//**
 * @brief Return Background Rate
 *
 * @return Background Rate.
 ***************************************************************************/
inline
const GArf& GCTAOnOffObservation::bgd(void) const
{
    return (m_bgd);
}


/***********************************************************************//**
 * @brief Return Redistribution Matrix File
 *
 * @return Redistribution Matrix File.
 ***************************************************************************/
inline
const GRmf& GCTAOnOffObservation::rmf(void) const
{
    return (m_rmf);
}

#endif /* GCTAONOFFOBSERVATION_HPP */
