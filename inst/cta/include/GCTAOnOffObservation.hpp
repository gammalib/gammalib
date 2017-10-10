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
 * @author Chia-Chun Lu & Christoph Deil & Pierrick Martin
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
#include "GSkyRegionMap.hpp"

/* __ Forward declarations _______________________________________________ */
class GModels;
class GOptimizerPars;
class GCTAObservation;
class GBounds;


/***********************************************************************//**
 * @class GCTAOnOffObservation
 *
 * @brief CTA On/Off observation class
 *
 * This class defines a CTA On/Off observation. An On/Off observation is
 * defined by two spectra, one for an On region including the source of
 * interest, and one for an Off region including only background. The
 * response of an On/Off observation is given by the Auxiliary Response File
 * (ARF) and the Redistribution Matrix File (RMF).
 *
 * The class uses GPha objects to store the On and Off spectra, and GArf and
 * GRmf objects to store the response information. On and Off regions are 
 * defined via GSkyRegionMap objects which are essentially finely-pixellized 
 * skymaps that allow handling arbitrarily complex shapes
 ***************************************************************************/
class GCTAOnOffObservation : public GObservation {

public:
    // Constructors and destructors
    GCTAOnOffObservation(void);
    GCTAOnOffObservation(const GCTAObservation& obs,
                         const GEbounds&        etrue,
                         const GEbounds&        ereco,
                         const GSkyRegionMap&     on,
                         const GSkyRegionMap&     off,
                         const GSkyDir            src_dir);
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
    virtual int    nobserved(void) const;

    // Other methods
    void               instrument(const std::string& instrument);
    void               on_regions(const GSkyRegionMap& regionmap);
    void               off_regions(const GSkyRegionMap& regionmap);
    const GSkyRegionMap& on_regions(void) const;
    const GSkyRegionMap& off_regions(void) const;
    const GPha&        on_spec(void) const;
    const GPha&        off_spec(void) const;
    const GArf&        arf(void) const;
    const GRmf&        rmf(void) const;
    void               src_dir(const GSkyDir& src_dir);
    const GSkyDir&     src_dir(void) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GCTAOnOffObservation& obs);
    void   free_members(void);
    void   set(const GCTAObservation& obs);
    void   compute_arf(const GCTAObservation& obs);
	void   compute_bgd(const GCTAObservation& obs);
	void   compute_alpha(const GCTAObservation& obs);
    void   compute_rmf(const GCTAObservation& obs);
    double N_gamma(const GModels& models, const int& ibin, GVector* grad) const;
	double N_bgd(const GModels& models, const int& ibin, GVector* grad) const;

    // Protected data members
    std::string         m_instrument;  //!< Instrument name
    GCTAResponse*       m_response;    //!< Pointer to IRFs
	double              m_ontime;      //!< Ontime (seconds)
    double              m_livetime;    //!< Livetime (seconds)
    double              m_deadc;       //!< Deadtime correction (livetime/ontime)
    GPha 		        m_on_spec;     //!< On counts spectrum
    GPha 		        m_off_spec;    //!< Off counts spectrum
    GArf                m_arf;         //!< Auxiliary Response Function vector
    GRmf                m_rmf;         //!< Redistribution matrix
    GSkyRegionMap       m_on_regions;  //!< Container of On region
    GSkyRegionMap       m_off_regions; //!< Container of Off regions
    GSkyDir             m_src_dir;     //!< Direction to source in On region
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
 * @brief Set On regions map
 *
 * @param[in] regions On regions map.
 ***************************************************************************/
inline
void GCTAOnOffObservation::on_regions(const GSkyRegionMap& regionmap)
{
    m_on_regions = regionmap;
    return;
}

/***********************************************************************//**
 * @brief Set Off regions map
 *
 * @param[in] regions Off regions map.
 ***************************************************************************/
inline
void GCTAOnOffObservation::off_regions(const GSkyRegionMap& regionmap)
{
    m_off_regions = regionmap;
    return;
}

/***********************************************************************//**
 * @brief Return On regions map
 *
 * @return On regions map.
 ***************************************************************************/
inline
const GSkyRegionMap& GCTAOnOffObservation::on_regions(void) const
{
    return (m_on_regions);
}

/***********************************************************************//**
 * @brief Return Off regions map
 *
 * @return Off regions map.
 ***************************************************************************/
inline
const GSkyRegionMap& GCTAOnOffObservation::off_regions(void) const
{
    return (m_off_regions);
}

/***********************************************************************//**
 * @brief Return On spectrum
 *
 * @return On spectrum.
 ***************************************************************************/
inline
const GPha& GCTAOnOffObservation::on_spec(void) const
{
    return (m_on_spec);
}

/***********************************************************************//**
 * @brief Return Off spectrum
 *
 * @return Off spectrum.
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
 * @brief Return Redistribution Matrix File
 *
 * @return Redistribution Matrix File.
 ***************************************************************************/
inline
const GRmf& GCTAOnOffObservation::rmf(void) const
{
    return (m_rmf);
}

/***********************************************************************//**
 * @brief Return number of observed events
 *
 * @return Number of observed events.
 ***************************************************************************/
inline
int GCTAOnOffObservation::nobserved(void) const
{
    return (m_on_spec.counts());
}

/***********************************************************************//**
 * @brief Set direction to source
 *
 * @param[in] direction to source.
 ***************************************************************************/
inline
void GCTAOnOffObservation::src_dir(const GSkyDir& src_dir)
{
    m_src_dir = src_dir;
    return;
}

/***********************************************************************//**
 * @brief Return direction to source
 *
 * @return direction to source
 ***************************************************************************/
inline
const GSkyDir& GCTAOnOffObservation::src_dir(void) const
{
    return (m_src_dir);
}

#endif /* GCTAONOFFOBSERVATION_HPP */