/***************************************************************************
 *          GCTAOnOffObservation.hpp - CTA On/Off observation class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2018 by Chia-Chun Lu & Christoph Deil               *
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
#include "GPha.hpp"
#include "GArf.hpp"
#include "GRmf.hpp"
#include "GFunction.hpp"
#include "GObservation.hpp"
#include "GNodeArray.hpp"

/* __ Forward declarations _______________________________________________ */
class GModels;
class GModelSpatial;
class GOptimizerPars;
class GObservations;
class GSkyRegions;
class GCTAObservation;
class GCTAResponse;
class GCTAResponseIrf;


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
    GCTAOnOffObservation(const bool& dummy, const std::string& instrument);
    GCTAOnOffObservation(const GObservations& obs);
    GCTAOnOffObservation(const GPha& pha_on,
                         const GPha& pha_off,
                         const GArf& arf,
                         const GRmf& rmf);
    GCTAOnOffObservation(const GCTAObservation& obs,
                         const GModels&         models,
                         const std::string&     srcname,
                         const GEbounds&        etrue,
                         const GEbounds&        ereco,
                         const GSkyRegions&     on,
                         const GSkyRegions&     off,
                         const bool&            use_model_bkg = true);
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
    void        instrument(const std::string& instrument);
    bool        has_response(void) const;
    const GPha& on_spec(void) const;
    const GPha& off_spec(void) const;
    const GArf& arf(void) const;
    const GRmf& rmf(void) const;
    GPha        model_gamma(const GModels& models) const;
    GPha        model_background(const GModels& models) const;

protected:
    // Protected methods
    void   init_members(void);
    void   copy_members(const GCTAOnOffObservation& obs);
    void   free_members(void);
    void   set_exposure(void);
    void   check_consistency(const std::string& method) const;
    GArf   arf_stacked(const GArf& arf, const GEnergy& emin, const GEnergy& emax) const;
    GRmf   rmf_stacked(const GRmf& rmf, const GEnergy& emin, const GEnergy& emax) const;
    void   set(const GCTAObservation& obs,
               const GModels&         models,
               const std::string&     srcname,
               const GSkyRegions&     on,
               const GSkyRegions&     off,
               const bool&            use_model_bkg);
    void   compute_arf(const GCTAObservation& obs,
                       const GModelSpatial&   spatial,
                       const GSkyRegions&     on);
    void   compute_arf_cut(const GCTAObservation& obs,
                           const GModelSpatial&   spatial,
                           const GSkyRegions&     on);
    void   compute_bgd(const GCTAObservation& obs,
                       const GSkyRegions&     off,
                       const GModels&         models,
                       const bool&            use_model_bkg);
    void   compute_alpha(const GCTAObservation& obs,
                         const GSkyRegions&     on,
                         const GSkyRegions&     off,
                         const GModels&         models,
                         const bool&            use_model_bkg);
    void   compute_rmf(const GCTAObservation& obs,
                       const GSkyRegions&     on);
    void   apply_ebounds(const GCTAObservation& obs);
    double arf_rad_max(const GCTAObservation& obs, const GSkyRegions& on) const;
    double N_gamma(const GModels& models, const int& ibin, GVector* grad) const;
    double N_bgd(const GModels& models, const int& ibin, GVector* grad) const;

    // Likelihood methods
    virtual double likelihood_cstat(const GModels& models,
                                    GVector*       gradient,
                                    GMatrixSparse* curvature,
                                    double*        npred) const;
    virtual double likelihood_wstat(const GModels& models,
                                    GVector*       gradient,
                                    GMatrixSparse* curvature,
                                    double*        npred) const;
    virtual double wstat_value(const double& non,
                               const double& noff,
                               const double& alpha,
                               const double& ngam,
                               double&       nonpred,
                               double&       nbgd,
                               double&       dlogLdsky,
                               double&       d2logLdsky2) const;

    // Protected data members
    std::string   m_instrument; //!< Instrument name
    GCTAResponse* m_response;   //!< Pointer to IRFs
    double        m_ontime;     //!< Ontime (seconds)
    double        m_livetime;   //!< Livetime (seconds)
    double        m_deadc;      //!< Deadtime correction (livetime/ontime)
    GPha          m_on_spec;    //!< On counts spectrum
    GPha          m_off_spec;   //!< Off counts spectrum
    GArf          m_arf;        //!< Auxiliary Response Function vector
    GRmf          m_rmf;        //!< Redistribution matrix
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
 * @brief Signal if observation contains response information
 *
 * @return True if observation contains response information.
 ***************************************************************************/
inline
bool GCTAOnOffObservation::has_response(void) const
{
    return (m_response != NULL);
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
    return (int(m_on_spec.counts()+0.5));
}

#endif /* GCTAONOFFOBSERVATION_HPP */
