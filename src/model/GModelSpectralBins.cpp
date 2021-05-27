/***************************************************************************
 *           GModelSpectralBins.cpp - Spectral bins model class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GModelSpectralBins.cpp
 * @brief Spectral bins model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GRan.hpp"
#include "GEbounds.hpp"
#include "GModelSpectralBins.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralBins     g_spectral_bins_seed;
const GModelSpectralRegistry g_spectral_bins_registry(&g_spectral_bins_seed);

/* __ Method name definitions ____________________________________________ */
#define G_FLUX                 "GModelSpectralBins::flux(GEnergy&, GEnergy&)"
#define G_EFLUX               "GModelSpectralBins::eflux(GEnergy&, GEnergy&)"
#define G_MC      "GModelSpectralBins::mc(GEnergy&, GEnergy&, GTime&, GRan&)"
#define G_READ                       "GModelSpectralBins::read(GXmlElement&)"
#define G_WRITE                     "GModelSpectralBins::write(GXmlElement&)"
#define G_APPEND    "GModelSpectralBins::append(GEnergy&, GEnergy&, double&)"
#define G_INSERT      "GModelSpectralBins::insert(int&, GEnergy&, GEnergy&, "\
                                                                   "double&)"
#define G_REMOVE                           "GModelSpectralBins::remove(int&)"
#define G_EMIN_GET                           "GModelSpectralBins::emin(int&)"
#define G_EMAX_GET                           "GModelSpectralBins::emax(int&)"
#define G_EMIN_SET                 "GModelSpectralBins::emin(int&, GEnergy&)"
#define G_EMAX_SET                 "GModelSpectralBins::emax(int&, GEnergy&)"
#define G_INTENSITY_GET                 "GModelSpectralBins::intensity(int&)"
#define G_INTENSITY_SET        "GModelSpectralBins::intensity(int&, double&)"
#define G_ERROR_GET                         "GModelSpectralBins::error(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty spectral bin model.
 ***************************************************************************/
GModelSpectralBins::GModelSpectralBins(void) : GModelSpectral()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Spectral model constructor
 *
 * @param[in] model Spectral model.
 * @param[in] ebounds Energy boundaries.
 * @param[in] index Spectral index.
 *
 * Constructs a spectral bins model from any spectral model using the
 * energy boundaries that are specified by @p ebounds. Within each spectral
 * bin the intensity follows a power law with spectral index @p index.
 ***************************************************************************/
GModelSpectralBins::GModelSpectralBins(const GModelSpectral& model,
                                       const GEbounds&       ebounds,
                                       const double&         index) :
                    GModelSpectral()
{
    // Initialise members
    init_members();

    // Append bins for all energies
    for (int i = 0; i < ebounds.size(); ++i) {
        append(ebounds.emin(i), ebounds.emax(i), model.eval(ebounds.elogmean(i)));
    }

    // Store index
    this->index(index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Construct spectral bins model by extracting information from an XML
 * element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpectralBins::GModelSpectralBins(const GXmlElement& xml) :
                    GModelSpectral()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Spectral bins model.
 ***************************************************************************/
GModelSpectralBins::GModelSpectralBins(const GModelSpectralBins& model) :
                    GModelSpectral(model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpectralBins::~GModelSpectralBins(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Spectral bins model.
 * @return Spectral bins model.
 ***************************************************************************/
GModelSpectralBins& GModelSpectralBins::operator=(const GModelSpectralBins& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpectral::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear spectral bins model
 ***************************************************************************/
void GModelSpectralBins::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpectral::free_members();

    // Initialise members
    this->GModelSpectral::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone spectral bins model
***************************************************************************/
GModelSpectralBins* GModelSpectralBins::clone(void) const
{
    // Clone spectral bins model
    return new GModelSpectralBins(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True photon energy.
 * @param[in] srcTime True photon arrival time.
 * @param[in] gradients Compute gradients?
 * @return Model value (ph/cm2/s/MeV).
 *
 * @todo Document method.
 ***************************************************************************/
double GModelSpectralBins::eval(const GEnergy& srcEng,
                                const GTime&   srcTime,
                                const bool&    gradients) const
{
    // Initialise value
    double value = 0.0;

    // Optionally initialise partial derivatives
    if (gradients) {
        m_index.factor_gradient(0.0);
        for (int i = 0; i < m_values.size(); ++i) {
            m_values[i].factor_gradient(0.0);
        }
    }

    // Get bin index
    int ibin = bin_index(srcEng);

    // Continue only if index is valid
    if (ibin >= 0) {

        // Get power-law index
        double index = m_index.value();

        // Compute and store value
        double eng    = srcEng.MeV();
        double e_norm = eng / m_epivot[ibin];
        double power  = (index != 0.0) ? std::pow(e_norm, index) : 1.0;

        // Compute function value
        value = m_values[ibin].value() * power;

        // Compute gradients
        if (gradients) {

            // Compute partial derivatives of the parameter values
            double g_norm  = (m_values[ibin].is_free())
                             ? m_values[ibin].scale() * power : 0.0;
            double g_index = (m_index.is_free())
                             ? value * m_index.scale() * std::log(e_norm) : 0.0;

            // Set gradients
            m_values[ibin].factor_gradient(g_norm);
            m_index.factor_gradient(g_index);

        } // endif: gradient computation was requested

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: GModelSpectralBins::eval";
            std::cout << "(srcEng=" << srcEng;
            std::cout << ", srcTime=" << srcTime << "):";
            std::cout << " NaN/Inf encountered";
            std::cout << " (value=" << value;
            std::cout << ", power=" << power;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: bin index was valid

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (units: ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Photon flux (ph/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 ***************************************************************************/
double GModelSpectralBins::flux(const GEnergy& emin,
                                const GEnergy& emax) const
{
    // Initialise flux
    double flux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();

        // Loop over all bins
        for (int ibin = 0; ibin < bins(); ++ibin) {

            // Get energy boundaires of bin in MeV
            double e_bin_min = m_emin[ibin].value();
            double e_bin_max = m_emax[ibin].value();

            // Skip bins that do not overlap with energy range
            if (e_min >= e_bin_max || e_max < e_bin_min) {
                continue;
            }

            // Raise lower boundary to minimum energy
            if (e_bin_min < e_min) {
                e_bin_min = e_min;
            }

            // Lower upper boundary to maximum energy
            if (e_bin_max > e_max) {
                e_bin_max = e_max;
            }

            // If energy interval is positive then compute flux
            if (e_bin_max > e_bin_min) {

                // Add photon flux
                flux += m_values[ibin].value() *
                        gammalib::plaw_photon_flux(e_bin_min,
                                                   e_bin_max,
                                                   m_epivot[ibin],
                                                   m_index.value());

            } // endif: energy interval was positive

        } // endfor: looped over all spectral bins

    } // endif: integration range was valid

    // Return
    return flux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (units: erg/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Energy flux (erg/cm2/s).
 *
 * Computes
 *
 * \f[
 *    \int_{\tt emin}^{\tt emax} S_{\rm E}(E | t) E \, dE
 * \f]
 *
 * where
 * - [@p emin, @p emax] is an energy interval, and
 * - \f$S_{\rm E}(E | t)\f$ is the spectral model (ph/cm2/s/MeV).
 ***************************************************************************/
double GModelSpectralBins::eflux(const GEnergy& emin,
                                 const GEnergy& emax) const
{
    // Initialise flux
    double eflux = 0.0;

    // Compute only if integration range is valid
    if (emin < emax) {

        
        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();

        // Loop over all bins
        for (int ibin = 0; ibin < bins(); ++ibin) {

            // Get energy boundaires of bin in MeV
            double e_bin_min = m_emin[ibin].value();
            double e_bin_max = m_emax[ibin].value();

            // Skip bins that do not overlap with energy range
            if (e_min >= e_bin_max || e_max < e_bin_min) {
                continue;
            }

            // Raise lower boundary to minimum energy
            if (e_bin_min < e_min) {
                e_bin_min = e_min;
            }

            // Lower upper boundary to maximum energy
            if (e_bin_max > e_max) {
                e_bin_max = e_max;
            }

            // If energy interval is positive then compute flux
            if (e_bin_max > e_bin_min) {

                // Add energy flux
                eflux += m_values[ibin].value() *
                         gammalib::plaw_energy_flux(e_bin_min,
                                                    e_bin_max,
                                                    m_epivot[ibin],
                                                    m_index.value());

            } // endif: energy interval was positive

        } // endfor: looped over all spectral bins

    } // endif: integration range was valid

    // Return flux
    return eflux;
}


/***********************************************************************//**
 * @brief Returns MC energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] time True photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Energy.
 *
 * @exception GException::invalid_return_value
 *            No valid Monte Carlo cache
 *
 * Returns Monte Carlo energy by randomly drawing from bin function.
 ***************************************************************************/
GEnergy GModelSpectralBins::mc(const GEnergy& emin,
                               const GEnergy& emax,
                               const GTime&   time,
                               GRan&          ran) const
{
    // Check energy interval
    gammalib::check_energy_interval(G_MC, emin, emax);

    // Allocate energy
    GEnergy energy;

    // Update cache
    mc_update(emin, emax);

    // Determine in which bin we reside. If this fails then thrown an
    // exception
    int inx = 0;
    if (m_mc_cum.size() > 1) {
        double u = ran.uniform();
        for (inx = m_mc_cum.size()-1; inx > 0; --inx) {
            if (m_mc_cum[inx-1] <= u) {
                break;
            }
        }
    }
    else if (m_mc_cum.size() == 0) {
        std::string msg = "No valid bins found for energy interval ["+
                          emin.print()+","+emax.print()+"]. Either restrict "
                          "the energy range to the one covered by the "
                          "spectral bins or extend the spectral bins "
                          "in energy.";
        throw GException::invalid_return_value(G_MC, msg);
    }

    // Get random energy for specific bin
    if (m_mc_exp[inx] != 0.0) {
        double e_min = m_mc_min[inx];
        double e_max = m_mc_max[inx];
        double u     = ran.uniform();
        double eng   = (u > 0.0)
                        ? std::exp(std::log(u * (e_max - e_min) + e_min) / m_mc_exp[inx])
                        : 0.0;
        energy.MeV(eng);
    }
    else {
        double e_min = m_mc_min[inx];
        double e_max = m_mc_max[inx];
        double u     = ran.uniform();
        double eng   = std::exp(u * (e_max - e_min) + e_min);
        energy.MeV(eng);
    }

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Lower energy limit larger than upper energy limit.
 *
 * Reads the spectral information from an XML element. The format of the XML
 * elements is
 *
 *     <spectrum type="BinFunction">
 *       <parameter name="Index"    ../>
 *       <bin>
 *         <parameter name="LowerLimit" ../>
 *         <parameter name="UpperLimit" ../>
 *         <parameter name="Intensity"  ../>
 *       </bin>
 *       ...
 *       <bin>
 *         <parameter name="LowerLimit" ../>
 *         <parameter name="UpperLimit" ../>
 *         <parameter name="Intensity"  ../>
 *       </bin>
 *     </spectrum>
 *
 * @todo Check that bins are ordered
 * @todo Check that energy boundaries are not overlapping
 ***************************************************************************/
void GModelSpectralBins::read(const GXmlElement& xml)
{
    // Free space for bins
    m_emin.clear();
    m_emax.clear();
    m_values.clear();

    // Get index parameter pointer
    const GXmlElement* index = gammalib::xml_get_par(G_READ, xml, m_index.name());

    // Read index parameter
    m_index.read(*index);

    // Get number of bins from XML file
    int bins = xml.elements("bin");

    // Loop over all bins
    for (int i = 0; i < bins; ++i) {

        // Allocate bin parameters
        GModelPar emin;
        GModelPar emax;
        GModelPar intensity;

        // Get bins
        const GXmlElement* bin = xml.element("bin", i);

        // Get parameters
        const GXmlElement* lpar = gammalib::xml_get_par(G_READ, *bin, "LowerLimit");
        const GXmlElement* upar = gammalib::xml_get_par(G_READ, *bin, "UpperLimit");
        const GXmlElement* ipar = gammalib::xml_get_par(G_READ, *bin, "Intensity");

        // Read parameters
        emin.read(*lpar);
        emax.read(*upar);
        intensity.read(*ipar);

        // Throw an exception if either LowerLimit is larger than UpperLimit
        if (emin.value() > emax.value()) {
            std::string msg = "Lower energy limit "+gammalib::str(emin.value())+
                              " MeV is larger than upper energy limit "+
                              gammalib::str(emax.value())+" MeV. Please "
                              "specify lower energy limits that are smaller "
                              "than the upper energy limits.";
            throw GException::invalid_value(G_READ, msg);
        }

        // Push bin parameters on list
        m_emin.push_back(emin);
        m_emax.push_back(emax);
        m_values.push_back(intensity);

    } // endfor: looped over bins

    // Update parameter mapping
    update_pars();

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::invalid_value
 *            Existing XML element is not of required type
 *            Invalid number of model parameters or bins found in XML element.
 *
 * Writes the spectral information into an XML element. The format of the XML
 * element is
 *
 *     <spectrum type="BinFunction">
 *       <parameter name="Index"    ../>
 *       <bin>
 *         <parameter name="LowerLimit" ../>
 *         <parameter name="UpperLimit" ../>
 *         <parameter name="Intensity"  ../>
 *       </bin>
 *       ...
 *       <bin>
 *         <parameter name="LowerLimit" ../>
 *         <parameter name="UpperLimit" ../>
 *         <parameter name="Intensity"  ../>
 *       </bin>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralBins::write(GXmlElement& xml) const
{
    // Verify model type
    gammalib::xml_check_type(G_WRITE, xml, type());

    // Get index parameter
    GXmlElement* index = gammalib::xml_need_par(G_WRITE, xml, m_index.name());

    // Write index parameter
    m_index.write(*index);

    // Determine number of bins
    int bins = this->bins();

    // If XML element has 0 bins then append bins
    if (xml.elements("bin") == 0) {
        for (int i = 0; i < bins; ++i) {
            xml.append(GXmlElement("bin"));
        }
    }

    // Verify that XML element has the required number of bins and elements.
    // Recall that the power-law index also counts as one element.
    if (xml.elements() != bins+1) {
        std::string msg = "Spectral bins model contains "+
                          gammalib::str(xml.elements())+
                          " elements but "+gammalib::str(bins+1)+" elements "
                          "were expected.";
        throw GException::invalid_value(G_WRITE, msg);
    }
    if (xml.elements("bin") != bins) {
        std::string msg = "Spectral bins model contains "+
                          gammalib::str(xml.elements("bin"))+
                          " \"bin\" elements but "+gammalib::str(bins)+
                          " \"bin\" elements were expected.";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Loop over all bins
    for (int i = 0; i < bins; ++i) {

        // Get bins
        GXmlElement* bin = xml.element("bin", i);

        // Get copy of parameters so that we can change their names
        GModelPar emin      = m_emin[i];
        GModelPar emax      = m_emax[i];
        GModelPar intensity = m_values[i];

        // Set names since we appended for internal usage the indices to the
        // parameter names, and we want to get rid of them for writing the
        // model into the XML element
        emin.name("LowerLimit");
        emax.name("UpperLimit");
        intensity.name("Intensity");

        // Get XML parameters
        GXmlElement* lpar = gammalib::xml_need_par(G_WRITE, *bin, emin.name());
        GXmlElement* upar = gammalib::xml_need_par(G_WRITE, *bin, emax.name());
        GXmlElement* ipar = gammalib::xml_need_par(G_WRITE, *bin, intensity.name());

        // Write parameters
        emin.write(*lpar);
        emax.write(*upar);
        intensity.write(*ipar);

    } // endfor: looped over bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append bin
 *
 * @param[in] emin Lower energy limit.
 * @param[in] emax Upper energy limit.
 * @param[in] intensity Bin intensity.
 *
 * @exception GException::invalid_argument
 *            Lower energy limit is larger than upper energy limit.
 *
 * Appends one bin to the bin function. By default the lower and upper energy
 * limit of the bin is fixed while the intensity of the bin is free.
 ***************************************************************************/
void GModelSpectralBins::append(const GEnergy& emin,
                                const GEnergy& emax,
                                const double&  intensity)
{
    // Throw an exception if lower limit is larger than upper limit
    if (emin > emax) {
        std::string msg = "Lower energy limit "+emin.print()+" is larger than "
                          "upper energy limit "+emax.print()+". Please "
                          "specify a lower energy limit that is smaller "
                          "than the upper energy limit.";
        throw GException::invalid_value(G_APPEND, msg);
    }

    // Allocate bin parameters
    GModelPar par_emin;
    GModelPar par_emax;
    GModelPar par_intensity;

    // Set emin attributes
    par_emin.name("LowerLimit");
    par_emin.value(emin.MeV());
    par_emin.unit("MeV");
    par_emin.has_grad(false);
    par_emin.fix();

    // Set emax attributes
    par_emax.name("UpperLimit");
    par_emax.value(emax.MeV());
    par_emax.unit("MeV");
    par_emax.has_grad(false);
    par_emax.fix();

    // Set intensity attributes
    par_intensity.name("Intensity");
    par_intensity.value(intensity);
    par_intensity.unit("ph/cm2/s/MeV");
    par_intensity.has_grad(true);
    par_intensity.free();

    // Append to bins
    m_emin.push_back(par_emin);
    m_emax.push_back(par_emax);
    m_values.push_back(par_intensity);

    // Update parameter mapping
    update_pars();

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert bin
 *
 * @param[in] index Bin index [0,...,bins()-1].
 * @param[in] emin Lower energy limit.
 * @param[in] emax Upper energy limit.
 * @param[in] intensity Bin intensity.
 *
 * @exception GException::out_of_range
 *            Bin index is out of range.
 * @exception GException::invalid_argument
 *            Lower energy limit is larger than upper energy limit.
 *
 * Inserts a bin into the bin function before the bin with the specified
 * @p index. By default the lower and upper energy limit of the bin is fixed
 * while the intensity of the bin is free.
 ***************************************************************************/
void GModelSpectralBins::insert(const int&     index,
                                const GEnergy& emin,
                                const GEnergy& emax,
                                const double&  intensity)
{
    // Raise exception if index is outside boundary
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_INSERT, "Bin index", index, bins());
    }
    #endif

    // Throw an exception if lower limit is larger than upper limit
    if (emin > emax) {
        std::string msg = "Lower energy limit "+emin.print()+" is larger than "
                          "upper energy limit "+emax.print()+". Please "
                          "specify a lower energy limit that is smaller "
                          "than the upper energy limit.";
        throw GException::invalid_value(G_INSERT, msg);
    }

    // Allocate bin parameters
    GModelPar par_emin;
    GModelPar par_emax;
    GModelPar par_intensity;

    // Set emin attributes
    par_emin.name("LowerLimit");
    par_emin.value(emin.MeV());
    par_emin.unit("MeV");
    par_emin.has_grad(false);
    par_emin.fix();

    // Set emax attributes
    par_emax.name("UpperLimit");
    par_emax.value(emax.MeV());
    par_emax.unit("MeV");
    par_emax.has_grad(false);
    par_emax.fix();

    // Set intensity attributes
    par_intensity.name("Intensity");
    par_intensity.value(intensity);
    par_intensity.unit("ph/cm2/s/MeV");
    par_intensity.has_grad(true);
    par_intensity.free();

    // Insert bin
    m_emin.insert(m_emin.begin()+index, par_emin);
    m_emax.insert(m_emax.begin()+index, par_emax);
    m_values.insert(m_values.begin()+index, par_intensity);

    // Update parameter mapping
    update_pars();

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove bin
 *
 * @param[in] index Bin index [0,...,bins()-1].
 *
 * @exception GException::out_of_range
 *            Bin index is out of range.
 *
 * Removes bin of specified @p index from the bin function.
 ***************************************************************************/
void GModelSpectralBins::remove(const int& index)
{
    // Raise exception if index is outside boundary
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_REMOVE, "Bin index", index, bins());
    }
    #endif

    // Erase energy and intensity
    m_emin.erase(m_emin.begin() + index);
    m_emax.erase(m_emax.begin() + index);
    m_values.erase(m_values.begin() + index);

    // Update parameter mapping
    update_pars();

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Reserve space for bins
 *
 * @param[in] num Number of reserved bins.
 *
 * Reserves space for @p num bins.
 ***************************************************************************/
void GModelSpectralBins::reserve(const int& num)
{
    // Reserve space
    m_emin.reserve(num);
    m_emax.reserve(num);
    m_values.reserve(num);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append bins from bin function
 *
 * @param[in] bins Bin function.
 *
 * Appends all bins from a bin function to current object.
 ***************************************************************************/
void GModelSpectralBins::extend(const GModelSpectralBins& bins)
{
    // Get number of bins in bin function. Note that we extract the size
    // first to avoid an endless loop that arises when a container is
    // appended to itself.
    int num = bins.bins();

    // Continue only if bin function is not empty
    if (num > 0) {

        // Reserve enough space
        reserve(this->bins() + num);

        // Loop over all bins and append them to the bin function
        for (int i = 0; i < num; ++i) {
            m_emin.push_back(bins.m_emin[i]);
            m_emax.push_back(bins.m_emax[i]);
            m_values.push_back(bins.m_values[i]);
        }

        // Update parameter mapping
        update_pars();

        // Set pre-computation cache
        set_cache();

    } // endif: bin function was not empty

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return spectral power law index
 *
 * @return Spectrel power law index.
 *
 * Returns the spectral power law index.
 ***************************************************************************/
double GModelSpectralBins::index(void) const
{
    // Return spectral power law index
    return (m_index.value());
}


/***********************************************************************//**
 * @brief Set spectral power law index
 *
 * @param[in] index Spectrel power law index.
 *
 * Set spectral power law index.
 ***************************************************************************/
void GModelSpectralBins::index(const double& index)
{
    // Set m_index
    m_index.value(index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return lower energy limit of bin
 *
 * @param[in] index Bin index [0,...,bins()-1].
 * @return Lower energy limit of bin @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Returns the lower energy limit of bin @p index.
 ***************************************************************************/
GEnergy GModelSpectralBins::emin(const int& index) const
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_EMIN_GET, "Bin index", index, bins());
    }
    #endif

    // Retrieve energy
    GEnergy energy;
    energy.MeV(m_emin[index].value());

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Return upper energy limit of bin
 *
 * @param[in] index Bin index [0,...,bins()-1].
 * @return Upper energy limit of bin @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Returns the upper energy limit of bin @p index.
 ***************************************************************************/
GEnergy GModelSpectralBins::emax(const int& index) const
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_EMAX_GET, "Bin index", index, bins());
    }
    #endif

    // Retrieve energy
    GEnergy energy;
    energy.MeV(m_emax[index].value());

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Set lower energy limit of bin
 *
 * @param[in] index Bin index [0,...,bins()-1].
 * @param[in] emin Lower energy limit of bin @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Sets the lower energy limit of bin @p index.
 ***************************************************************************/
void GModelSpectralBins::emin(const int& index, const GEnergy& emin)
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_EMIN_SET, "Bin index", index, bins());
    }
    #endif

    // Set energy
    m_emin[index].value(emin.MeV());

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set upper energy limit of bin
 *
 * @param[in] index Bin index [0,...,bins()-1].
 * @param[in] emax Upper energy limit of bin @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Sets the upper energy limit of bin @p index.
 ***************************************************************************/
void GModelSpectralBins::emax(const int& index, const GEnergy& emax)
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_EMAX_SET, "Bin index", index, bins());
    }
    #endif

    // Set energy
    m_emax[index].value(emax.MeV());

    // Set pre-computation cache
    set_cache();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return bin intensity
 *
 * @param[in] index Bin index [0,...,bins()-1].
 * @return Intensity of bin @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Returns the intensity of bin @p index.
 ***************************************************************************/
double GModelSpectralBins::intensity(const int& index) const
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_INTENSITY_GET, "Bin index", index, bins());
    }
    #endif

    // Return intensity
    return (m_values[index].value());
}


/***********************************************************************//**
 * @brief Return intensity error of bin
 *
 * @param[in] index Bin index [0,...,bins()-1].
 * @return Intensity error of bin @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Returns the intensity error of bin @p index.
 ***************************************************************************/
double GModelSpectralBins::error(const int& index) const
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_ERROR_GET, "Bin index", index, bins());
    }
    #endif

    // Return intensity error
    return (m_values[index].error());
}


/***********************************************************************//**
 * @brief Set bin intensity
 *
 * @param[in] index Bin index [0,...,bins()-1].
 * @param[in] intensity Intensity of bin @p index.
 *
 * @exception GException::out_of_range
 *            Index is out of range.
 *
 * Set the intensity of bin @p index.
 ***************************************************************************/
void GModelSpectralBins::intensity(const int& index, const double& intensity)
{
    // Raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= bins()) {
        throw GException::out_of_range(G_INTENSITY_SET, "Bin index", index, bins());
    }
    #endif

    // Set intensity
    m_values[index].value(intensity);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print bin function information
 *
 * @param[in] chatter Chattiness.
 * @return String containing bin function information.
 ***************************************************************************/
std::string GModelSpectralBins::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralBins ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of bins"));
        result.append(gammalib::str(bins()));
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralBins::init_members(void)
{
    // Initialise powerlaw index
    m_index.clear();
    m_index.name("Index");
    m_index.scale(1.0);
    m_index.value(-2.0);        // default: -2.0
    m_index.range(-10.0,+10.0); // range:   [-10,+10]
    m_index.fix();
    m_index.gradient(0.0);
    m_index.has_grad(true);

    // Initialise bin parameter vectors
    m_emin.clear();
    m_emax.clear();
    m_values.clear();

    // Initialise evaluation cache
    m_epivot.clear();

    // Initialise MC cache
    m_mc_emin.clear();
    m_mc_emax.clear();
    m_mc_cum.clear();
    m_mc_min.clear();
    m_mc_max.clear();
    m_mc_exp.clear();

    // Update parameter mapping
    update_pars();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Spectral bin function.
 ***************************************************************************/
void GModelSpectralBins::copy_members(const GModelSpectralBins& model)
{
    // Copy bin energies and values
    m_index    = model.m_index;
    m_emin     = model.m_emin;
    m_emax     = model.m_emax;
    m_values   = model.m_values;

    // Copy evluation cache
    m_epivot  = model.m_epivot;

    // Copy MC cache
    m_mc_emin = model.m_mc_emin;
    m_mc_emax = model.m_mc_emax;
    m_mc_cum  = model.m_mc_cum;
    m_mc_min  = model.m_mc_min;
    m_mc_max  = model.m_mc_max;
    m_mc_exp  = model.m_mc_exp;

    // Update parameter mapping
    update_pars();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralBins::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return bin index for energy
 *
 * @param[in] energy Energy.
 * @return Bin index (-1 if energy bin was not found).
 ***************************************************************************/
int GModelSpectralBins::bin_index(const GEnergy& energy) const
{
    // Initialise bin index
    int bin = -1;

    // Get energy in MeV
    double energy_MeV = energy.MeV();

    // Search bin index
    for (int i = 0; i < bins(); ++i) {
        if (energy_MeV >= m_emin[i].value() && energy_MeV < m_emax[i].value()) {
            bin = i;
            break;
        }
    }

    // Return bin index
    return bin;
}


/***********************************************************************//**
 * @brief Update parameter mapping
 *
 * Sets the parameter pointers in the m_pars array, enabling iterating over
 * all model parameters. This method needs to be called after changing the
 * number of bins in the spectral model. The method needs not to be called
 * after value update.
 ***************************************************************************/
void GModelSpectralBins::update_pars(void)
{
    // Clear parameter pointer(s)
    m_pars.clear();

    // Push back common spectral index
    m_pars.push_back(&m_index);

    // Set parameter pointers for all bins
    for (int i = 0; i < bins(); ++i) {

        // Set parameter names
        std::string lowerlimit_name = "LowerLimit" + gammalib::str(i);
        std::string upperlimit_name = "UpperLimit" + gammalib::str(i);
        std::string intensity_name  = "Intensity"  + gammalib::str(i);

        // Set lower limit attributes
        m_emin[i].name(lowerlimit_name);
        m_emin[i].unit("MeV");
        m_emin[i].has_grad(false);

        // Set upper limit attributes
        m_emax[i].name(upperlimit_name);
        m_emax[i].unit("MeV");
        m_emax[i].has_grad(false);

        // Set intensity attributes
        m_values[i].name(intensity_name);
        m_values[i].unit("ph/cm2/s/MeV");
        m_values[i].has_grad(true);

        // Set pointer
        m_pars.push_back(&(m_emin[i]));
        m_pars.push_back(&(m_emax[i]));
        m_pars.push_back(&(m_values[i]));

    } // endfor: looped over bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set cache
 *
 * Sets the computation cache. The method computes for the time being the
 * pivot energies for all bins in MeV.
 ***************************************************************************/
void GModelSpectralBins::set_cache(void) const
{
    // Clear cache
    m_epivot.clear();

    // Compute pivot energies for all bins
    for (int i = 0; i < bins(); ++i) {
        double epivot = std::sqrt(m_emin[i].value() * m_emax[i].value());
        m_epivot.push_back(epivot);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set MC pre-computation cache
 *
 * @param[in] emin Minimum energy.
 * @param[in] emax Maximum energy.
 *
 * This method sets up an array of indices and the cumulative distribution
 * function needed for MC simulations.
 ***************************************************************************/
void GModelSpectralBins::mc_update(const GEnergy& emin,
                                   const GEnergy& emax) const
{
    // Check if we need to update the cache
    if (emin != m_mc_emin || emax != m_mc_emax) {

        // Store new energy interval
        m_mc_emin = emin;
        m_mc_emax = emax;

        // Initialise cache
        m_mc_cum.clear();
        m_mc_min.clear();
        m_mc_max.clear();
        m_mc_exp.clear();

        // Get energy range in MeV
        double e_min = emin.MeV();
        double e_max = emax.MeV();

        // Continue only if e_max > e_min
        if (e_max > e_min) {

            // Loop over all bins
            for (int ibin = 0; ibin < bins(); ++ibin) {

                // Get energy boundaires of bin in MeV
                double e_bin_min = m_emin[ibin].value();
                double e_bin_max = m_emax[ibin].value();

                // Skip bins that do not overlap with energy range
                if (e_min >= e_bin_max || e_max < e_bin_min) {
                    continue;
                }

                // Raise lower boundary to minimum energy
                if (e_bin_min < e_min) {
                    e_bin_min = e_min;
                }

                // Lower upper boundary to maximum energy
                if (e_bin_max > e_max) {
                    e_bin_max = e_max;
                }

                // If energy interval is positive then compute flux
                if (e_bin_max > e_bin_min) {

                    // Get index
                    double index = m_index.value();

                    // Compute flux
                    double flux = m_values[ibin].value() *
                                  gammalib::plaw_photon_flux(e_bin_min,
                                                             e_bin_max,
                                                             m_epivot[ibin],
                                                             index);

                    // Push results in cache
                    m_mc_cum.push_back(flux);
                    m_mc_min.push_back(e_bin_min);
                    m_mc_max.push_back(e_bin_max);
                    m_mc_exp.push_back(index);

                } // endif: energy interval was positive

            } // endfor: looped over all spectral bins

            // Build cumulative distribution
            for (int i = 1; i < m_mc_cum.size(); ++i) {
                m_mc_cum[i] += m_mc_cum[i-1];
            }
            double norm = m_mc_cum[m_mc_cum.size()-1];
            if (norm > 0.0) {
                for (int i = 0; i < m_mc_cum.size(); ++i) {
                    m_mc_cum[i] /= norm;
                }
            }

            // Set MC values
            for (int i = 0; i < m_mc_cum.size(); ++i) {

                // Compute exponent
                double exponent = m_mc_exp[i] + 1.0;

                // Exponent dependend computation
                if (std::abs(exponent) > 1.0e-11) {

                    // If the exponent is too small then use lower energy
                    // boundary
                    if (exponent < -50.0) {
                        m_mc_exp[i] = 0.0;
                        m_mc_min[i] = std::log(m_mc_min[i]);
                        m_mc_max[i] = m_mc_min[i];
                    }

                    // ... otherwise if exponent is too large then use
                    // upper energy boundary
                    else if (exponent > +50.0) {
                        m_mc_exp[i] = 0.0;
                        m_mc_min[i] = std::log(m_mc_max[i]);
                        m_mc_max[i] = m_mc_min[i];
                    }

                    // ... otherwise use transformation formula
                    else {
                        m_mc_exp[i] = exponent;
                        m_mc_min[i] = std::pow(m_mc_min[i], exponent);
                        m_mc_max[i] = std::pow(m_mc_max[i], exponent);
                    }
                }
                else {
                    m_mc_exp[i] = 0.0;
                    m_mc_min[i] = std::log(m_mc_min[i]);
                    m_mc_max[i] = std::log(m_mc_max[i]);
                }

            } // endfor: set MC values

        } // endif: e_max > e_min

    } // endif: Update was required

    // Return
    return;
}
