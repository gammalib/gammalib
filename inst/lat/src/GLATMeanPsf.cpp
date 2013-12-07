/***************************************************************************
 *                  GLATMeanPsf.cpp - Fermi/LAT mean PSF                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GLATMeanPsf.cpp
 * @brief Fermi/LAT mean PSF class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GLATMeanPsf.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GLATAeff.hpp"
#include "GLATPsf.hpp"
#include "GLATObservation.hpp"
#include "GLATEventCube.hpp"
#include "GLATException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SET                  "GLATMeanPsf::set(GSkyDir&, GLATObservation&)"
#define G_EXPOSURE                              "GLATMeanPsf::exposure(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_SIGNAL_NEGATIVE_MEAN_PSF 1     //!< Signal negative mean PSF value
#define G_DEBUG_INTEGRAL 0               //!< Debug mean PSF integration

/* __ Constants __________________________________________________________ */

/* __ Typedefs ___________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATMeanPsf::GLATMeanPsf(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief PSF constructor
 *
 * @param[in] dir Sky direction.
 * @param[in] obs LAT observation.
 ***************************************************************************/
GLATMeanPsf::GLATMeanPsf(const GSkyDir& dir, const GLATObservation& obs)
{
    // Initialise class members
    init_members();

    // Setup PSF
    set(dir, obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] psf Mean PSF.
 ***************************************************************************/
GLATMeanPsf::GLATMeanPsf(const GLATMeanPsf& psf)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(psf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATMeanPsf::~GLATMeanPsf(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] psf Mean PSF.
 ***************************************************************************/
GLATMeanPsf& GLATMeanPsf::operator= (const GLATMeanPsf& psf)
{
    // Execute only if object is not identical
    if (this != &psf) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(psf);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return mean PSF*exposure value
 *
 * @param[in] offset Angular distance from true source direction (degrees).
 * @param[in] logE log10 of energy in MeV.
 *
 * Computes
 * \f[{\rm PSF}(\delta, \log E) \times {\rm exposure}(\log E)\f]
 * by bilinear interpolation, where
 * \f$\delta\f$ is the offset angle from the source direction in degrees, and
 * \f$\log E\f$ is the logarithm of base 10 of the energy in MeV.
 * A zero value is returned if the offset angle is equal or larger than
 * 70 degrees or if \f$\log E\f$ is not positive.
 ***************************************************************************/
double GLATMeanPsf::operator() (const double& offset, const double& logE)
{
    // Initialise response
    double value = 0.0;

    // Continue only if arguments are within valid range
    if (offset < 70.0 && logE > 0.0) {

        // Flag no change of values
        //bool change = false;

        // Set offset interpolation
        //if (offset != m_last_offset) {
            m_offset.set_value(offset);
        //    m_last_offset = offset;
        //    change        = true;
        //}

        // Set energy interpolation
        //if (logE != m_last_energy) {
            m_energy.set_value(logE);
        //    m_last_energy = logE;
        //    change        = true;
        //}

        // If change occured then update interpolation indices and weighting
        // factors
        //if (change) {

            // Set energy indices for exposure computation
            m_inx1_exp = m_energy.inx_left();
            m_inx2_exp = m_energy.inx_right();

            // Set energy indices for PSF computation
            int inx_energy_left  = m_inx1_exp  * noffsets();
            int inx_energy_right = m_inx2_exp * noffsets();
            
            // Set array indices for bi-linear interpolation
            m_inx1 = m_offset.inx_left()  + inx_energy_left;
            m_inx2 = m_offset.inx_left()  + inx_energy_right;
            m_inx3 = m_offset.inx_right() + inx_energy_left;
            m_inx4 = m_offset.inx_right() + inx_energy_right;

            // Set weighting factors for bi-linear interpolation
            m_wgt1 = m_offset.wgt_left()  * m_energy.wgt_left();
            m_wgt2 = m_offset.wgt_left()  * m_energy.wgt_right();
            m_wgt3 = m_offset.wgt_right() * m_energy.wgt_left();
            m_wgt4 = m_offset.wgt_right() * m_energy.wgt_right();

        //} // endif: logE or ctheta changed

        // Compute energy dependent exposure and map corrections
        double fac_left  = m_exposure[m_inx1_exp] * m_mapcorr[m_inx1_exp];
        double fac_right = m_exposure[m_inx2_exp] * m_mapcorr[m_inx2_exp];

        // Perform bi-linear interpolation
        value = m_wgt1 * m_psf[m_inx1] * fac_left  +
                m_wgt2 * m_psf[m_inx2] * fac_right +
                m_wgt3 * m_psf[m_inx3] * fac_left  +
                m_wgt4 * m_psf[m_inx4] * fac_right;

        // Optionally check for negative values
        #if G_SIGNAL_NEGATIVE_MEAN_PSF
        if (value < 0.0) {
            std::cout << "GLATMeanPsf::operator()(offset=" << offset;
            std::cout << ", logE=" << logE << ") = " << value << std::endl;
        }
        #endif

    } // endif: arguments were valid

    // Return value
    return value;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATMeanPsf::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATMeanPsf* GLATMeanPsf::clone(void) const
{
    return new GLATMeanPsf(*this);
}


/***********************************************************************//**
 * @brief Return size of mean PSF
***************************************************************************/
int GLATMeanPsf::size(void) const
{
    // Compute size
    int size = m_energy.size()*m_offset.size();

    // Return size
    return size;
}


/***********************************************************************//**
 * @brief Compute mean PSF and exposure
 *
 * @param[in] dir Source location.
 * @param[in] obs LAT observation.
 *
 * @exception GLATException::no_ltcube
 *            Livetime cube has not been defined.
 *
 * Computes the mean PSF and the energy dependent exposure for a source at
 * a given sky location. The PSF is computed for all bin boundaries of the
 * observation, hence if there are N energy bins there will be N+1 energies
 * at which the mean PSF is computed.
 ***************************************************************************/
void GLATMeanPsf::set(const GSkyDir& dir, const GLATObservation& obs)
{
    // Clear PSF, exposure and energy arrays
    m_psf.clear();
    m_exposure.clear();
    m_energy.clear();

    // Get pointers on response, livetime cube and energy boundaries
    const GLATResponse& rsp = obs.response();

    // Get pointer on livetime cube
    GLATLtCube* ltcube = obs.ltcube();
    if (ltcube == NULL)
        throw GLATException::no_ltcube(G_SET);
    
    // Get energy boundaries
    GEbounds ebds = obs.events()->ebounds();

    // Store source direction
    m_dir = dir;

    // Limit computation to zenith angles < m_theta_max (typically 70
    // degrees - this is the hardwired value in the ST). For this purpose
    // set the costhetamin parameter of Aeff to m_theta_max. Store also the
    // original values in a vector for later restoration.
    std::vector<double> save_costhetamin;
    for (int i = 0; i < rsp.size(); ++i) {
        save_costhetamin.push_back(rsp.aeff(i)->costhetamin());
        rsp.aeff(i)->costhetamin(cos(m_theta_max*gammalib::deg2rad));
    }
    
    // Allocate room for arrays
    m_psf.reserve(size());
    m_exposure.reserve(m_energy.size());

    // Set energy nodes from the bin boundaries of the observations energy
    // boundaries. Store the energy nodes locally as GEnergy objects and
    // save them also in the class as log10 of energy in MeV.
    std::vector<GEnergy> energy;
    energy.reserve(ebds.size()+1);
    energy.push_back(ebds.emin(0));
    m_energy.append(ebds.emin(0).log10MeV());
    for (int i = 0; i < ebds.size(); ++i) {
        m_energy.append(ebds.emax(i).log10MeV());
        energy.push_back(ebds.emax(i));
    }

    // Loop over energies
    for (int ieng = 0; ieng < energy.size(); ++ieng) {

        // Compute exposure by looping over the responses
        double exposure = 0.0;
        for (int i = 0; i < rsp.size(); ++i) {
            exposure += (*ltcube)(dir, energy[ieng], *rsp.aeff(i));
        }

        // Set exposure
        m_exposure.push_back(exposure);

        // Loop over all offset angles
        for (int ioffset = 0; ioffset < m_offset.size(); ++ioffset) {

            // Compute point spread function by looping over the responses
            double psf = 0.0;
            for (int i = 0; i < rsp.size(); ++i)
                psf += (*ltcube)(dir, energy[ieng], m_offset[ioffset],
                                 *rsp.psf(i), *rsp.aeff(i));

            // Normalize PSF by exposure and clip when exposure drops to 0
            psf = (exposure > 0.0) ? psf/exposure : 0.0;

            // Set PSF value
            m_psf.push_back(psf);

        } // endfor: looped over offsets
    } // endfor: looped over energies

    // Restore initial Aeff zenith angle restriction
    for (int i = 0; i < rsp.size(); ++i) {
        rsp.aeff(i)->costhetamin(save_costhetamin[i]);
    }

    // Compute map corrections
    set_map_corrections(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return mean PSF value
 *
 * @param[in] offset Angular distance from true source direction (degrees).
 * @param[in] logE log10 of energy in MeV.
 *
 * Computes
 * \f[{\rm PSF}(\delta, \log E)\f]
 * by bilinear interpolation, where
 * \f$\delta\f$ is the offset angle from the source direction in degrees, and
 * \f$\log E\f$ is the logarithm of base 10 of the energy in MeV.
 * A zero value is returned if the offset angle is equal or larger than
 * 70 degrees or if \f$\log E\f$ is not positive.
 ***************************************************************************/
double GLATMeanPsf::psf(const double& offset, const double& logE)
{
    // Initialise response
    double value = 0.0;

    // Continue only if arguments are within valid range
    if (offset < 70.0 && logE > 0.0) {

        // Flag no change of values
        //bool change = false;

        // Set offset interpolation
        //if (offset != m_last_offset) {
            m_offset.set_value(offset);
        //    m_last_offset = offset;
        //    change        = true;
        //}

        // Set energy interpolation
        //if (logE != m_last_energy) {
            m_energy.set_value(logE);
        //    m_last_energy = logE;
        //    change        = true;
        //}

        // If change occured then update interpolation indices and weighting
        // factors
        //if (change) {

            // Set energy indices for exposure computation
            m_inx1_exp = m_energy.inx_left();
            m_inx2_exp = m_energy.inx_right();

            // Set energy indices
            int inx_energy_left  = m_energy.inx_left()  * noffsets();
            int inx_energy_right = m_energy.inx_right() * noffsets();
            
            // Set array indices for bi-linear interpolation
            m_inx1 = m_offset.inx_left()  + inx_energy_left;
            m_inx2 = m_offset.inx_left()  + inx_energy_right;
            m_inx3 = m_offset.inx_right() + inx_energy_left;
            m_inx4 = m_offset.inx_right() + inx_energy_right;

            // Set weighting factors for bi-linear interpolation
            m_wgt1 = m_offset.wgt_left()  * m_energy.wgt_left();
            m_wgt2 = m_offset.wgt_left()  * m_energy.wgt_right();
            m_wgt3 = m_offset.wgt_right() * m_energy.wgt_left();
            m_wgt4 = m_offset.wgt_right() * m_energy.wgt_right();

        //} // endif: logE or ctheta changed

        // Compute energy dependentmap corrections
        double fac_left  = m_mapcorr[m_inx1_exp];
        double fac_right = m_mapcorr[m_inx2_exp];

        // Perform bi-linear interpolation
        value = m_wgt1 * m_psf[m_inx1] * fac_left  +
                m_wgt2 * m_psf[m_inx2] * fac_right +
                m_wgt3 * m_psf[m_inx3] * fac_left  +
                m_wgt4 * m_psf[m_inx4] * fac_right;

        // Optionally check for negative values
        #if G_SIGNAL_NEGATIVE_MEAN_PSF
        if (value < 0.0) {
            std::cout << "GLATMeanPsf::psf(offset=" << offset;
            std::cout << ", logE=" << logE << ") = " << value << std::endl;
        }
        #endif

    } // endif: arguments were valid

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return exposure value
 *
 * @param[in] logE log10 of energy in MeV.
 *
 * Computes
 * \f[{\rm exposure}(\log E)\f]
 * by linear interpolation, where
 * \f$\log E\f$ is the logarithm of base 10 of the energy in MeV.
 * A zero value is returned if \f$\log E\f$ is not positive.
 ***************************************************************************/
double GLATMeanPsf::exposure(const double& logE)
{
    // Initialise response
    double value = 0.0;

    // Continue only if arguments are within valid range
    if (logE > 0.0) {

        // Set energy interpolation
        //if (logE != m_last_energy) {
            m_energy.set_value(logE);
        //    m_last_energy = logE;
            m_inx1_exp    = m_energy.inx_left();
            m_inx2_exp    = m_energy.inx_right();
        //}

        // Perform linear interpolation
        value = m_energy.wgt_left()  * m_exposure[m_inx1_exp] +
                m_energy.wgt_right() * m_exposure[m_inx2_exp];

    } // endif: arguments were in valid range

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Print lifetime cube information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing ROI information.
 ***************************************************************************/
std::string GLATMeanPsf::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute exposure range
        double min_exposure = 0.0;
        double max_exposure = 0.0;
        for (int i = 0; i < nenergies(); ++i) {
            if (i == 0) {
                min_exposure = m_exposure[i];
                max_exposure = m_exposure[i];
            }
            else {
                if (m_exposure[i] < min_exposure) {
                    min_exposure = m_exposure[i];
                }
                if (m_exposure[i] > max_exposure) {
                    max_exposure = m_exposure[i];
                }
            }
        }
    
        // Append header
        result.append("=== GLATMeanPsf ===");

        // Append information
        result.append("\n"+gammalib::parformat("Source name")+name());
        result.append("\n"+gammalib::parformat("Source direction"));
        result.append(gammalib::str(m_dir.ra_deg()));
        result.append(", ");
        result.append(gammalib::str(m_dir.dec_deg()));
        result.append("\n"+gammalib::parformat("Offset angles")+gammalib::str(noffsets()));
        result.append("\n"+gammalib::parformat("Energy values")+gammalib::str(nenergies()));
        result.append("\n"+gammalib::parformat("Exposure range"));
        result.append(gammalib::str(min_exposure)+" - "+gammalib::str(max_exposure)+" s cm2");
        for (int i = 0; i < nenergies(); ++i) {
            GEnergy energy;
            energy.log10MeV(m_energy[i]);
            result.append("\n"+gammalib::parformat(energy.print()));
            result.append(gammalib::str(m_exposure[i])+" s cm2");
            if (m_mapcorr.size() == nenergies()) {
                result.append("   (map correction="+gammalib::str(m_mapcorr[i])+")");
            }
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
void GLATMeanPsf::init_members(void)
{
    // Initialise members
    m_name.clear();
    m_dir.clear();
    m_psf.clear();
    m_exposure.clear();
    m_mapcorr.clear();
    m_energy.clear();
    m_offset.clear();
    m_theta_max   = 70.0;  //!< Maximum zenith angle
    m_last_energy = -1.0;
    m_last_offset = -1.0;
    m_inx1_exp    = 0;
    m_inx2_exp    = 0;
    m_inx1        = 0;
    m_inx2        = 0;
    m_inx3        = 0;
    m_inx4        = 0;
    m_wgt1        = 0.0;
    m_wgt2        = 0.0;
    m_wgt3        = 0.0;
    m_wgt4        = 0.0;

    // Set offset array
    set_offsets();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Mean PSF.
 ***************************************************************************/
void GLATMeanPsf::copy_members(const GLATMeanPsf& psf)
{
    // Copy members
    m_name        = psf.m_name;
    m_dir         = psf.m_dir;
    m_psf         = psf.m_psf;
    m_exposure    = psf.m_exposure;
    m_mapcorr     = psf.m_mapcorr;
    m_energy      = psf.m_energy;
    m_offset      = psf.m_offset;
    m_theta_max   = psf.m_theta_max;
    m_last_energy = psf.m_last_energy;
    m_last_offset = psf.m_last_offset;
    m_inx1_exp    = psf.m_inx1_exp;
    m_inx2_exp    = psf.m_inx2_exp;
    m_inx1        = psf.m_inx1;
    m_inx2        = psf.m_inx2;
    m_inx3        = psf.m_inx3;
    m_inx4        = psf.m_inx4;
    m_wgt1        = psf.m_wgt1;
    m_wgt2        = psf.m_wgt2;
    m_wgt3        = psf.m_wgt3;
    m_wgt4        = psf.m_wgt4;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATMeanPsf::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set array of offset values in degrees
 *
 * The array of offset values defines the points at with the mean PSF will
 * be evaluated and stored. The array is logarithmically spaced.
 ***************************************************************************/
void GLATMeanPsf::set_offsets(void)
{
    // Set array parameters
    const double offset_min = 1.0e-4;
    const double offset_max = 70.0;
    const int    offset_num = 200;

    // Clear offset array
    m_offset.clear();

    // Set array
    double step = log(offset_max/offset_min)/(offset_num - 1.0);
    for (int i = 0; i < offset_num; ++i)
        m_offset.append(offset_min*exp(i*step));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute map corrections
 *
 * @param[in] obs LAT observation.
 *
 * Sets up a vector of energy dependent corrections that assures that the
 * integral over the PSF in a specific event cube is correctly normalised.
 * For this purpose the PSF is integrated over a radius that fully lies
 * within the event cube and devided by the PSF pixel sum within this radius.
 * The map corrections are mainly important at high energies where the
 * PSF size can become comparable to the pixel size.
 *
 * @todo We can also implement a method for event atoms, yet it is not clear
 *       whether we really need this.
 ***************************************************************************/
void GLATMeanPsf::set_map_corrections(const GLATObservation& obs)
{
    // Initialise map corrections to unity vector
    m_mapcorr.assign(m_energy.size(), 1.0);

    // Get pointer on event cube
    const GLATEventCube* cube = dynamic_cast<const GLATEventCube*>(obs.events());

    // Continue only if we have an event cube
    if (cube != NULL) {

        // Get maximum PSF radius
        double radius = cube->maxrad(m_dir);
        if (radius > 0.0) {

            // Initialise sum over pixels
            std::vector<double> sum(m_energy.size(), 0.0);

            // Loop over all event cube spatial pixels
            for (int iy = 0; iy < cube->ny(); ++iy) {
                for (int ix = 0; ix < cube->nx(); ++ix) {

                    // Compute offset angle in degrees
                    GSkyPixel pixel  = GSkyPixel(double(ix), double(iy));
                    double    offset = cube->map().pix2dir(pixel).dist_deg(m_dir);
                    double    omega  = cube->map().omega(pixel);

                    // Use only pixels within maximum PSF radius
                    if (offset <= radius) {

                        // Accumulate energy dependent pixel sum
                        for (int ieng = 0; ieng < m_energy.size(); ++ieng)
                            sum[ieng] += psf(offset, m_energy[ieng]) * omega;
                        
                    } // endif: pixel was within maximum PSF radius

                } // endfor: looped over latitude pixels
            } // endfor: looped over longitude pixels

            // Compute map correction
            for (int ieng = 0; ieng < m_energy.size(); ++ieng) {
                if (sum[ieng] > 0.0) {
                    m_mapcorr[ieng] = integral(radius, m_energy[ieng]) / sum[ieng];
                }
            }

        } // endif: radius was positive
    } // endif: observation had an event cube

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute integral over PSF
 *
 * @param[in] offsetmax Maximum offset angle.
 * @param[in] logE log10 of energy in MeV.
 ***************************************************************************/
double GLATMeanPsf::integral(const double& offsetmax, const double& logE)
{
    // Set energy and offset interpolation weigths
    m_energy.set_value(logE);
    m_offset.set_value(offsetmax);

    // Get PSF array offsets
    int inx_energy_left  = m_energy.inx_left()  * noffsets();
    int inx_energy_right = m_energy.inx_right() * noffsets();

    // Initialise integrals
    double int_left  = 0.0;
    double int_right = 0.0;
    
    // Loop over all offset angles
    for (int i = 0; i < m_offset.size()-1; ++i) {

        // If interval is fully contained within offset then simply use
        // tabulated values for integration
        if (m_offset[i+1] <= offsetmax) {
            double theta_min = m_offset[i]   * gammalib::deg2rad;
            double theta_max = m_offset[i+1] * gammalib::deg2rad;
            int_left  += 0.5 * (m_psf[inx_energy_left+i]   * sin(theta_min) +
                                m_psf[inx_energy_left+i+1] * sin(theta_max)) *
                               (theta_max - theta_min);
            int_right += 0.5 * (m_psf[inx_energy_right+i]   * sin(theta_min) +
                                m_psf[inx_energy_right+i+1] * sin(theta_max)) *
                               (theta_max - theta_min);
        }

        // ... otherwise interpolate the PSF value for theta_max. We can
        // then exit the loop since we're done with the integration
        else {
            double theta_min = m_offset[i] * gammalib::deg2rad;
            double theta_max = offsetmax   * gammalib::deg2rad;
            double psf_left  = m_psf[inx_energy_left+i]    * m_offset.wgt_left() +
                               m_psf[inx_energy_left+i+1]  * m_offset.wgt_right();
            double psf_right = m_psf[inx_energy_right+i]   * m_offset.wgt_left() +
                               m_psf[inx_energy_right+i+1] * m_offset.wgt_right();                               
            int_left  += 0.5 * (m_psf[inx_energy_left+i] * sin(theta_min) +
                                psf_left                 * sin(theta_max)) *
                               (theta_max - theta_min);
            int_right += 0.5 * (m_psf[inx_energy_right+i] * sin(theta_min) +
                                psf_right                 * sin(theta_max)) *
                               (theta_max - theta_min);

            // Debug option: Dump integral computation results
            #if G_DEBUG_INTEGRAL
            std::cout << "offsetmax=" << offsetmax;
            std::cout << " offset.left=" << m_offset.inx_left();
            std::cout << " offset.right=" << m_offset.inx_right();
            std::cout << " i=" << i;
            std::cout << " psf_left=" << psf_left;
            std::cout << " psf_right=" << psf_right;
            std::cout << " int_left=" << int_left*gammalib::twopi;
            std::cout << " int_right=" << int_right*gammalib::twopi;
            #endif
            break;
        }

    } // endfor: looped over offset angles

    // Interpolate now in energy
    double integral = gammalib::twopi * (m_energy.wgt_left()  * int_left +
                                         m_energy.wgt_right() * int_right);

    // Debug option: Dump integral computation results
    #if G_DEBUG_INTEGRAL
    std::cout << " integral=" << integral << std::endl;
    #endif

    // Return integral
    return integral;
}
