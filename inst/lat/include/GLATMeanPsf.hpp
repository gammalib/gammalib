/***************************************************************************
 *                GLATMeanPsf.hpp - Fermi/LAT mean PSF class               *
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
 * @file GLATMeanPsf.hpp
 * @brief Fermi/LAT mean PSF class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATMEANPSF_HPP
#define GLATMEANPSF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GSkyDir.hpp"
#include "GNodeArray.hpp"

/* __ Forward declarations _______________________________________________ */
class GLATObservation;


/***********************************************************************//**
 * @class GLATMeanPsf
 *
 * @brief Fermi/LAT mean PSF class
 *
 * The position-dependent mean PSF is the point spread function that has
 * been averaged over the zenith and azimuth angles of an observation. The
 * averaging is done using the livetime cube which holds the lifetime as
 * function and zenith and azimuth angles for an observation.
 ***************************************************************************/
class GLATMeanPsf : public GBase {

public:
    // Constructors and destructors
    GLATMeanPsf(void);
    GLATMeanPsf(const GSkyDir& dir, const GLATObservation& obs);
    GLATMeanPsf(const GLATMeanPsf& cube);
    virtual ~GLATMeanPsf(void);

    // Operators
    GLATMeanPsf& operator=(const GLATMeanPsf& cube);
    double       operator()(const double& offset, const double& logE);

    // Methods
    void               clear(void);
    GLATMeanPsf*       clone(void) const;
    int                size(void) const;
    void               set(const GSkyDir& dir, const GLATObservation& obs);
    int                noffsets(void) const;
    int                nenergies(void) const;
    const double&      offset(const int& inx) const;
    const double&      energy(const int& inx) const;
    const GSkyDir&     dir(void) const;
    const std::string& name(void) const;
    void               name(const std::string& name);
    const double&      thetamax(void) const;
    void               thetamax(const double& value);
    double             psf(const double& offset, const double& logE);
    double             exposure(const double& logE);
    std::string        print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void   init_members(void);
    void   copy_members(const GLATMeanPsf& psf);
    void   free_members(void);
    void   set_offsets(void);
    void   set_map_corrections(const GLATObservation& obs);
    double integral(const double& radmax, const double& logE);
    
    // Protected members
    std::string          m_name;         //!< Source name for mean PSF
    GSkyDir              m_dir;          //!< Source direction for mean PSF
    std::vector<double>  m_psf;          //!< Mean PSF values
    std::vector<double>  m_exposure;     //!< Mean exposure
    std::vector<double>  m_mapcorr;      //!< Map corrections
    GNodeArray           m_offset;       //!< Offsets of mean PSF
    GNodeArray           m_energy;       //!< log10(energy) of mean PSF
    double               m_theta_max;    //!< Maximum inclination angle (default 70 deg)

    // Bi-linear interpolation data
    double               m_last_energy;  //!< Last requested logE value
    double               m_last_offset;  //!< Last requested offset value
    int                  m_inx1_exp;     //!< Exposure index 1
    int                  m_inx2_exp;     //!< Exposure index 2
    int                  m_inx1;         //!< Index 1
    int                  m_inx2;         //!< Index 2
    int                  m_inx3;         //!< Index 3
    int                  m_inx4;         //!< Index 4
    double               m_wgt1;         //!< Weighting factor 1
    double               m_wgt2;         //!< Weighting factor 2
    double               m_wgt3;         //!< Weighting factor 3
    double               m_wgt4;         //!< Weighting factor 4
};


/***********************************************************************//**
 * @brief Return number of bins in mean PSF
 *
 * @return Number of bins in mean PSF
 ***************************************************************************/
inline
int GLATMeanPsf::size(void) const
{
    // Compute size
    int size = m_energy.size()*m_offset.size();

    // Return size
    return size;
}


/***********************************************************************//**
 * @brief Return number of offset bins
 *
 * @return Number of offset bins.
 ***************************************************************************/
inline
int GLATMeanPsf::noffsets(void) const
{
    return m_offset.size();
}


/***********************************************************************//**
 * @brief Return number of energy bins
 *
 * @return Number of energy bins.
 ***************************************************************************/
inline
int GLATMeanPsf::nenergies(void) const
{
    return m_energy.size();
}


/***********************************************************************//**
 * @brief Return offset angle for given bin
 *
 * @param[in] inx Bin index [0,...,noffsets()-1]
 * @return Offset angle.
 ***************************************************************************/
inline
const double& GLATMeanPsf::offset(const int& inx) const
{
    return m_offset[inx];
}


/***********************************************************************//**
 * @brief Return energy for given bin
 *
 * @param[in] inx Bin index [0,...,nenergies()-1]
 * @return Energy.
 ***************************************************************************/
inline
const double& GLATMeanPsf::energy(const int& inx) const
{
    return m_energy[inx];
}


/***********************************************************************//**
 * @brief Return sky direction for mean PSF
 *
 * @return Sky direction.
 ***************************************************************************/
inline
const GSkyDir& GLATMeanPsf::dir(void) const
{
    return m_dir;
}


/***********************************************************************//**
 * @brief Return source name for mean PSF
 *
 * @return Source name.
 ***************************************************************************/
inline
const std::string& GLATMeanPsf::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Set source name for mean PSF
 *
 * @param[in] name Source name.
 ***************************************************************************/
inline
void GLATMeanPsf::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Return maximum theta angle
 *
 * @return Maximum theta angle.
 ***************************************************************************/
inline
const double& GLATMeanPsf::thetamax(void) const
{
    return m_theta_max;
}


/***********************************************************************//**
 * @brief Set maximum theta angle
 *
 * @param[in] value Maximum theta angle.
 ***************************************************************************/
inline
void GLATMeanPsf::thetamax(const double& value)
{
    m_theta_max = value;
    return;
}

#endif /* GLATMEANPSF_HPP */
