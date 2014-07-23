/***************************************************************************
 *        GCTASourceCubeDiffuse.hpp - CTA diffuse source cube class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTASourceCubeDiffuse.hpp
 * @brief CTA diffuse source cube class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTASOURCECUBEDIFFUSE_HPP
#define GCTASOURCECUBEDIFFUSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GCTASourceCube.hpp"
#include "GSkymap.hpp"
#include "GFunction.hpp"
#include "GCTAResponseCube.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declarations _______________________________________________ */
class GModelSpatial;
class GObservation;
class GModelSpatialDiffuse;
class GSkyDir;
class GEnergy;
class GTime;


/***********************************************************************//**
 * @class GCTASourceCubeDiffuse
 *
 * @brief CTA diffuse source cube class
 *
 * This class handles pre-computed response information for a diffuse source
 * in a stacked cube analysis. It derives from the abstract GCTASourceCube
 * class.
 ***************************************************************************/
class GCTASourceCubeDiffuse : public GCTASourceCube {

public:
    // Constructors and destructors
    GCTASourceCubeDiffuse(void);
    GCTASourceCubeDiffuse(const GCTASourceCubeDiffuse& cube);
    virtual ~GCTASourceCubeDiffuse(void);

    // Operators
    GCTASourceCubeDiffuse& operator=(const GCTASourceCubeDiffuse & cube);

    // Implemented pure virtual methods
    void                   clear(void);
    GCTASourceCubeDiffuse* clone(void) const;
    virtual GCTAClassCode  code(void) const;
    void                   set(const std::string&   name,
                               const GModelSpatial& model,
                               const GObservation&  obs);
    std::string            print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double irf(const int& pixel, const int& iebin) const;
    double psf(const GCTAResponseCube*     rsp,
               const GModelSpatialDiffuse* model,
               const GSkyDir&              srcDir,
               const GEnergy&              srcEng,
               const GTime&                srcTime) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTASourceCubeDiffuse& cube);
    void free_members(void);

    // PSF delta integration kernel
    class psf_kern_delta : public GFunction {
    public:
        psf_kern_delta(const GCTAResponseCube*     rsp,
                       const GModelSpatialDiffuse* model,
                       const GSkyDir&              srcDir,
                       const GEnergy&              srcEng,
                       const GTime&                srcTime,
                       const GMatrix&              rot) :
                       m_rsp(rsp),
                       m_model(model),
                       m_srcDir(srcDir),
                       m_srcEng(srcEng),
                       m_srcTime(srcTime),
                       m_rot(rot),
                       m_psf_max(rsp->psf()(srcDir, 0.0, srcEng)) { }
        double eval(const double& delta);
    protected:
        const GCTAResponseCube*     m_rsp;     //!< Response cube
        const GModelSpatialDiffuse* m_model;   //!< Spatial model
        const GSkyDir&              m_srcDir;  //!< True photon arrival direction
        const GEnergy&              m_srcEng;  //!< True photon energy
        const GTime&                m_srcTime; //!< True photon arrival time
        const GMatrix&              m_rot;     //!< Rotation matrix
        double                      m_psf_max; //!< Maximum PSF value
    };

    // PSF phi integration kernel
    class psf_kern_phi : public GFunction {
    public:
        psf_kern_phi(const GModelSpatialDiffuse* model,
                     const GEnergy&              srcEng,
                     const GTime&                srcTime,
                     const GMatrix&              rot,
                     const double&               sin_delta,
                     const double&               cos_delta) :
                     m_model(model),
                     m_srcEng(srcEng),
                     m_srcTime(srcTime),
                     m_rot(rot),
                     m_sin_delta(sin_delta),
                     m_cos_delta(cos_delta) { }
        double eval(const double& phi);
    protected:
        const GModelSpatialDiffuse* m_model;     //!< Spatial model
        const GEnergy&              m_srcEng;    //!< True photon energy
        const GTime&                m_srcTime;   //!< True photon arrival time
        const GMatrix&              m_rot;       //!< Rotation matrix
        const double&               m_sin_delta; //!< sin(delta)
        const double&               m_cos_delta; //!< cos(delta)
    };

    // Data members
    GSkymap m_cube;  //!< Diffuse map convolved with IRF
};


/***********************************************************************//**
 * @brief Return class type code
 *
 * @return GCTA_SOURCE_CUBE_DIFFUSE.
 *
 * Returns the class type code GCTA_SOURCE_CUBE_DIFFUSE.
 ***************************************************************************/
inline
GCTAClassCode GCTASourceCubeDiffuse::code(void) const
{
    return (GCTA_SOURCE_CUBE_DIFFUSE);
}


/***********************************************************************//**
 * @brief Return instrument response function
 *
 * @param[in] pixel Spatial pixel index [0,...,???]
 * @param[in] iebin Energy index [0,...,???]
 * @return Instrument response function.
 *
 * Returns the instrument response function.
 ***************************************************************************/
inline
double GCTASourceCubeDiffuse::irf(const int& pixel, const int& iebin) const
{
    return (m_cube(pixel, iebin));
}

#endif /* GCTASOURCECUBEDIFFUSE_HPP */
