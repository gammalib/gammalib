/***************************************************************************
 *      GCTAPsfKing.hpp - King profile CTA point spread function class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2017 by Michael Mayer                               *
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
 * @file GCTAPsfKing.hpp
 * @brief King profile CTA point spread function class definition
 * @author Michael Mayer
 */

#ifndef GCTAPsfKing_HPP
#define GCTAPsfKing_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFilename.hpp"
#include "GCTAPsf.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GFitsTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_cta_psfking = "POINT SPREAD FUNCTION";
}


/***********************************************************************//**
 * @class GCTAPsfKing
 *
 * @brief CTA point spread function class with a King profile
 *
 * This class implements the CTA point spread function response as function
 * of energy as determined from a FITS table.
 ***************************************************************************/
class GCTAPsfKing : public GCTAPsf {

public:
    // Constructors and destructors
    GCTAPsfKing(void);
    explicit GCTAPsfKing(const GFilename& filename);
    GCTAPsfKing(const GCTAPsfKing& psf);
    virtual ~GCTAPsfKing(void);

    // Operators
    GCTAPsfKing& operator=(const GCTAPsfKing& psf);
    double operator()(const double& delta,
                      const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void         clear(void);
    GCTAPsfKing* clone(void) const;
    std::string  classname(void) const;
    void         load(const GFilename& filename);
    GFilename    filename(void) const;
    double       mc(GRan&         ran,
                    const double& logE, 
                    const double& theta = 0.0, 
                    const double& phi = 0.0,
                    const double& zenith = 0.0,
                    const double& azimuth = 0.0,
                    const bool&   etrue = true) const;
    double       delta_max(const double& logE, 
                           const double& theta = 0.0, 
                           const double& phi = 0.0,
                           const double& zenith = 0.0,
                           const double& azimuth = 0.0,
                           const bool&   etrue = true) const;
    double       containment_radius(const double& fraction,
                                    const double& logE,
                                    const double& theta = 0.0,
                                    const double& phi = 0.0,
                                    const double& zenith = 0.0,
                                    const double& azimuth = 0.0,
                                    const bool&   etrue = true) const;

    std::string  print(const GChatter& chatter = NORMAL) const;

    // Methods
    const GCTAResponseTable&   table(void) const;
    void                       table(const GCTAResponseTable& table);
    void                       read(const GFitsTable& table);
    void                       write(GFitsBinTable& table) const;
    void                       save(const GFilename& filename,
                                    const bool&      clobber = false) const;
  


private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAPsfKing& psf);
    void free_members(void);
    void update(const double& logE, const double& theta) const;

    // Members
    GFilename         m_filename;   //!< Name of Aeff response file
    GCTAResponseTable m_psf;        //!< PSF response table
    int               m_inx_energy; //!< Energy index
    int               m_inx_theta;  //!< Theta index
    int               m_inx_gamma;  //!< Gamma
    int               m_inx_sigma;  //!< Sigma

    // Evaluation cache
    mutable double m_par_logE;   //!< Cache energy
    mutable double m_par_theta;  //!< Cache offset angle
    mutable double m_par_norm;   //!< King profile normalization
    mutable double m_par_sigma;  //!< King profile sigma (radians)
    mutable double m_par_sigma2; //!< King profile sigma squared
    mutable double m_par_gamma;  //!< King profile gamma parameter
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAPsfKing").
 ***************************************************************************/
inline
std::string GCTAPsfKing::classname(void) const
{
    return ("GCTAPsfKing");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Filename from which point spread function was loaded.
 *
 * Returns filename from which point spread function was loaded.
 ***************************************************************************/
inline
GFilename GCTAPsfKing::filename(void) const
{
    return m_filename;
}

/***********************************************************************//**
 * @brief Return response table
 *
 * @return Response table.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTAPsfKing::table(void) const
{
    return m_psf;
}

/***********************************************************************//**
 * @brief Assign response table
 *
 * @param[in] table Response table.
 ***************************************************************************/
inline
void GCTAPsfKing::table(const GCTAResponseTable& table)
{
     m_psf = table;
}


#endif /* GCTAPsfKing_HPP */
