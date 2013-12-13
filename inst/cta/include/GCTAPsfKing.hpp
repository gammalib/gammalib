/***************************************************************************
 *       GCTAPsfKing.hpp - CTA point spread function vector class        *
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
 * @file GCTAPsfKing.hpp
 * @brief CTA point spread function class definition with a King profile
 * @author Michael Mayer
 */

#ifndef GCTAPsfKing_HPP
#define GCTAPsfKing_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFits.hpp"
#include "GRan.hpp"
#include "GCTAPsf.hpp"
#include "GCTAResponseTable.hpp"


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
    GCTAPsfKing(const std::string& filename);
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
    void           clear(void);
    GCTAPsfKing* clone(void) const;
    void           load(const std::string& filename);
    std::string    filename(void) const;
    double         mc(GRan&         ran,
                      const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;
    double         delta_max(const double& logE, 
                             const double& theta = 0.0, 
                             const double& phi = 0.0,
                             const double& zenith = 0.0,
                             const double& azimuth = 0.0,
                             const bool&   etrue = true) const;
    std::string print(const GChatter& chatter = NORMAL) const;


private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAPsfKing& psf);
    void free_members(void);
    void update(const double& logE, const double& theta) const;

    // Members
   std::string       m_filename;   //!< Name of Aeff response file
   GCTAResponseTable m_psf;        //!< PSF response table

   mutable double m_par_logE; //!< Cache energy
   mutable double m_par_theta; //!< Cache offset angle
   mutable double m_par_norm; //!< King profile normalization
   mutable double m_par_sigma; //!< King profile sigma (radians)
   mutable double m_par_gamma; //!< King profile gamma parameter

};

#endif /* GCTAPsfKing_HPP */
