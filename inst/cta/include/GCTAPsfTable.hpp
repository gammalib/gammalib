/***************************************************************************
 *         GCTAPsfTable.hpp - CTA point spread function table class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2017 by Juergen Knoedlseder                         *
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
 * @file GCTAPsfTable.hpp
 * @brief CTA point spread function table class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAPSFTABLE_HPP
#define GCTAPSFTABLE_HPP

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
    const std::string extname_cta_psftable = "POINT SPREAD FUNCTION";
}


/***********************************************************************//**
 * @class GCTAPsfTable
 *
 * @brief CTA point spread function table class
 *
 * This class implements the CTA point spread function response as function
 * of energy as determined from a FITS table.
 ***************************************************************************/
class GCTAPsfTable : public GCTAPsf {

public:
    // Constructors and destructors
    GCTAPsfTable(void);
    explicit GCTAPsfTable(const GFilename& filename);
    GCTAPsfTable(const GCTAPsfTable& psf);
    virtual ~GCTAPsfTable(void);

    // Operators
    GCTAPsfTable& operator=(const GCTAPsfTable& psf);
    double operator()(const double& delta,
                      const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void          clear(void);
    GCTAPsfTable* clone(void) const;
    std::string   classname(void) const;
    void          load(const GFilename& filename);
    GFilename     filename(void) const;
    double        mc(GRan&         ran,
                     const double& logE,
                     const double& theta = 0.0,
                     const double& phi = 0.0,
                     const double& zenith = 0.0,
                     const double& azimuth = 0.0,
                     const bool&   etrue = true) const;
    double        delta_max(const double& logE,
                            const double& theta = 0.0,
                            const double& phi = 0.0,
                            const double& zenith = 0.0,
                            const double& azimuth = 0.0,
                            const bool&   etrue = true) const;
    double        containment_radius(const double& fraction,
                                     const double& logE,
                                     const double& theta = 0.0,
                                     const double& phi = 0.0,
                                     const double& zenith = 0.0,
                                     const double& azimuth = 0.0,
                                     const bool&   etrue = true) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

    // Methods
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFitsTable& table);
    void                     write(GFitsBinTable& table) const;
    void                     save(const GFilename& filename,
                                  const bool&      clobber = false) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAPsfTable& psf);
    void free_members(void);
    void precompute(void);
    int  element(const int& ieng, const int& itheta, const int& idelta);

    // Members
    GFilename         m_filename;   //!< Name of Aeff response file
    GCTAResponseTable m_psf;        //!< PSF response table
    int               m_inx_energy; //!< Energy index
    int               m_inx_theta;  //!< Theta index
    int               m_inx_delta;  //!< Delta index
    int               m_inx_rpsf;   //!< PSF histogram
    double            m_delta_max;  //!< Maximum delta angle (radians)
    double            m_psf_max;    //!< Maximum PSF value
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAPsfTable").
 ***************************************************************************/
inline
std::string GCTAPsfTable::classname(void) const
{
    return ("GCTAPsfTable");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Filename from which point spread function was loaded.
 *
 * Returns filename from which point spread function was loaded.
 ***************************************************************************/
inline
GFilename GCTAPsfTable::filename(void) const
{
    return m_filename;
}

/***********************************************************************//**
 * @brief Return response table
 *
 * @return Response table.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTAPsfTable::table(void) const
{
    return m_psf;
}

/***********************************************************************//**
 * @brief Assign response table
 *
 * @param[in] table Response table.
 ***************************************************************************/
inline
void GCTAPsfTable::table(const GCTAResponseTable& table)
{
     m_psf = table;
}

#endif /* GCTAPSFTABLE_HPP */
