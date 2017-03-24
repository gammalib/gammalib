/***************************************************************************
 *                 GCTAAeff2D.hpp - CTA 2D effective area class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file GCTAAeff2D.hpp
 * @brief CTA 2D effective area class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAAEFF2D_HPP
#define GCTAAEFF2D_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GFilename.hpp"
#include "GCTAAeff.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFitsTable;
class GFitsBinTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_cta_aeff2d = "EFFECTIVE AREA";
}


/***********************************************************************//**
 * @class GCTAAeff2D
 *
 * @brief CTA 2D effective area class
 *
 * This class implements the CTA effective area response as function of
 * energy and offset angle.
 ***************************************************************************/
class GCTAAeff2D : public GCTAAeff {

public:
    // Constructors and destructors
    GCTAAeff2D(void);
    explicit GCTAAeff2D(const GFilename& filename);
    GCTAAeff2D(const GCTAAeff2D& cta);
    virtual ~GCTAAeff2D(void);

    // Operators
    GCTAAeff2D& operator=(const GCTAAeff2D& aeff);
    double operator()(const double& logE, 
                      const double& theta = 0.0, 
                      const double& phi = 0.0,
                      const double& zenith = 0.0,
                      const double& azimuth = 0.0,
                      const bool&   etrue = true) const;

    // Implemented pure virtual methods
    void        clear(void);
    GCTAAeff2D* clone(void) const;
    std::string classname(void) const;
    void        load(const GFilename& filename);
    GFilename   filename(void) const;
    double      max(const double& logE,
                    const double& zenith,
                    const double& azimuth,
                    const bool&   etrue = true) const;
    std::string print(const GChatter& chatter = NORMAL) const;

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
    void copy_members(const GCTAAeff2D& aeff);
    void free_members(void);
    void set_indices(void);
    void set_boundaries(void);

    // Members
    GFilename         m_filename;      //!< Name of Aeff response file
    GCTAResponseTable m_aeff;          //!< Aeff response table
    int               m_inx_energy;    //!< Energy index
    int               m_inx_theta;     //!< Theta index
    int               m_inx_aeff;      //!< Effective area (true energy)
    int               m_inx_aeff_reco; //!< Effective area (reconstructed energy)
    double            m_logE_min;      //!< Minimum logE (log10(E/TeV))
    double            m_logE_max;      //!< Maximum logE (log10(E/TeV))
    double            m_theta_min;     //!< Minimum theta (radians)
    double            m_theta_max;     //!< Maximum theta (radians)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAAeff2D").
 ***************************************************************************/
inline
std::string GCTAAeff2D::classname(void) const
{
    return ("GCTAAeff2D");
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which effective area was loaded
 ***************************************************************************/
inline
GFilename GCTAAeff2D::filename(void) const
{
    return m_filename;
}


/***********************************************************************//**
 * @brief Return response table
 *
 * @return Response table.
 *
 * Returns the response table of the effective area. The effective area
 * values are given in units of cm2.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTAAeff2D::table(void) const
{
    return m_aeff;
}

#endif /* GCTAAEFF2D_HPP */
