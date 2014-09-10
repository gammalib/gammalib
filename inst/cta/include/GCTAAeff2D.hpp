/***************************************************************************
 *                 GCTAAeff2D.hpp - CTA 2D effective area class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
#include "GFits.hpp"
#include "GCTAAeff.hpp"
#include "GCTAResponseTable.hpp"


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
    explicit GCTAAeff2D(const std::string& filename);
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
    void        load(const std::string& filename);
    std::string filename(void) const;
    std::string print(const GChatter& chatter = NORMAL) const;

    // Methods
    const GCTAResponseTable& table(void) const;
    void                     table(const GCTAResponseTable& table);
    void                     read(const GFits& file);
    void                     write(GFitsBinTable& hdu) const;
    void                     save(const std::string& filename,
                                  const bool& clobber = false) const;
    
private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAAeff2D& aeff);
    void free_members(void);

    // Members
    std::string       m_filename;  //!< Name of Aeff response file
    GCTAResponseTable m_aeff;      //!< Aeff response table
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
std::string GCTAAeff2D::filename(void) const
{
    // Return filename
    return m_filename;
}


/***********************************************************************//**
 * @brief Return response table
 *
 * @return Response table.
 ***************************************************************************/
inline
const GCTAResponseTable& GCTAAeff2D::table(void) const
{
    return m_aeff;
}


/***********************************************************************//**
 * @brief Assign response table
 *
 * @param[in] table Response table.
 ***************************************************************************/
inline
void GCTAAeff2D::table(const GCTAResponseTable& table)
{
     m_aeff = table;
}

#endif /* GCTAAEFF2D_HPP */
