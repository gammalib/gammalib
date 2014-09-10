/***************************************************************************
 *                  GLATAeff.hpp - Fermi/LAT effective area                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GLATAeff.hpp
 * @brief Fermi/LAT effective area class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATAEFF_HPP
#define GLATAEFF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GLATResponseTable.hpp"
#include "GLATEfficiency.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GLATAeff
 *
 * @brief Interface for the Fermi/LAT effective area
 *
 * This class handles the effective area information for Fermi/LAT. It also
 * handles the IRF efficiency information that has been introduced for
 * Pass 7 data analysis.
 ***************************************************************************/
class GLATAeff : public GBase {

public:
    // Constructors and destructors
    GLATAeff(void);
    explicit GLATAeff(const std::string& filename);
    GLATAeff(const GLATAeff& aeff);
    virtual ~GLATAeff(void);

    // Operators
    GLATAeff& operator=(const GLATAeff& aeff);
    double    operator()(const double& logE, const double& ctheta);
    double    operator()(const double& logE, const double& ctheta,
                         const double& phi);

    // Methods
    void          clear(void);
    GLATAeff*     clone(void) const;
    std::string   classname(void) const;
    void          load(const std::string& filename);
    void          save(const std::string& filename,
                       const bool &clobber = false);
    void          read(const GFits& file);
    void          write(GFits& file) const;
    int           size(void) const;
    int           nenergies(void) const;
    int           ncostheta(void) const;
    const double& costhetamin(void) const;
    void          costhetamin(const double& ctheta);
    bool          has_phi(void) const;
    bool          has_efficiency(void) const;
    double        efficiency_factor1(const GEnergy& srcEng) const;
    double        efficiency_factor2(const GEnergy& srcEng) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATAeff& aeff);
    void free_members(void);
    void read_aeff(const GFitsTable& hdu);
    void read_efficiency(const GFitsTable& hdu);
    void write_aeff(GFits& file) const;
    void write_efficiency(GFits& file) const;
    
    // Protected members
    GLATResponseTable   m_aeff_bins;    //!< Aeff energy and cos theta binning
    std::vector<double> m_aeff;         //!< Aeff array
    double              m_min_ctheta;   //!< Minimum valid cos(theta)
    bool                m_front;        //!< Response is for front section
    bool                m_back;         //!< Response is for back section
    GLATEfficiency*     m_eff_func1;    //!< Efficiency functor 1
    GLATEfficiency*     m_eff_func2;    //!< Efficiency functor 2
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATAeff").
 ***************************************************************************/
inline
std::string GLATAeff::classname(void) const
{
    return ("GLATAeff");
}


/***********************************************************************//**
 * @brief Return number of bins in effective area response
 *
 * @return Number of bins in effective area response.
 ***************************************************************************/
inline
int GLATAeff::size(void) const
{
    return nenergies()*ncostheta();
}


/***********************************************************************//**
 * @brief Return number of energies in effective area response
 *
 * @return Number of energies in effective area response.
 ***************************************************************************/
inline
int GLATAeff::nenergies(void) const
{
    return m_aeff_bins.nenergies();
}


/***********************************************************************//**
 * @brief Return number of cosine theta bins in effective area response
 *
 * @return Number of cosine theta bins in effective area response.
 ***************************************************************************/
inline
int GLATAeff::ncostheta(void) const
{
    return m_aeff_bins.ncostheta();
}


/***********************************************************************//**
 * @brief Return cosine theta minimum
 *
 * @return Cosine theta minimum.
 ***************************************************************************/
inline
const double& GLATAeff::costhetamin(void) const
{
    return m_min_ctheta;
}


/***********************************************************************//**
 * @brief Set minimum cos(theta) angle for effective area access
 *
 * @param[in] ctheta Cosine of maximum zenith angle.
 ***************************************************************************/
inline
void GLATAeff::costhetamin(const double& ctheta)
{
    m_min_ctheta = ctheta;
    return;
}


/***********************************************************************//**
 * @brief Signal that effective area has Phi dependence
 *
 * @return True if effective area has Phi dependence.
 ***************************************************************************/
inline
bool GLATAeff::has_phi(void) const
{
    return false;
}

#endif /* GLATAEFF_HPP */
