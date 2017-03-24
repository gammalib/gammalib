/***************************************************************************
 *              GLATEdisp.hpp  -  Fermi LAT energy dispersion              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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
 * @file GLATEdisp.hpp
 * @brief Fermi LAT energy dispersion class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GLATEDISP_HPP
#define GLATEDISP_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GLATResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFits;
class GFitsTable;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_lat_edisp       = "ENERGY DISPERSION";
    const std::string extname_lat_edisp_scale = "EDISP_SCALING_PARAMS";
}


/***********************************************************************//**
 * @class GLATEdisp
 *
 * @brief Interface for the Fermi LAT energy dispersion.
 *
 * @todo Implement support for older response functions?
 ***************************************************************************/
class GLATEdisp : public GBase {

public:
    // Constructors and destructors
    GLATEdisp(void);
    GLATEdisp(const GFilename& filename, const std::string& evtype);
    GLATEdisp(const GLATEdisp& edisp);
    virtual ~GLATEdisp(void);

    // Operators
    GLATEdisp& operator=(const GLATEdisp& edisp);
    //double    operator()(const double& logE, const double& ctheta);
    //double    operator()(const GSkyDir& srcDir, const GEnergy& srcEng,
    //                     const GTime& srcTime, const GLATPointing& pnt);

    // Methods
    void               clear(void);
    GLATEdisp*         clone(void) const;
    std::string        classname(void) const;
    const std::string& evtype(void) const;
    void               load(const GFilename&   filename,
                            const std::string& evtype);
    void               save(const GFilename& filename,
                            const bool&      clobber = false);
    void               read(const GFits& file);
    void               write(GFits& file) const;
    int                size(void) const;
    int                nenergies(void) const;
    int                ncostheta(void) const;
    //double             costhetamin(void) const { return m_min_ctheta; }
    //void               costhetamin(const double& ctheta);
    bool               has_phi(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATEdisp& edisp);
    void free_members(void);
    void read_edisp(const GFitsTable& hdu);
    void write_edisp(GFits& file) const;
    
    // Protected members
    std::string         m_evtype;     //!< Event type
    GLATResponseTable   m_edisp_bins; //!< Energy dispersion energy and cos theta binning
    std::vector<double> m_norm;       //!< Energy dispersion normalization
    std::vector<double> m_ls1;        //!< Energy dispersion ...
    std::vector<double> m_scale;      //!< Energy dispersion scaling parameters
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATEdisp").
 ***************************************************************************/
inline
std::string GLATEdisp::classname(void) const
{
    return ("GLATEdisp");
}


/***********************************************************************//**
 * @brief Return event type
 *
 * @return Event type.
 ***************************************************************************/
inline
const std::string& GLATEdisp::evtype(void) const
{
    return (m_evtype);
}


/***********************************************************************//**
 * @brief Return number of bins in energy dispersion response
 *
 * @return Number of bins in energy dispersion response.
 ***************************************************************************/
inline
int GLATEdisp::size(void) const
{
    return nenergies()*ncostheta();
}


/***********************************************************************//**
 * @brief Return number of energies in energy dispersion response
 *
 * @return Number of energies in energy dispersion response.
 ***************************************************************************/
inline
int GLATEdisp::nenergies(void) const
{
    return m_edisp_bins.nenergies();
}


/***********************************************************************//**
 * @brief Return number of cosine theta bins in energy dispersion response
 *
 * @return Number of cosine theta bins in energy dispersion response.
 ***************************************************************************/
inline
int GLATEdisp::ncostheta(void) const
{
    return m_edisp_bins.ncostheta();
}


/***********************************************************************//**
 * @brief Signal that energy dispersion has Phi dependence
 *
 * @return True if energy dispersion has Phi dependence.
 ***************************************************************************/
inline
bool GLATEdisp::has_phi(void) const
{
    return false;
}

#endif /* GLATEDISP_HPP */
