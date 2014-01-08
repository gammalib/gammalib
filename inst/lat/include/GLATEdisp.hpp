/***************************************************************************
 *              GLATEdisp.hpp  -  Fermi LAT energy dispersion              *
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
#include "GFits.hpp"
#include "GFitsTable.hpp"


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
    explicit GLATEdisp(const std::string& filename);
    GLATEdisp(const GLATEdisp& edisp);
    virtual ~GLATEdisp(void);

    // Operators
    GLATEdisp& operator=(const GLATEdisp& edisp);
    //double    operator()(const double& logE, const double& ctheta);
    //double    operator()(const GSkyDir& srcDir, const GEnergy& srcEng,
    //                     const GTime& srcTime, const GLATPointing& pnt);

    // Methods
    void         clear(void);
    GLATEdisp*   clone(void) const;
    void         load(const std::string& filename);
    void         save(const std::string& filename,
                      const bool& clobber = false);
    void         read(const GFits& file);
    void         write(GFits& file) const;
    int          size(void) const;
    int          nenergies(void) const;
    int          ncostheta(void) const;
    //double       costhetamin(void) const { return m_min_ctheta; }
    //void         costhetamin(const double& ctheta);
    bool         has_phi(void) const { return false; }
    std::string  print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATEdisp& edisp);
    void free_members(void);
    void read_edisp(const GFitsTable& hdu);
    void write_edisp(GFits& file) const;
    
    // Protected members
    GLATResponseTable   m_edisp_bins;   //!< Energy dispersion energy and cos theta binning
    std::vector<double> m_norm;         //!< Energy dispersion normalization
    std::vector<double> m_ls1;          //!< Energy dispersion ...
    std::vector<double> m_scale;        //!< Energy dispersion scaling parameters
};


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

#endif /* GLATEDISP_HPP */
