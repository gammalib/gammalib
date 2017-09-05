/***************************************************************************
 *                  GCOMDri.hpp - COMPTEL Data Space class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMDri.hpp
 * @brief COMPTEL Data Space class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMDRI_HPP
#define GCOMDRI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GSkyMap.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFits;
class GFitsImage;

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const std::string extname_dri = "DRI";
}


/***********************************************************************//**
 * @class GCOMDri
 *
 * @brief COMPTEL Data Space class
 ***************************************************************************/
class GCOMDri : public GBase {

public:
    // Constructors and destructors
    GCOMDri(void);
    explicit GCOMDri(const GFilename& filename);
    GCOMDri(const GCOMDri& dri);
    virtual ~GCOMDri(void);

    // Operators
    GCOMDri&      operator=(const GCOMDri& dri);
    double&       operator[](const int& index);
    const double& operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMDri*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    int             size(void) const;
    int             nchi(void) const;
    int             npsi(void) const;
    int             nphibar(void) const;
    void            load(const GFilename& filename);
    void            save(const GFilename& filename, const bool& clobber = false) const;
    void            read(const GFitsImage& image);
    void            write(GFits& fits, const std::string& extname = "") const;
    const GSkyMap&  map(void) const;
    const GEbounds& ebounds(void) const;
    const GGti&     gti(void) const;
    const double&   phimin(void) const;
    const double&   phibin(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMDri& dri);
    void free_members(void);
    void read_attributes(const GFitsHDU* hdu);
    void write_attributes(GFitsHDU* hdu) const;

    // Protected members
    GSkyMap  m_dri;      //!< Data cube
    GEbounds m_ebounds;  //!< Energy boundaries of data cube
    GGti     m_gti;      //!< Good Time Intervals of data cube
    double   m_phimin;   //!< Phibar minimum (deg)
    double   m_phibin;   //!< Phibar binsize (deg)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMDri").
 ***************************************************************************/
inline
std::string GCOMDri::classname(void) const
{
    return ("GCOMDri");
}


/***********************************************************************//**
 * @brief DRI bin access operators
 *
 * @param[in] index DRI bin index [0,...,size()-1].
 * @return Reference to DRI bin.
 ***************************************************************************/
inline
double& GCOMDri::operator[](const int& index)
{
    return (const_cast<double&>((m_dri.pixels()[index])));
}


/***********************************************************************//**
 * @brief DRI bin access operators (const version)
 *
 * @param[in] index DRI bin index [0,...,size()-1].
 * @return Reference to DRI bin.
 ***************************************************************************/
inline
const double& GCOMDri::operator[](const int& index) const
{
    return (m_dri.pixels()[index]);
}


/***********************************************************************//**
 * @brief Return number of bins
 *
 * @return Number of bins.
 ***************************************************************************/
inline
int GCOMDri::size(void) const
{
    return (m_dri.npix()*m_dri.nmaps());
}


/***********************************************************************//**
 * @brief Return number of Chi bins
 *
 * @return Number of Chi bins.
 ***************************************************************************/
inline
int GCOMDri::nchi(void) const
{
    return (m_dri.nx());
}


/***********************************************************************//**
 * @brief Return number of Psi bins
 *
 * @return Number of Psi bins.
 ***************************************************************************/
inline
int GCOMDri::npsi(void) const
{
    return (m_dri.ny());
}


/***********************************************************************//**
 * @brief Return number of Phibar bins
 *
 * @return Number of Phibar bins.
 ***************************************************************************/
inline
int GCOMDri::nphibar(void) const
{
    return (m_dri.nmaps());
}


/***********************************************************************//**
 * @brief Return DRI sky map
 *
 * @return Sky map containing DRI data.
 ***************************************************************************/
inline
const GSkyMap& GCOMDri::map(void) const
{
    return (m_dri);
}


/***********************************************************************//**
 * @brief Return energy boundaries of DRI data
 *
 * @return Energy boundaries of DRI data.
 ***************************************************************************/
inline
const GEbounds& GCOMDri::ebounds(void) const
{
    return (m_ebounds);
}


/***********************************************************************//**
 * @brief Return Good Time Intervals of DRI data
 *
 * @return Good Time Intervals of DRI data.
 ***************************************************************************/
inline
const GGti& GCOMDri::gti(void) const
{
    return (m_gti);
}


/***********************************************************************//**
 * @brief Return minimum Compton scatter angle of DRI data
 *
 * @return Minimum Compton scatter angle of DRI data.
 ***************************************************************************/
inline
const double& GCOMDri::phimin(void) const
{
    return (m_phimin);
}


/***********************************************************************//**
 * @brief Return Compton scatter angle bin of DRI data
 *
 * @return Compton scatter angle bin of DRI data.
 ***************************************************************************/
inline
const double& GCOMDri::phibin(void) const
{
    return (m_phibin);
}

#endif /* GCOMDRI_HPP */
