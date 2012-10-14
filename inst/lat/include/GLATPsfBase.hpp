/***************************************************************************
 * GLATPsfBase.hpp  -  Fermi/LAT point spread function abstract base class *
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
 * @file GLATPsfBase.hpp
 * @brief Fermi/LAT point spread function abstract base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATPSFBASE_HPP
#define GLATPSFBASE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GLATResponseTable.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GLATPsfBase
 *
 * @brief Abstract interface for the Fermi/LAT point spread function
 *
 * This class defines the abstract interface for the Fermi/LAT point spread
 * function (PSF). Classes that are derived from GLATPsfBase implement the
 * version dependent LAT PSFs.
 ***************************************************************************/
class GLATPsfBase : public GBase {

public:
    // Constructors and destructors
    GLATPsfBase(void);
    GLATPsfBase(const GLATPsfBase& psf);
    virtual ~GLATPsfBase(void);

    // Operators
    GLATPsfBase& operator= (const GLATPsfBase& psf);

    // Pure virtual methods
    virtual void         clear(void) = 0;
    virtual GLATPsfBase* clone(void) const = 0;
    virtual void         read(const GFitsTable* hdu) = 0;
    virtual void         write(GFits& file) const = 0;
    virtual double       psf(const double& offset, const double& logE,
                             const double& ctheta) = 0;
    virtual int          version(void) const = 0;
    virtual std::string  print(void) const = 0;

    // Other methods
    void   read_scale(const GFitsTable* hdu);
    void   write_scale(GFits& file) const;
    int    size(void) const { return nenergies()*ncostheta(); }
    int    nenergies(void) const { return m_rpsf_bins.nenergies(); }
    int    ncostheta(void) const { return m_rpsf_bins.ncostheta(); }
    double costhetamin(void) const { return m_min_ctheta; }
    void   costhetamin(const double& ctheta) { m_min_ctheta = ctheta; }
    bool   hasphi(void) const { return false; }
    bool   front(void) const { return m_front; }
    void   front(const bool& front) { m_front = front; }
    double scale_par1(void) const { return m_scale_par1; }
    double scale_par2(void) const { return m_scale_par2; }
    double scale_index(void) const { return m_scale_index; }

protected:
    // Methods
    void   init_members(void);
    void   copy_members(const GLATPsfBase& psf);
    void   free_members(void);
    double scale_factor(const double& energy) const;

    // Protected members
    bool              m_front;        //!< PSF is for front section?
    GLATResponseTable m_rpsf_bins;    //!< PSF energy and cos theta binning
    double            m_scale_par1;   //!< PSF scaling parameter 1
    double            m_scale_par2;   //!< PSF scaling parameter 2
    double            m_scale_index;  //!< PSF scaling index
    double            m_min_ctheta;   //!< Minimum valid cos(theta)
};

#endif /* GLATPSFBASE_HPP */
