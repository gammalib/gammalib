/***************************************************************************
 *  GLATPsfBase.hpp - Abstract Fermi/LAT point spread function base class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @brief Abstract Fermi/LAT point spread function base class definition
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
 * @brief Abstract Fermi/LAT point spread function base class
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
    GLATPsfBase& operator=(const GLATPsfBase& psf);

    // Pure virtual methods
    virtual void         clear(void) = 0;
    virtual GLATPsfBase* clone(void) const = 0;
    virtual void         read(const GFitsTable& table) = 0;
    virtual void         write(GFits& file) const = 0;
    virtual double       psf(const double& offset, const double& logE,
                             const double& ctheta) = 0;
    virtual int          version(void) const = 0;
    virtual std::string  print(const GChatter& chatter = NORMAL) const = 0;

    // Other methods
    void          read_scale(const GFitsTable& hdu);
    void          write_scale(GFits& file) const;
    int           size(void) const;
    int           nenergies(void) const;
    int           ncostheta(void) const;
    const double& costhetamin(void) const;
    void          costhetamin(const double& ctheta);
    bool          has_phi(void) const;
    const bool&   front(void) const;
    void          front(const bool& front);
    const double& scale_par1(void) const;
    const double& scale_par2(void) const;
    const double& scale_index(void) const;

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


/***********************************************************************//**
 * @brief Return number of bins in point spread function
 *
 * @return Number of bins in point spread function.
 ***************************************************************************/
inline
int GLATPsfBase::size(void) const
{
    return nenergies()*ncostheta();
}


/***********************************************************************//**
 * @brief Return number of energies in point spread function
 *
 * @return Number of energies in point spread function.
 ***************************************************************************/
inline
int GLATPsfBase::nenergies(void) const
{
    return m_rpsf_bins.nenergies();
}


/***********************************************************************//**
 * @brief Return number of cosine theta bins in point spread function
 *
 * @return Number of cosine theta bins in point spread function.
 ***************************************************************************/
inline
int GLATPsfBase::ncostheta(void) const
{
    return m_rpsf_bins.ncostheta();
}


/***********************************************************************//**
 * @brief Return cosine theta minimum
 *
 * @return Cosine theta minimum.
 ***************************************************************************/
inline
const double& GLATPsfBase::costhetamin(void) const
{
    return m_min_ctheta;
}


/***********************************************************************//**
 * @brief Set minimum cos(theta) angle for point spread function
 *
 * @param[in] ctheta Cosine of maximum zenith angle.
 ***************************************************************************/
inline
void GLATPsfBase::costhetamin(const double& ctheta)
{
    m_min_ctheta = ctheta;
    return;
}


/***********************************************************************//**
 * @brief Signal that point spread function has Phi dependence
 *
 * @return True if point spread function has Phi dependence.
 ***************************************************************************/
inline
bool GLATPsfBase::has_phi(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal that point spread function is for front section
 *
 * @return True if point spread function is for front section.
 ***************************************************************************/
inline
const bool& GLATPsfBase::front(void) const
{
    return m_front;
}


/***********************************************************************//**
 * @brief Set if point spread function is for front section
 *
 * @param[in] front True if point spread function is for front section.
 ***************************************************************************/
inline
void GLATPsfBase::front(const bool& front)
{
    m_front = front;
    return;
}


/***********************************************************************//**
 * @brief Return first scaling parameter
 *
 * @return First scaling parameter.
 ***************************************************************************/
inline
const double& GLATPsfBase::scale_par1(void) const
{
    return m_scale_par1;
}


/***********************************************************************//**
 * @brief Return second scaling parameter
 *
 * @return Second scaling parameter.
 ***************************************************************************/
inline
const double& GLATPsfBase::scale_par2(void) const
{
    return m_scale_par2;
}


/***********************************************************************//**
 * @brief Return scaling index
 *
 * @return Scaling index.
 ***************************************************************************/
inline
const double& GLATPsfBase::scale_index(void) const
{
    return m_scale_index;
}

#endif /* GLATPSFBASE_HPP */
