/***************************************************************************
 *              GLATPsf.hpp - Fermi LAT point spread function              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
 * @file GLATPsf.hpp
 * @brief Fermi LAT point spread function class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATPSF_HPP
#define GLATPSF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GLATPsfBase.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GFits;


/***********************************************************************//**
 * @class GLATPsf
 *
 * @brief Interface for the Fermi LAT point spread function
 *
 * This class provides the interface to the Fermi LAT point spread function
 * (PSF). An instance of the class holds the PSF for a specific detector
 * section (front or back). The class handles different PSF versions by
 * holding a pointer to the versioned PSF, which is defined by the abstract
 * base class GLATPsfBase. On loading (or reading) of the PSF information
 * from the FITS file, the appropriate versioned PSF will be allocated.
 *
 * @todo Implement Phi dependence
 ***************************************************************************/
class GLATPsf : public GBase {

public:
    // Constructors and destructors
    GLATPsf(void);
    GLATPsf(const GFilename& filename, const std::string& evtype);
    GLATPsf(const GLATPsf& psf);
    virtual ~GLATPsf(void);

    // Operators
    GLATPsf& operator= (const GLATPsf& psf);
    double   operator()(const double& offset, const double& logE,
                        const double& ctheta);

    // Methods
    void        clear(void);
    GLATPsf*    clone(void) const;
    std::string classname(void) const;
    void        load(const GFilename& filename, const std::string& evtype);
    void        save(const GFilename& filename,
                     const bool&      clobber = false);
    void        read(const GFits& file);
    void        write(GFits& file) const;
    int         size(void) const;
    int         nenergies(void) const;
    int         ncostheta(void) const;
    double      costhetamin(void) const;
    void        costhetamin(const double& ctheta);
    bool        has_phi(void) const;
    bool        is_front(void) const;
    bool        is_back(void) const;
    int         version(void) const;
    std::string print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATPsf& psf);
    void free_members(void);
    
    // Members
    std::string  m_evtype; //!< Event type
    GLATPsfBase* m_psf;    //!< Pointer to versioned point spread function
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATPsf").
 ***************************************************************************/
inline
std::string GLATPsf::classname(void) const
{
    return ("GLATPsf");
}

#endif /* GLATPSF_HPP */
