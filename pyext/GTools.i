/***************************************************************************
 *                        GTools.i - GammaLib tools                        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2018 by Juergen Knoedlseder                         *
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
 * @file GTools.i
 * @brief Gammalib tools definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTools.hpp"
%}

/* __ Constants __________________________________________________________ */
namespace gammalib {
    const double MeV2erg    =  1.6021765e-6;
    const double erg2MeV    =  624150.96;
    const double pc2cm      =  3.08568025e18;
    const double sec_in_day = 86400.0;
}

/* __ Prototypes ________________________________________________________ */
namespace gammalib {
    std::string              strip_whitespace(const std::string& arg);
    std::string              strip_chars(const std::string& arg,
                                         const std::string& chars);
    std::string              rstrip_chars(const std::string& arg,
                                          const std::string& chars);
    std::string              replace_segment(const std::string& arg,
                                             const std::string& segment,
                                             const std::string& replacement);
    std::string              expand_env(const std::string& arg);
    std::string              filepath(const std::string& pathname,
                                      const std::string& filename);
    std::string              toupper(const std::string& s);
    std::string              tolower(const std::string& s);
    std::vector<std::string> split(const std::string& s, const std::string& sep);
    std::string              fill(const std::string& s, int n);
    std::string              left(const std::string& s, int n, char c = ' ');
    std::string              right(const std::string& s, int n, char c = ' ');
    std::string              centre(const std::string& s, int n, char c = ' ');
    std::string              parformat(const std::string& s, const int& indent = 0);
    std::string              number(const std::string& noun, const int& number);
    double                   plaw_photon_flux(const double& emin,
                                              const double& emax,
                                              const double& epivot,
                                              const double& gamma);
    double                   plaw_energy_flux(const double& emin,
                                              const double& emax,
                                              const double& epivot,
                                              const double& gamma);
    GEnergy                  elogmean(const GEnergy& a, const GEnergy& b);
    bool                     dir_exists(const std::string& dirname);
    bool                     is_infinite(const double& x);
    bool                     is_notanumber(const double& x);
    bool                     contains(const std::string& str,
                                      const std::string& substring);
    bool                     contains(const std::vector<std::string>& strings,
                                      const std::string& string);
    void                     warning(const std::string& origin,
                                     const std::string& message);
    std::string              xml2str(const std::string& arg);
    std::string              str2xml(const std::string& arg);
    bool                     xml_has_par(const GXmlElement& xml,
                                         const std::string& name);
    GXmlElement*             xml_need_par(const std::string& origin,
                                          GXmlElement&       xml,
                                          const std::string& name);
    const GXmlElement*       xml_get_par(const std::string& origin,
                                         const GXmlElement& xml,
                                         const std::string& name);
    std::string              xml_get_attr(const std::string& origin,
                                          const GXmlElement& xml,
                                          const std::string& name,
                                          const std::string& attribute);
    void                     xml_check_par(const std::string& origin,
                                           const std::string& name,
                                           const int&         number);
    GFilename                xml_file_expand(const GXmlElement& xml,
                                             const std::string& filename);
    GFilename                xml_file_reduce(const GXmlElement& xml,
                                             const std::string& filename);
    int                      recv(int fd, char *buffer, int len, int flags,
                                  int timeout);
    double                   roi_arclength(const double& rad,
                                           const double& dist,
                                           const double& cosdist,
                                           const double& sindist,
                                           const double& roi,
                                           const double& cosroi);
}
