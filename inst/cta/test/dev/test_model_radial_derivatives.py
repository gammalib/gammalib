#! /usr/bin/env python
# ==========================================================================
# Tests the formulae for the derivatives of a radial model.
#
# Copyright (C) 2020 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
import math
import gammalib


# ============ #
# Compute zeta #
# ============ #
def compute_zeta(dir_true, dir_reco, dir_model):
    """
    Compute zeta
    """
    # Compute zeta
    zeta = math.acos(math.sin(dir_model.dec()) * math.sin(dir_reco.dec()) +
                     math.cos(dir_model.dec()) * math.cos(dir_reco.dec()) *
                     math.cos(dir_model.ra() - dir_reco.ra()))

    # Compute gradients
    sin_zeta = math.sin(zeta)
    if sin_zeta != 0.0:
        dra  = (math.cos(dir_model.dec()) * math.cos(dir_reco.dec()) *
                math.sin(dir_model.ra() - dir_reco.ra())) / sin_zeta
        ddec = (math.sin(dir_model.dec()) * math.cos(dir_reco.dec()) *
                math.cos(dir_model.ra() - dir_reco.ra()) -
                math.cos(dir_model.dec()) * math.sin(dir_reco.dec())) / sin_zeta
    else:
        dra  = math.cos(dir_model.dec()) * math.cos(dir_reco.dec())
        ddec = 1.0
        print('dzeta/dra:  sin(zeta)=0', dra)
        print('dzeta/ddec: sin(zeta)=0', ddec)

    # Return zeta and gradients
    return zeta, dra, ddec


# =========== #
# Compute phi #
# =========== #
def compute_phi(dir_true, dir_reco, dir_model):
    """
    Compute phi
    """
    # Precompute values
    sin_beta_reco = math.sin(dir_reco.dec())
    cos_beta_reco = math.cos(dir_reco.dec())
    tan_beta_0    = math.tan(dir_model.dec())
    sin_dalpha    = math.sin(dir_model.ra() - dir_reco.ra())
    cos_dalpha    = math.cos(dir_model.ra() - dir_reco.ra())
    denom         = cos_beta_reco * tan_beta_0 - sin_beta_reco * cos_dalpha

    # Compute phi
    phi_true = math.atan2(math.sin(dir_true.ra() - dir_reco.ra()),
                          math.cos(dir_reco.dec()) * math.tan(dir_true.dec()) -
                          sin_beta_reco * math.cos(dir_true.ra() - dir_reco.ra()))
    phi_reco = math.atan2(math.sin(dir_model.ra() - dir_reco.ra()), denom)
    phi      = phi_true - phi_reco

    # Compute gradients
    denom2    = denom * denom
    dra_nom   = sin_beta_reco - cos_dalpha * cos_beta_reco * tan_beta_0
    ddec_nom  = sin_dalpha * cos_beta_reco * (1.0 + tan_beta_0*tan_beta_0)
    dra_denom = sin_dalpha * sin_dalpha + denom2
    if dra_denom != 0.0:
        dra  = dra_nom  / dra_denom
        ddec = ddec_nom / dra_denom
    else:
        print('dphi/dra:   denominator=0')
        print('dphi/ddec:  denominator=0')
        dra  = 0.0
        ddec = 0.0

    # Return phi and gradients
    return phi, dra, ddec


# ============= #
# Compute theta #
# ============= #
def compute_theta(dir_true, dir_reco, dir_model):
    """
    Compute theta
    """
    # Precompute values
    phi,  dphi_ra,  dphi_dec  = compute_phi(dir_true, dir_reco, dir_model)
    zeta, dzeta_ra, dzeta_dec = compute_zeta(dir_true, dir_reco, dir_model)

    # Compute theta
    theta = math.acos(math.sin(dir_true.dec()) * math.sin(dir_model.dec()) +
                      math.cos(dir_true.dec()) * math.cos(dir_model.dec()) *
                      math.cos(dir_true.ra() - dir_model.ra()))

    # Compute delta
    delta = math.acos(math.sin(dir_true.dec()) * math.sin(dir_reco.dec()) +
                      math.cos(dir_true.dec()) * math.cos(dir_reco.dec()) *
                      math.cos(dir_true.ra() - dir_reco.ra()))

    # Compute partial derivatives of theta with respect to zeta and phi
    sin_zeta  = math.sin(zeta)
    cos_zeta  = math.cos(zeta)
    sin_delta = math.sin(delta)
    cos_delta = math.cos(delta)
    sin_phi   = math.sin(phi)
    cos_phi   = math.cos(phi)

    #
    sin_theta = math.sin(theta)
    if sin_theta != 0.0:
        dzeta = (cos_delta * sin_zeta - sin_delta * cos_zeta * cos_phi) / sin_theta
        dphi  = (sin_delta * sin_zeta * sin_phi) / sin_theta
        print(dzeta, dphi, dzeta_ra, dphi_ra, dzeta * dzeta_ra  + dphi * dphi_ra)
        #print('beta ', phi, dphi,  dzeta * dzeta_dec,dphi * dphi_dec)
    else:
        print('dtheta/dzeta: sin(theta)=0')
        print('dtheta/dphi:  sin(theta)=0')
        dzeta = 0.0
        dphi  = 0.0
    
    # Compute partial derivaties of theta with respect to RA and Dec
    if zeta != 0.0:
        dra  = dzeta * dzeta_ra  + dphi * dphi_ra
        ddec = dzeta * dzeta_dec + dphi * dphi_dec
    else:
        dra  = 0.0 # ???? What value
        ddec = dzeta

    # Return theta and gradients
    return theta, dra, ddec


# ================ #
# Test derivatives #
# ================ #
def test_derivatives(function, dir_true, dir_reco, dir_model):
    """
    Test derivative
    """
    # Get modified model directions for numerical gradient
    dh         = 0.000001
    dra_plus   = gammalib.GSkyDir()
    dra_minus  = gammalib.GSkyDir()
    ddec_plus  = gammalib.GSkyDir()
    ddec_minus = gammalib.GSkyDir()
    dra_plus.radec(dir_model.ra()+dh,dir_model.dec())
    dra_minus.radec(dir_model.ra()-dh,dir_model.dec())
    ddec_plus.radec(dir_model.ra(),dir_model.dec()+dh)
    ddec_minus.radec(dir_model.ra(),dir_model.dec()-dh)

    # Test function
    f, dfda_ana, dfdb_ana = function(dir_true, dir_reco, dir_model)
    fa_plus, _, _         = function(dir_true, dir_reco, dra_plus)
    fa_minus, _, _        = function(dir_true, dir_reco, dra_minus)
    fb_plus, _, _         = function(dir_true, dir_reco, ddec_plus)
    fb_minus, _, _        = function(dir_true, dir_reco, ddec_minus)

    # Compute numerical gradients
    dfda_num = 0.5 * (fa_plus - fa_minus) / dh
    dfdb_num = 0.5 * (fb_plus - fb_minus) / dh

    # Collect result
    result = {'f': f,
              'dfda_ana': dfda_ana, 'dfdb_ana': dfdb_ana,
              'dfda_num': dfda_num, 'dfdb_num': dfdb_num}

    # Return
    return result


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':

    # Dump header
    print("")
    print("*********************************")
    print("* Test radial model derivatives *")
    print("*********************************")

    # Set test directions
    dir_true  = gammalib.GSkyDir()
    dir_reco  = gammalib.GSkyDir()
    dir_model = gammalib.GSkyDir()
    dir_true.radec_deg(85.0, +25.0)
    dir_reco.radec_deg(82.0, +23.0)
    dir_model.radec_deg(83.0, +22.0)
    
    # Test theta=0: dir_true = dir_model
    #dir_true.radec_deg(83.0, +22.0)
    
    # Test zeta=0: dir_reco = dir_model
    dir_reco.radec_deg(83.0, +22.0)
    #dir_reco.radec_deg(dir_reco.ra_deg(),dir_reco.dec_deg()-0.0001)
    
    #dir_true.radec_deg(82.0, +23.0)
    #dir_reco.radec_deg(82.0, +23.0)
    #dir_reco.radec_deg(82.001, +23.001)
    #dir_reco.radec_deg(85.0, +25.0)
    #dir_model.radec_deg(83.0, +22.0)
    #dir_model.radec_deg(85.0, +25.0)
    #dir_model.radec_deg(82.0, +23.0)

    # Test
    zeta  = test_derivatives(compute_zeta,  dir_true, dir_reco, dir_model)
    phi   = test_derivatives(compute_phi,   dir_true, dir_reco, dir_model)
    theta = test_derivatives(compute_theta, dir_true, dir_reco, dir_model)

    # Print results
    print('Zeta ...............: %f deg' % (zeta['f'] * gammalib.rad2deg))
    print('d(Zeta)/d(alpha_0) .: %9.6f (ana) %9.6f (num)' % (zeta['dfda_ana'], zeta['dfda_num']))
    print('d(Zeta)/d(beta_0) ..: %9.6f (ana) %9.6f (num)' % (zeta['dfdb_ana'], zeta['dfdb_num']))
    print('Phi ................: %f deg' % (phi['f'] * gammalib.rad2deg))
    print('d(Phi)/d(alpha_0) ..: %9.6f (ana) %9.6f (num)' % (phi['dfda_ana'], phi['dfda_num']))
    print('d(Phi)/d(beta_0) ...: %9.6f (ana) %9.6f (num)' % (phi['dfdb_ana'], phi['dfdb_num']))
    print('Theta ..............: %f deg' % (theta['f'] * gammalib.rad2deg))
    print('d(Theta)/d(alpha_0) : %9.6f (ana) %9.6f (num)' % (theta['dfda_ana'], theta['dfda_num']))
    print('d(Theta)/d(beta_0) .: %9.6f (ana) %9.6f (num)' % (theta['dfdb_ana'], theta['dfdb_num']))
