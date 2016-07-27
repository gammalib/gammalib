#! /usr/bin/env python
# ==========================================================================
# This script tests the HealPix boundaries method.
#
# --------------------------------------------------------------------------
#
# Copyright (C) 2015 Juergen Knoedlseder
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
import gammalib
import matplotlib.pyplot as plt

map   = gammalib.GSkyMap("GAL", 1, "RING")
print(map)

# Create figures
plt.figure(1)

# Loop over indices
for index in range(4,5):

    centre = map.inx2dir(index)
    proj   = map.projection()
    dirs   = proj.boundaries(gammalib.GSkyPixel(index), 3)
    print len(dirs)
    print dirs[0], dirs[1], dirs[2], dirs[3], dirs[4], dirs[5], dirs[6], dirs[7]

    l0 = centre.l_deg()
    if l0 > 180:
        l0 -= 360.0
    plt.plot([l0], [centre.b_deg()], "bo")
    plt.plot([dirs[0].l_deg(), dirs[1].l_deg(), \
              dirs[2].l_deg(), dirs[3].l_deg(), dirs[0].l_deg()], \
             [dirs[0].b_deg(), dirs[1].b_deg(), \
              dirs[2].b_deg(), dirs[3].b_deg(), dirs[0].b_deg()], "ro-")
    plt.plot([dirs[4].l_deg(), dirs[5].l_deg(), \
              dirs[6].l_deg(), dirs[7].l_deg(), dirs[4].l_deg()], \
             [dirs[4].b_deg(), dirs[5].b_deg(), \
              dirs[6].b_deg(), dirs[7].b_deg(), dirs[4].b_deg()], "go-")
    plt.plot([dirs[8].l_deg(), dirs[9].l_deg(), \
              dirs[10].l_deg(), dirs[11].l_deg(), dirs[8].l_deg()], \
             [dirs[8].b_deg(), dirs[9].b_deg(), \
              dirs[10].b_deg(), dirs[11].b_deg(), dirs[8].b_deg()], "bo-")

plt.xlim([ 180,-180])
plt.ylim([-90,90])

dir = gammalib.GSkyDir()
dir.lb_deg(1.0,1.0)
bilinear = proj.interpolator(dir)
print bilinear.index1(), bilinear.index2(), \
      bilinear.index3(), bilinear.index4()
print bilinear.weight1(), bilinear.weight2(), \
      bilinear.weight3(), bilinear.weight4()
point1 = map.inx2dir(bilinear.index1())
point2 = map.inx2dir(bilinear.index2())
point3 = map.inx2dir(bilinear.index3())
point4 = map.inx2dir(bilinear.index4())
l1 = point1.l_deg()
if l1 > 180:
    l1 -= 360.0
l2 = point2.l_deg()
if l2 > 180:
    l2 -= 360.0
l3 = point3.l_deg()
if l3 > 180:
    l3 -= 360.0
l4 = point4.l_deg()
if l4 > 180:
    l4 -= 360.0
plt.plot([l1], [point1.b_deg()], "go")
plt.plot([l2], [point2.b_deg()], "go")
plt.plot([l3], [point3.b_deg()], "go")
plt.plot([l4], [point4.b_deg()], "go")

plt.show()
