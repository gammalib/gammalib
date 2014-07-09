from gammalib import *
from ctools import *
from userfunc import *
import math

obs = GObservations("/lfs/l2/hess/users/chiachun/CubeAnalysis/cta1dc/crab/Avg/fit_binned_crab_1run.xml")
cobs = GCTAObservation(obs[0])
rsp = cobs.response();
pnt = cobs.pointing().dir();
x = 83.63
y = 22.01
dx = 0.1
dy = 0.1
nx = 40
ny =40
nebins = 15
emin = 0.1
emax = 100
min = 0
max = 0.2
nbins = 50

psfcube = GCTAMeanPsf(obs, x, y, dx, dy, nx, ny, emin, emax, nebins, min, max, nbins)
psfmap = psfcube.map()
ctrmap = GSkymap('CAR','CEL', x, y, dx, dy, nx, ny, nebins)


# test psf
for j in range(0, psfcube.ebounds().size()):
    for i in range(0, psfmap.npix()):
        deltasqs = []
        psfvals = []
        for k in range(0, psfcube.deltas().size()):
            idx = k + j*psfcube.deltas().size()
            deltasqs.append( psfcube.deltas()[k]*psfcube.deltas()[k] )
            psfvals.append( psfmap[i,idx] )
        ctr = searchctr(deltasqs, psfvals)
        ctrmap[i,j] = ctr
ctrmap.save('ctr_psfmap.fits',True) 

for ie in range(0, psfcube.ebounds().size()):
    for pix in range(0, psfmap.npix()):
        deltasqs = []
        psfvals = []
        for k in range(0, psfcube.deltas().size()):
            dir = psfmap.inx2dir(pix);
            theta = pnt.dist(dir); # radian
            idx = k + ie*psfcube.deltas().size()
            deltasqs.append( psfcube.deltas()[k]*psfcube.deltas()[k] )
            psfvals.append(rsp.psf(psfcube.deltas()[k]*gammalib.deg2rad, theta, 
                                   0.0, 0.0, 0.0, 
                                   psfcube.ebounds().emean(ie).log10TeV()) )
        ctr = searchctr(deltasqs, psfvals)
        ctrmap[pix,ie] = ctr
ctrmap.save('ctr_psf2d.fits',True) 

expcube = GCTAExposure(obs, x, y, dx, dy, nx, ny, emin, emax, nebins)
expmap = expcube.map()
arfmap = GSkymap('CAR','CEL', x, y, dx, dy, nx, ny, nebins)
# test arf
for ie in range(0, expcube.ebounds().size()):
    for pix in range(0, expmap.npix()):
        dir = expmap.inx2dir(pix);
        theta = pnt.dist(dir); # radian
        aeff = rsp.aeff(theta, 
                        0.0, 0.0, 0.0, 
                        expcube.ebounds().emean(ie).log10TeV())
        arfmap[pix,ie] = aeff
arfmap.save('arfmap.fits',True) 
