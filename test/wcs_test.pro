;==============================================================================
; Script: wcs_test
; ----------------------------------------------------------------------------
; Create WCS test map
; ----------------------------------------------------------------------------
; Parameters:
;  ra       - Right Ascension of map centre (deg)
;  dec      - Declination of map centre (deg)
;  pixsize  - Pixelsize for FITS map
;  proj     - Map projection
;  filename - Filename
;==============================================================================
pro wcs_test, ra=ra, dec=dec, pixsize=pixsize, proj=proj, filename=filename


;====================
; Default parameters
;=====================
if n_elements(ra)       eq 0 then ra       = 83.63
if n_elements(dec)      eq 0 then dec      = 22.0
if n_elements(pixsize)  eq 0 then pixsize  = 0.02
if n_elements(size)     eq 0 then size     = 200
if n_elements(proj)     eq 0 then proj     = 'CAR'
if n_elements(filename) eq 0 then filename = 'map.fits'


;==============================
; Build astronometry structure
;==============================
galactic = 0
equinox  = 2000.0
cdelt    = [pixsize, pixsize]
crpix    = [double(size)/2.0+0.5, double(size)/2.0+0.5]
crval    = [ra, dec]
ctype    = ['RA---'+proj,'DEC--'+proj]
make_astr, astr, delt=cdelt, crpix=crpix, crval=crval, ctype=ctype


;=========
; Set map
;=========
map = dblarr(size, size)
for ix=0,size-1 do begin
  x = replicate(double(ix),size)
  y = dindgen(size)
  xy2ad, x, y, astr, a, d
  gcirc, 0, 83.63*!dtor, 22.0*!dtor, a*!dtor, d*!dtor, distance
  distance = distance/!dtor
  inx = where(distance le 0.1, nvalues)
  if (nvalues gt 0) then map[ix,inx] = 1.0
endfor


;==========
; Save map
;==========
mkhdr, hdr, map
putast, hdr, astr ;, cd_type=1
;sxdelpar, hdr, 'LONPOLE'
;sxdelpar, hdr, 'LATPOLE'
;sxdelpar, hdr, 'PV2_1'
mwrfits, map, filename, hdr, /create

end
