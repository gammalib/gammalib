# Make sky maps with pixel solid angle as value
# (for debugging)
import numpy as np
import gammalib

projection = 'AIT'

x, y = 0, 0
dx, dy = 1, 1
nx, ny = 360, 180
nmaps = 1
image = gammalib.GSkymap(projection, "CEL", x, y,
                         -dx, dy, nx, ny, nmaps)


#dir = gammalib.GSkyDir()
#dir.lb_deg(0, 0)
#pix = image.dir2pix(dir)

#solidangle = image.solidangle(pix)

#print(solidangle)



# Fill the sky map with the model image
for pix in range(image.npix()):
    try:
        STERADIAN_TO_DEG2 = np.degrees(1) ** 2
        image[pix] = STERADIAN_TO_DEG2 * image.solidangle(pix)
    except RuntimeError:
        pass

#import IPython; IPython.embed()

# Save the image to a FITS file
filename = 'test_{0}.fits'.format(projection)
print('Writing {0}'.format(filename))
image.save(filename, True)
