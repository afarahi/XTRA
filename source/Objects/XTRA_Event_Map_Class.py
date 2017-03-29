import os
import json
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from scipy.ndimage.filters import gaussian_filter
from Objects import FDIR_SB_MAP_SAVE, FDIR_EVENT_MAP


def dist_circle(n, xcen=0 ,ycen=0):

    if type(n)==int or type(n)==float or type(n)==np.float64 or type(n)==np.int64:
        # Rectangle
        ny = n
        nx = n
    else:
        # Square
        nx = n[0]
        ny = n[1] 

    if xcen==0 and ycen==0:
        # Putting it in the middle
        xcen = (nx)/2.
        ycen = (ny)/2.
    
    x_2 = (np.arange(0, nx, 1.0) - xcen)**2     # X distances (squared)
    y_2 = (np.arange(0, ny, 1.0) - ycen)**2     # Y distances (squared)
    im = np.zeros((nx, ny))      # Make uninitialized output array
    for i in range(nx):              # Row loop
        for j in range(ny):
            im[i,j] = np.sqrt(x_2[i] + y_2[j])     # Euclidian distance
    return im


def xy2polar(x_pnt, y_pnt, x, y):

    pixelsize = 1.20833333333333E-03 # in degrees
    pixelsize = pixelsize * 60.0 #in arcmin

    # find theta
    dist_in_pixels = np.sqrt((y_pnt-y)**2 + (x_pnt-x)**2)
    theta = dist_in_pixels * pixelsize # pixelsize in arcmin 

    return theta


def ensure_dir(f):
    if not os.path.exists(f):
        os.makedirs(f)


def save_img_fit(img, fname, header):

    # remove existing file
    if os.path.exists(fname):  os.remove(fname)

    # save new image
    hdu = pyfits.PrimaryHDU(img, header)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(fname)


class Event_Map_Class:

    def __init__(self):

        with open('parameters/Event_Map_Parameters.json') as fp:
            _param = json.load(fp)

        # Grid Size
        self.grid_size = _param["grid_size"]
        # Exposure Time
        self.ks = _param["exposure_time"]
        del _param

    def saveEventMapPic(self, img, fname, smoothed='False'):

        sizes = np.shape(img)
        fig = plt.figure(figsize=(1, 1))
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        if smoothed:
            ax.imshow(np.log10(img + 1.0).T, aspect='auto', origin='lower',
                      vmin=0, vmax=1, cmap='bone', interpolation='gaussian')
            plt.savefig(fname[:-4] + '_smoothed.png', dpi=sizes[0])
        else:
            ax.imshow(np.log10(img+1.0).T, aspect='auto', origin='lower',
                       vmin=0, vmax=1, cmap='bone', interpolation='nearest')
            plt.savefig(fname, dpi=sizes[0])
        plt.close()

    def make_event_map(self, i, j, ks):

        self.obsid = j + i * self.grid_size
        header = self.header

        # 2) Cutout overlapping regions.
        centre_x = j * 256
        centre_y = 256 + i * 256

        image = self.whole_map[centre_y - 256:centre_y + 256, centre_x - 256:centre_x + 256]

        string_number = self.fname[:-4]
        ensure_dir(FDIR_EVENT_MAP + string_number)
        # print FDIR_EVENT_MAP + string_number
        string_number += '/%010d' % self.obsid
        ensure_dir(FDIR_EVENT_MAP + string_number)
        # print FDIR_EVENT_MAP + string_number
        string_number2 = '%010d' % self.obsid

        save_img_fit(image, FDIR_EVENT_MAP + string_number + '/origImage.fits', header)

        # 3) Create header
        ra_pnt = (j * 256) * 0.017917
        dec_pnt = (256 + i * 256) * 0.001194

        header['NAXIS1'] = 512
        header['NAXIS2'] = 512

        header['CRPIX1'] = 256.5
        header['CRPIX2'] = 256.5
        header['CRVAL1'] = ra_pnt
        header['CRVAL2'] = dec_pnt
        header['CTYPE1'] = 'RA---TAN'
        header['CTYPE2'] = 'DEC--TAN'
        header['CUNIT1'] = 'deg'
        header['CUNIT2'] = 'deg'
        header['CROTA2'] = 0.00000000000000e0
        header['CDELT1'] = 1.20833333333333e-3
        header['CDELT2'] = 1.20833333333333e-3
        header['RA'] = 1.20833333333333e-3
        header['DEC'] = 1.20833333333333e-3
        header['PA_PNT'] = 0
        header['EMSCI001'] = 'EMOS1'
        header['EXP_TIME'] = ks

        # 4) Blur with psf.
        # - This actually blurs the psf into an 11x11 square
        #      but as it drops off radially quite quickly this should be fine

        a = 5.074
        b = -0.236
        c = 0.002
        d = -0.0180
        x = 1.472
        y = -0.010
        z = -0.001
        w = -0.0016
        E = 1.5

        image_psf = 0.0 * gaussian_filter(image, sigma=0.005)
        arr = dist_circle(11, 5, 5)
        arr = arr * 0.5  # ; distarr in arcsec

        for k in range(0, 512):
            for l in range(0, 512):

                if ((k - 255) ** 2 + long(l - 255) ** 2 < long(239) ** 2):
                    theta_in = xy2polar(255, 255, k, l)

                    r_c = a + b * E + c * theta_in + d * E * theta_in
                    alpha = x + y * E + z * theta_in + w * E * theta_in

                    psf = np.power((1. + (arr / r_c) ** 2), (-1. * alpha))

                    image_psf[k - 5:k + 6, l - 5:l + 6] += (psf * image[k, l] / np.sum(psf))

        save_img_fit(image_psf, FDIR_EVENT_MAP + string_number + '/imagepsf.fits', header)

        # 4) Vignette,
        #    Multiply by exposure time,
        #    Draw value per pixel from poisson distribution.
        #    Add flat noise across field,
        #    Add vignetted noise

        exposure_map = pyfits.open('./parameters/expmap.fits')[0].data

        image_psf = image_psf * exposure_map * 2.4966 * 1e14 * ks
        image_psf = np.random.poisson(lam=image_psf)

        save_img_fit(image_psf, FDIR_EVENT_MAP + string_number + '/imageprenoise.fits', header)

        # 5) Add noise
        blanker = np.zeros([512, 512])
        blanker[exposure_map > 0] = 1

        image_psf = image_psf + \
                    np.random.poisson(lam=0.0015 * ks * exposure_map) + \
                    np.random.poisson(lam=0.0008 * ks, size=(512, 512)) * blanker

        # 6) Save file
        ensure_dir(FDIR_EVENT_MAP + string_number + '/images/')
        save_img_fit(image_psf,
                     FDIR_EVENT_MAP + string_number + '/images/' + string_number2 +
                     '-0.50-2.00keVmerged_img.fits', header)

        save_img_fit(exposure_map,
                     FDIR_EVENT_MAP + string_number + '/images/' + string_number2 +
                     '-0.50-2.00keVmerged_expmap.fits', header)

        save_img_fit(image_psf,
                     FDIR_EVENT_MAP + string_number + '/images/' + string_number2 +
                     '-0.50-2.00keV-mos1_merged_img.fits', header)

        save_img_fit(exposure_map,
                     FDIR_EVENT_MAP + string_number + '/images/' + string_number2 +
                     '-0.50-2.00keV-mos1_merged_expmap.fits', header)

        # 7) save image
        sub_dir = self.fname[:-4] + '/img'
        ensure_dir(FDIR_EVENT_MAP + sub_dir)
        self.saveEventMapPic(image_psf, FDIR_EVENT_MAP + sub_dir + '/%010d_event_img.png'%self.obsid)

    def run_tiled_map(self, fname):

        with open('parameters/Map_Parameters.json') as fp:
            _param = json.load(fp)

        XrayBandstr = _param["Flux_energy_band"]

        self.fname = fname[:-4] + '_' + XrayBandstr + '.fit'

        # 1) Read in initial fits file.
        # whole_map=pyfits.open(fdir+'maps/'+'%i'%realization_id+'/'+fname+'.fit')[0].data
        self.whole_map = pyfits.open(FDIR_SB_MAP_SAVE + self.fname)[0].data
        self.header = pyfits.open(FDIR_SB_MAP_SAVE + self.fname)[0].header

        for i in range(0, self.grid_size):
            for j in range(1, self.grid_size+1):
                print i, j, j + i * self.grid_size
                self.make_event_map(i, j, self.ks)



