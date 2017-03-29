from Objects import FDIR_SB_MAP_SAVE
import numpy as np
import json
import os


class Map_Parameters_Class:

    def __init__(self):

        with open('parameters/Map_Parameters.json') as fp:
            _param = json.load(fp)

        # Energy Band
        self.XrayBandstr = _param["Flux_energy_band"]
        # Map Intensity Unit
        self.mapUnit = _param["Unit"]
        # Map size in degree (it can not be lager than 5 degree!!)
        self.DECMapSize = _param["DEC_map_size"]  # in degree
        self.RAMapSize = _param["RA_map_size"]  # in degree
        # Number of pixels in each side
        self.PixSide = _param["Num_pixels"]


class Map_Class:

    def __init__(self):

        self.MPa = Map_Parameters_Class()
        # Map is nxn array which firs element is the RA and second one is DEC
        self.MAP = np.zeros([self.MPa.PixSide, self.MPa.PixSide])
        # Map file name
        self.fname = 'sample.fits'
        from Surface_Brightness_Models import Surface_Brightness_Class
        SBC = Surface_Brightness_Class()
        self.flux2SurfaceBrightness = SBC.F2SB_Class.flux2SurfaceBrightness
        del SBC

        # Number of pixels = nside x nside
        # ( DECmax - DECmin ) / dDEC = nside
        # ( RAmax  - RAmin  ) / dRA  = nside

        # Center of the map
        self.DECc = 0.0
        self.RAc = 0.0

        self.DECmin = 0.0
        self.RAmin = 0.0
        self.DECmax = 0.0
        self.RAmax = 0.0

        self.dDEC = 0.0
        self.dRA = 0.0

    def update(self):

        self.i, self.j = np.mgrid[0:self.MPa.PixSide, 0:self.MPa.PixSide]
        self.DECmin = - self.MPa.DECMapSize / 2.
        self.RAmin = - self.MPa.RAMapSize / 2.
        self.DECmax = self.MPa.DECMapSize / 2.
        self.RAmax = self.MPa.RAMapSize / 2.
        self.nside = self.MPa.PixSide
        self.dDEC = self.MPa.DECMapSize / float(self.nside)
        self.dRA = self.MPa.RAMapSize / float(self.nside)
        self.DECc = 0.
        self.RAc = 0.
        self.fname = 'sample'

    def ang2pix(self, RA, DEC):
        # RA  component
        i = int((RA - self.RAmin) / self.dRA)
        # DEC component
        j = int((DEC - self.DECmin) / self.dDEC)
        return (i, j)

    def pix2ang(self, i, j):
        RA = self.RAmin + self.dRA / 2.0 + float(i) * self.dRA
        DEC = self.DECmin + self.dDEC / 2.0 + float(j) * self.dDEC
        return (RA, DEC)

    # Exact formula for angular distance
    # (Note that RA and DEC should be in radian unit)
    def angularDistance(self, RA1, DEC1, RA2, DEC2):
        return np.arccos(np.sin(DEC1) * np.sin(DEC2) + np.cos(DEC1) * np.cos(DEC2) * np.cos(RA1 - RA2))

    def pixInsideArea(self, RA, DEC, angularDis):
        pixi = []
        pixj = []
        angularDis2 = angularDis * angularDis
        (imin, jmin) = self.ang2pix(RA - angularDis, DEC - angularDis)
        (imax, jmax) = self.ang2pix(RA + angularDis, DEC + angularDis)
        for i in range(imin - 1, imax + 1):
            for j in range(jmin - 1, jmax + 1):
                (newRA, newDEC) = self.pix2ang(i, j)
                if (((RA - newRA) ** 2 + (DEC - newDEC) ** 2) < angularDis2):
                    pixi.append(i)
                    pixj.append(j)
        return (pixi, pixj)

    # Working in real space (not log space)
    def fluxCalculator(self, RA, DEC, angularDis):
        (pixi, pixj) = self.pix_inside_halo(RA, DEC, angularDis)
        flux = 0.0
        for i in range(len(pixi)):
            if (pixi[i] < 0 or pixj[i] < 0
                or pixi[i] > self.nside - 1
                or pixj[i] > self.nside - 1):
                pass
            else:
                flux += self.Map[pixi[i], pixj[i]]
        return flux

    def showMap(self):
        import matplotlib.pyplot as plt
        plt.figure()
        plt.clf()
        plt.imshow(self.MAP.T, aspect='auto', origin='lower'
                   , extent=(self.RAmin, self.RAmax, self.DECmin, self.DECmax))
        plt.colorbar()
        plt.xlabel('RA', {'fontsize': 20})
        plt.ylabel('DEC', {'fontsize': 20})
        plt.title('DEC =%.2f , RA =%.2f ' % (self.DECc, self.RAc), {'fontsize': 20})
        plt.grid()
        plt.show()

    def saveMapPic(self):
        import matplotlib.pyplot as plt

        path = FDIR_SB_MAP_SAVE
        if not os.path.exists(path): os.makedirs(path)
        fname = path + self.fname + r'_' + self.MPa.XrayBandstr + r'.pdf'

        plt.figure()
        plt.clf()
        plt.imshow(np.log10(self.MAP / 6.24e8 + 1e-20).T,
                   aspect='auto', origin='lower', vmin=-18, vmax=-14,
                   extent=(self.RAmin, self.RAmax, self.DECmin, self.DECmax))
        plt.colorbar()
        plt.xlabel('RA', {'fontsize': 20})
        plt.ylabel('DEC', {'fontsize': 20})
        plt.title('DEC =%.2f , RA =%.2f ' % (self.DECc, self.RAc), {'fontsize': 20})
        plt.grid()
        plt.savefig(fname)
        plt.close()

    def saveMapFits(self):

        import astropy.io.fits as pyfits
        from Objects import FDIR_CLUSTERS

        path = FDIR_SB_MAP_SAVE
        if not os.path.exists(path): os.makedirs(path)
        fdir = path + self.fname + r'_' + self.MPa.XrayBandstr + r'.fit'

        print("Removing old halo/map catalog (if exist) : %s"%fdir)
        os.system('rm -r -f %s' % fdir)
        print('Saving data %s'%self.fname)

        hduListObj = pyfits.open(FDIR_CLUSTERS + self.fname + '.fit')[1:]
        hdr = pyfits.open(FDIR_CLUSTERS + self.fname + '.fit')[0].header
        hdr.add_comment('------- MAP INFO -------- ')
        hdr.set('XRAY_MAP', 'TRUE')
        hdr.set('DEC_SIZE', str(self.MPa.DECMapSize), 'unit:degree')
        hdr.set('RA_SIZE', str(self.MPa.RAMapSize), 'unit:degree')
        hdr.set('PIX_NUM', str(self.MPa.PixSide))
        hdr.set('PIX_UNIT', self.MPa.mapUnit)
        hdr.add_comment(' --------------------------- ')
        hdr.add_comment('This map is created by XTRA version 1.0.0')

        hdu = pyfits.PrimaryHDU(data=self.MAP / 6.24e8 + 1e-20, header=hdr)
        hdu = [hdu]
        for i in range(len(hduListObj)):
            hdu.append(hduListObj[i])
        thdulist = pyfits.HDUList(hdu)
        thdulist.writeto(fdir)
        thdulist.close()

        print("Map file is generated.")

    def addClusters2Map(self, clusters):

        self.MPa.XrayBandstr = clusters.XrayBandstr

        for i in range(len(clusters.RA[::])):
            if clusters.F[i] == 0.:
                pass
            else:
                self.flux2SurfaceBrightness(self, clusters, i, self.MPa)

        self.fname = clusters.fname




