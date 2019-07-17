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

    def saveMapPic(self,draw_halos=False,Halos=False):
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        # imports 
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        # check path; create if not exant
        path = FDIR_SB_MAP_SAVE
        if not os.path.exists(path): os.makedirs(path)
        fname = path + self.fname + r'_' + self.MPa.XrayBandstr + r'.pdf'
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        # set pixilation
        my_dpi = 96. # for WKB screen
        pxs = 750 # pixilation
        plt.figure(figsize=(pxs/my_dpi,pxs/my_dpi), dpi=my_dpi)
        # plt.clf()
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        # draw flux picture
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        
        # using vmax=None increases color contrast
        # but limits comparison between plots (formerly: vmax=-14)
        plt.imshow(np.log10(self.MAP / 6.24e8 + 1e-20).T,
                   aspect='auto', origin='lower', vmin=-18, vmax=-15.25,
                   extent=(self.RAmin, self.RAmax, self.DECmin, self.DECmax),
                   cmap='afmhot') # 'gist_heat') #'afmhot')
        
        # set tight colorbar for flux
        ax1 = plt.gca()
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(label=r'$\log_{10}$ flux (erg/s/cm$^2$/px)',cax=cax)
        plt.sca(ax1) # return to original axis
        
        # set axis labels 
        plt.xlabel('RA', {'fontsize': 20})
        plt.ylabel('DEC', {'fontsize': 20})
        
        if not draw_halos: # add a grid onto the data
          plt.grid()
          plt.title('DEC =%.2f , RA =%.2f ' % (self.DECc, self.RAc), {'fontsize': 20})
        else: # WKB
          ax1.axhline(0,linestyle='--',color='k')
          ax1.axvline(0,linestyle='--',color='k')
          print "Drawing R500 for each halo"
          
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          # import the halo set
          
          # read in xc_bar:
          with open('parameters/Input_Parameters.json') as fp:
            _inputParam = json.load(fp) # read in all input parameters
          xc_bar = _inputParam["xc_bar"]
          
          # pull in halo attributes: 
          radii = Halos.Rc/xc_bar
          if 0: print 'radii:',radii[-10:],max(radii)
          ra,dec = np.array(Halos.RA),np.array(Halos.DEC)
          # there's no "Halos.M500" attribute, so color by redshift
          # nor {M200, M200b, M2500, MVIR, mass, M}!
          z = Halos.Z
          
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          # for each halo in the set
          
          # create colormapping
          from matplotlib import cm
          colormap = 'coolwarm'
          cw = cm.get_cmap(colormap)
          if 0: # use dynamic mapping (different for each plot)
            colors = cw((z-min(z))/(max(z)-min(z)))
          else: # use static mapping (0.1--0.3 for each plot) 
            min_z,max_z = .1,.3
            colors = cw((z-min_z)/(max_z-min_z))
          
          # create list of circles to draw
          output_list = []
          for i in range(len(ra)):
              output_list.append(plt.Circle((ra[i],dec[i]),radii[i],
                                 color=colors[i],fill=False))
          ax = plt.gca() #aspect='equal') # for some reason, equal messes things up...
          L = self.MPa.DECMapSize/2
          # ax.set_xlim(-L,L); ax.set_ylim(-L,L)
          for circle in output_list:    
            ax.add_artist(circle)
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        # draw clusters using the argv for HID1
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        from sys import argv
        if len(argv)>3: # draw clusters using the argv for an HID
          HID = int(argv[3])
          print "Drawing cluster radius using halo %i" % HID
          
          if 0: HID = 10450059 # 13133502 # HARDCODE
          
          # scan in from xtools bicycle output
          xtools = '/home/wkblack/projects/xproject/xtools/'
          halo_fname = xtools + 'output_bicycle.csv'
          dat = np.loadtxt(halo_fname,delimiter=',',skiprows=1)
          HID1 = dat[:,0]; idx_mask = HID1==int(HID)
          cluster_lam = dat[idx_mask,4]; cluster_z = dat[idx_mask,5]
          cluster_ra = dat[idx_mask,8]; cluster_dec = dat[idx_mask,9]
          
          # check if halo is duplicate; if so, read option from argv
          if len(cluster_lam)>1: # select
            print "Halo is degenerate; reading in selection from argv."
            i = int(argv[4]) # choice
            cluster_lam = cluster_lam[i]; cluster_z = cluster_z[i]
            cluster_ra = cluster_ra[i]; cluster_dec = cluster_dec[i]
          
          if 1: # grab nearby clusters
            thetamax = .2 # maximum R_lambda
            width = self.MPa.DECMapSize/2 + thetamax
            min_ra = cluster_ra-width
            max_ra = cluster_ra+width
            min_dec = cluster_ra-width
            max_dec = cluster_ra+width
            proximity_mask = (dat[:,6]>min_ra)*(dat[:,6]<max_ra) \
                            *(dat[:,7]>min_dec)*(dat[:,7]<max_dec)
            print 'proximity count :',len(proximity_mask)
          
          # calculate thetas using cosmology of Aardvark
          from sys import path
          path.insert(0,xtools)
          from cosmology import theta_R_lambda, to_degrees
          
          theta_R_lambda = theta_R_lambda(cluster_lam,cluster_z) # 0.086697
          theta_500kpc = to_degrees(.5,cluster_z) # 500 kpc in degrees
          
          # draw on circles (background & foreground)
          ax.add_artist(plt.Circle((0,0),theta_R_lambda,
                         color='k',fill=False,lw=3))
          ax.add_artist(plt.Circle((0,0),theta_R_lambda,
                         color='lime',fill=False))
          
          ax.add_artist(plt.Circle((0,0),theta_500kpc,
                         color='k',fill=False,lw=3))
          ax.add_artist(plt.Circle((0,0),theta_500kpc,
                         color='w',fill=False))
          plt.title(r'Cluster HID$_1$ : %i' '\n'
                    r'(RA,DEC,$z$,$\lambda$)=(%.2f,%.2f,%.2f,%.2f)' \
                    % (HID,cluster_ra,cluster_dec,cluster_z,cluster_lam),{'fontsize': 20})
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        if 0: # set tight colorbar for redshift
          # divider = make_axes_locatable(ax1)
          # cax = divider.append_axes("top", size="5%", pad=0.05)
          plt.scatter(colors,colors,c=z3,cmap='coolwarm') 
          plt.gca().set_visible(False)
          plt.colorbar(label=r'Redshift $z$',orientation='horizontal')
          plt.sca(ax1) # return to original axis
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
        # finish up and exit 
        # plt.tight_layout()
        plt.savefig(fname,dpi=my_dpi)
        plt.close()

    def saveMapFits(self):

        import astropy.io.fits as pyfits
        from Objects import FDIR_CLUSTERS

        path = FDIR_SB_MAP_SAVE
        if not os.path.exists(path): os.makedirs(path)
        fdir = path + self.fname + r'_' + self.MPa.XrayBandstr + r'.fit'

        print("Removing old halo/map catalog (if exists) : %s"%fdir)
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


    def drawHaloBounds(self, Halos): 
        print "Drawing..."
        # get input ra+dec and R500 scaled to an angular size (Rc/xc_bar)
        
        # read in xc_bar:
        with open('parameters/Input_Parameters.json') as fp:
          _inputParam = json.load(fp) # read in all input parameters
        xc_bar = _inputParam["xc_bar"]
        
        # pull in halo attributes: 
        radii = Halos.Rc/xc_bar
        if 1: print 'radii:',radii[-10:],max(radii)
        ra,dec = np.array(Halos.RA),np.array(Halos.DEC)
        z = Halos.Z
        try: 
          m500 = Halos.M500
        except AttributeError: 
          print "no mass!"
        
        from matplotlib import cm
        # nor {M200, M200b, M2500, MVIR, mass, M}!
        colormap = 'coolwarm'
        cw = cm.get_cmap(colormap)
        colors = cw((z-min(z))/(max(z)-min(z))) # (0.1,0.2,0.3)
        
        output_list = []
        import matplotlib.pyplot as plt
        for i in range(len(ra)):
            output_list.append(plt.Circle((ra[i],dec[i]),radii[i],
                               color=colors[i],fill=False))
        ax = plt.gca(aspect='equal')
        ax.cla()
        L = self.MPa.DECMapSize/2
        ax.set_xlim(-L,L); ax.set_ylim(-L,L)
        for circle in output_list:    
           ax.add_artist(circle)
        
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        outname = "halo_circles.png"
        plt.savefig(outname)
        # plt.show() 


