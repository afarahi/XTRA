import os
import json
import numpy as np
import astropy.io.fits as pyfits
from XTRA_pkg import __version__

def Read_Halos_Catalog_fit_file(fname, fdir):

    fdir = fdir + fname
    print "Reading Halos Catalog :", fdir

    with open('parameters/Input_Parameters.json') as fp:
        _param = json.load(fp)
    mlim = _param['Mass_limit']

    try:
        Halos_list = pyfits.open(fdir)
        Halos_head = Halos_list[1].header
        Halos_data = Halos_list[1].data[::]
        Halos_data = Halos_data[ Halos_data.M500 > mlim ]
        Halos_list.close()
    except IOError:
        print("Error: can\'t find file or read file")
        exit(1)

    return (fname[:-4], Halos_data, Halos_list[0].header)


class Halos_Class():

   def __init__(self):

       from Objects import Input_Parameters, Halo_Status_Class, Solvers_Class

       self.InputParam = Input_Parameters()
       self.Status = Halo_Status_Class()
       self.Solvers = Solvers_Class()
       # self.Solvers.ProperDistanceTabulate(self.InputParam, 4.0)

       self.ID   = []; self.pd   = []; self.RA   = []; self.DEC   = []
       self.R500 = []; self.M500 = []; self.Z_red = []
       self.lgT  = []; self.lgLx = []; self.lgFx = []
       self.Rc   = []; self.beta = []
       self.XCatRA  = [];  self.XCatDEC = []
       self.RAref  = 0.0
       self.DECref = 0.0
       # self.OutputParam = Output_Parameters()

   def readPDistance(self, Z): # read proper distance

       ztab, rtab = np.loadtxt("./Output/tabulated_data/Proper_Distance.txt", unpack=True)
       lenztab = len(ztab); maxztab = max(ztab)
       PDis = np.zeros(len(Z))
       for i in range(len(Z)):
           zloc= int(lenztab*Z[i]/maxztab)
           # interpolate:
           PDis[i] = rtab[zloc] + (Z[i]-ztab[zloc])*(rtab[zloc+1]-rtab[zloc])/(ztab[zloc+1]-rtab[zloc])
       return PDis

   def addHalosCatalog(self, fname=None):

       from numpy import zeros
       from Objects import FDIR_HALOS

       self.fname = fname

       fname, Halos, hdr = Read_Halos_Catalog_fit_file(fname=self.fname, fdir=FDIR_HALOS)
       
       if len(Halos)<1:
           print "ERROR: halo dataset empty! Perhaps mass limit is too high?"
           # vetting occurs here, at: XTRA/source/Objects/XTRA_Halos_Class.py +20
           raise SystemExit
       
       if fname:
           self.RA = Halos['RA']
           self.DEC = Halos['DEC']
           self.R500 = Halos['R500']
           self.M500 = Halos['M500']
           self.Z_red = Halos['Z']
           self.pd = self.readPDistance(Halos['Z'])

           self.ID = Halos['HALOID']
           self.XCatRA = Halos['RA']
           self.XCatDEC = Halos['DEC']

           n = len(self.RA[:])
           self.Rc = list(zeros(n))
           self.lgT = list(zeros(n))
           self.lgLx = list(zeros(n))
           self.lgFx = list(zeros(n))
           self.beta = list(zeros(n))
           self.Z_min = min(self.Z_red)
           self.Z_max = max(self.Z_red)
           self.fname = fname
           self.number_of_halos = n
           self.Status.update(self)
           self.Status.XCatPreferedCoordinate = True

           del Halos, n

           print("Number of Halos = %i" % self.number_of_halos)
           print("Halos class initialized successfully.")

   def SolveLxTxFlux(self):

       if self.Status.HalosDataExist:
           self.Solvers.LxTxSolver(self)
           self.Status.LxTxSolved = True
       else:
           print("ERROR: Halos class is not initialized! ")
           print("Please add a halos catalog! ")
           exit(1)

   def SaveHalosData(self):

       from Objects import FDIR_HALOS, FDIR_CLUSTERS

       # Cheking for error
       if not self.Status.HalosDataExist:
           print("ERROR: Halos class is not initialized! ")
           print("Please add a halos catalog! ")
           exit(1)

       if not self.Status.LxTxSolved:
           print("ERROR: Lx, Tx, and flux was not created! ")
           print("Please solve Lx, Tx, and flux for the halo catalog! ")
           exit(1)

       fdir = FDIR_CLUSTERS + '%s.fit' % (self.fname)
       print("Removing old halo catalog (if exist) : %s" % fdir)
       os.system('rm -r -f %s' % fdir)
       print("Saving new halo catalog : %s" % fdir)

       hdr = pyfits.Header()
       self._setHeader(hdr)
       PrimaryHDUList = np.arange(1)

       # HALOS INFO
       T = np.array(self.lgT)
       Lx = np.array(self.lgLx)
       F = self.InputParam.xray_band.fluxFac * np.power(10., np.array(self.lgFx))
       FU = self.InputParam.xray_band.fluxFacStr
       col1 = pyfits.Column(name='HALOID',
                            format='K', array=np.array(self.ID))
       col2 = pyfits.Column(name='RA', unit='degree',
                            format='E', array=np.array(self.RA))
       col3 = pyfits.Column(name='DEC', unit='degree',
                            format='E', array=np.array(self.DEC))
       col4 = pyfits.Column(name='R500', unit='Mpc/h',
                            format='E', array=np.array(self.R500))
       col5 = pyfits.Column(name='M500', unit='Msun/h',
                            format='E', array=np.array(self.M500))
       col6 = pyfits.Column(name='Z',
                            format='E', array=np.array(self.Z_red))
       col7 = pyfits.Column(name='T', unit='keV',
                            format='E', array=T)
       col8 = pyfits.Column(name='Lx', unit='ergs/s',
                            format='D', array=Lx)
       col9 = pyfits.Column(name='FLUX', unit=FU,
                            format='E', array=F)
       col10 = pyfits.Column(name='DIS', unit='Mpc/h',
                             format='E', array=np.array(self.pd))
       col11 = pyfits.Column(name='XCAT_RA', unit='degree',
                             format='E', array=np.array(self.XCatRA))
       col12 = pyfits.Column(name='XCAT_DEC', unit='degree',
                             format='E', array=np.array(self.XCatDEC))
       col13 = pyfits.Column(name='RC', unit='degree',
                             format='E', array=np.array(self.Rc))
       col14 = pyfits.Column(name='BETA',
                             format='E', array=np.array(self.beta))
       cols = pyfits.ColDefs([col1, col2, col3, col4, col5, col6,
                              col7, col8, col9, col10, col11, col12,
                              col13, col14])

       hdu = pyfits.PrimaryHDU(data=PrimaryHDUList, header=hdr)
       tbhdu_halo = pyfits.BinTableHDU.from_columns(cols)

       hdu = [hdu, tbhdu_halo]
       for i in range(2, len(pyfits.open(FDIR_HALOS + '%s.fit' % self.fname))):
           tbhdu_extra = pyfits.open(FDIR_HALOS + '%s.fit' % self.fname)[2]
           hdu += [tbhdu_extra]

       thdulist = pyfits.HDUList(hdu)

       thdulist.writeto(fdir)

       return 0

   def _setHeader(self, hdr):

       if self.Status.XCatPreferedCoordinate:
           hdr.set('XCATPREF', 'TRUE', 'XTRA prefered coordinate ')
           hdr.set('RA_REF', str(self.RAref), 'XTRA prefered coordinate ')
           hdr.set('DEC_REF', str(self.DECref), 'XTRA prefered coordinate ')
       else:
           hdr.set('XCATPREF', 'FALSE', 'XTRA prefered coordinate transformation')
       if self.Status.LxTxSolved:
           hdr.set('XRAYINFO', 'TRUE', 'whether xray information exist.')
           hdr.set('XRAYBAND', self.InputParam.xray_band.Str, 'unit:keV')
       else:
           hdr.set('XRAYINFO', 'FALSE', 'whether xray information exist.')
       hdr.set('XRAY_MAP', 'FALSE')
       hdr.set('SZ_MAP', 'FALSE')
       # COSMOLOGY
       hdr.add_comment('------- COSMOLOGY --------')
       hdr.set('H0', str(self.InputParam.h_0))
       hdr.set('OMEGA_DE', str(self.InputParam.Omega_DE))
       hdr.set('OMEGA_M', str(self.InputParam.Omega_M))
       hdr.set('OMEGA_b', str(self.InputParam.Omega_b))
       hdr.set('OMEGA_R', str(self.InputParam.Omega_R))
       hdr.set('sigma_8', str(self.InputParam.sigma_8))
       hdr.set('w', str(self.InputParam.w))
       hdr.set('ns', str(self.InputParam.ns))
       hdr.add_comment('This table is created by XTRA version %s' % __version__)


class Map_Halos_Class:

   def __init__(self):

       self.DIS  = []; self.DEC = []; self.F = []
       self.R500 = []; self.RA  = []; self.Z = []
       self.Rc   = []; self.beta = []
       self.RA_ref = 0.; self.DEC_ref = 0.
       self.XCat_prefered = False; self.halosDataExist= False
 
   def update(self, RACAT=0.0, DECCAT=0.0, Index=None, fname=None):

       from index2coordinateReference import index2coordinateReference
       from Objects import FDIR_CLUSTERS

       if fname is None:
           if Index is not None:
               RACAT, DECCAT = index2coordinateReference(index=Index)
           fname = 'XTRA_RA_%0.2f_DEC_%0.2f_'%(RACAT, DECCAT) + 'Xray'

       clusters = pyfits.open(FDIR_CLUSTERS+fname)[1].data
       hdr = pyfits.open(FDIR_CLUSTERS+fname)[0].header

       try:
             self.RA = clusters['XCAT_RA']
             self.DEC = clusters['XCAT_DEC']
             self.R500 = clusters['R500']
             self.Z = clusters['Z']
             self.F = clusters['FLUX']
             self.Rc = clusters['RC']
             self.beta = clusters['BETA']
             self.RA_ref = float(hdr['RA_REF'])
             self.DEC_ref = float(hdr['DEC_REF'])
             self.fname = fname[:-4]
             self.halosDataExist = True
             self.XrayBandstr = hdr['XRAYBAND']
             print("Number of clusters: %i"%len(clusters[:]['Z']))
       except KeyError:
             print("The input table has some problem. ")
             print("It is not consistent with XTRA table structure.")
             self.halosDataExist= False

