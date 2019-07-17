import numpy as np
from Objects import DtoR, RtoD

class Surface_Brightness_Model():

   def __init__(self):
      self.clusterEffRadCOEFF = 20.0
      pass


   def flux2SurfaceBrightness(self,Map,Halos,samplei,MPa):
 
      pixelSide = MPa.PixSide 
      mapCenterShift = MPa.DECMapSize/2. # in degree
      mapSize = MPa.DECMapSize # in degree

      RApix   = (Halos.RA[samplei]+mapCenterShift)*pixelSide/mapSize
      DECpix  = (Halos.DEC[samplei]+mapCenterShift)*pixelSide/mapSize
      rad     = Halos.Rc[samplei]*pixelSide/mapSize
      effRad  = self.clusterEffRadCOEFF*rad
      thetac2 = rad**2

      bpower = -3.0*Halos.beta[samplei]+0.5

      iMax = int(RApix+effRad);   iMin = int(RApix-effRad)
      jMax = int(DECpix+effRad);  jMin = int(DECpix-effRad)
      # if the halo wouldn't show on the plot; don't bother. 
      if (iMin>=pixelSide-1 or jMin>=pixelSide-1): return 0
      if (iMax<=0 or jMax<=0): return 0

      i,j = np.mgrid[iMin:iMax,jMin:jMax]
      if 1: # use original method of all circular halos
        mask = (i-RApix)**2 + (j-DECpix)**2 < effRad**2
        mapA = mask*( 1. + ((RApix - i)**2 + (DECpix - j)**2)/thetac2 )**bpower
      else: # generate random orientation and ellipticity for each halo
        thetaOrient = np.random.random()*2*np.pi
        # ellipticity = np.random.random()*0.5
        ellipticity = 0.5
        x = i-RApix; y = j-DECpix; 
        coords = (x*np.cos(thetaOrient)+y*np.sin(thetaOrient))**2 \
               + (x*np.sin(thetaOrient)-y*np.cos(thetaOrient))**2/(1-ellipticity**2)
        mask = coords < effRad**2
        mapA = mask*( 1. + coords/thetac2)**bpower
        del thetaOrient, ellipticity, x, y, coords
      
      mapA *= Halos.F[samplei] / sum(sum(mapA))

      iMaxA = min(iMax,pixelSide)
      iMinA = max(iMin,0)
      jMaxA = min(jMax,pixelSide)
      jMinA = max(jMin,0)

      Map.MAP[iMinA:iMaxA,jMinA:jMaxA] += \
            mapA[(iMinA-iMin):(iMaxA-iMin),(jMinA-jMin):(jMaxA-jMin)]
