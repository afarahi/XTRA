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
      if (iMin>=pixelSide-1 or jMin>=pixelSide-1): return 0
      if (iMax<=0 or jMax<=0): return 0

      i,j = np.mgrid[iMin:iMax,jMin:jMax]
      mask = (i-RApix)**2 + (j-DECpix)**2 < effRad**2
      mapA = mask*( 1. + ((RApix - i)**2 + (DECpix - j)**2)/thetac2 )**bpower
      mapA *= Halos.F[samplei] / sum(sum(mapA))

      iMaxA = min(iMax,pixelSide)
      iMinA = max(iMin,0)
      jMaxA = min(jMax,pixelSide)
      jMinA = max(jMin,0)

      Map.MAP[iMinA:iMaxA,jMinA:jMaxA] += \
            mapA[(iMinA-iMin):(iMaxA-iMin),(jMinA-jMin):(jMaxA-jMin)]
