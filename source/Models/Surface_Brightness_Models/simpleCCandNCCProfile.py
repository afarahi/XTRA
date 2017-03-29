import numpy as np
import random as rd
from XCat_Objects import DtoR, RtoD

class Surface_Brightness_Model():

   def __init__(self):

      # For CC clusters
      self.beta = 0.666
      self.power = -3.0*self.beta+0.5
      self.fc = 0.2 # fc = Rc/R500 

      # For NCC clusters
      self.beta_1 = 0.666
      self.power_1 = -3.0*self.beta_1+0.5
      self.fc_1 = 0.2 # fc = Rc/R500 

      self.beta_2 = 0.666
      self.power_2 = -3.0*self.beta_2+0.5
      self.fc_2 = 0.2 # fc = Rc/R500 

      self.fS = 1.0 # fS = S2/S1 

      # Fraction of CC clusters
      self.fCC = 0.5 

   def flux2SurfaceBrightness(self,Map,Halos,samplei):

      if (rd.random() < fCC):

         RA    = Halos.RA[samplei]
         DEC   = Halos.DEC[samplei]
  
         Rc_1 = self.fc_1*Halos.R500[samplei]
         thetac2_1 = ( RtoD * Rc_1 * (1.0+Halos.Z[samplei]) / Halos.DIS[samplei] )**2
         Rc_2 = self.fc_2*Halos.R500[samplei]
         thetac2_2 = ( RtoD * Rc_2 * (1.0+Halos.Z[samplei]) / Halos.DIS[samplei] )**2

         (pixi,pixj) = Map.pixInsideArea(RA,DEC,10.0*np.sqrt(thetac2_1))
         pixValue = np.zeros(len(pixi))
         summ = 0.0
         sigmaRA = Map.dRA/2.0
         sigmaDEC= Map.dDEC/2.0
  
         print "Sample Number : ", len(pixi)
         for i in range(len(pixi)):
            (sRA,sDEC) = Map.pix2ang(pixi[i],pixj[i]) 
            stheta2  = ((sRA - RA)**2 + (sDEC-DEC)**2)
            pixValue[i] = ( 1.0 + stheta2/thetac2_1 )**self.power_1
                        + fS * ( 1.0 + stheta2/thetac2_2 )**self.power_2
         summ = np.sum(pixValue[:])
         pixValue = Halos.F[samplei] * pixValue / summ / Map.MPa.effArea
         for i in range(len(pixValue[:])): 
            if (pixi[i] < 0 or pixj[i] < 0 
               or pixi[i] > Map.nside-1 
               or pixj[i] > Map.nside-1): pass
            else: Map.MAP[pixi[i],pixj[i]] += pixValue[i] 

      else:
         RA    = Halos.RA[samplei]
         DEC   = Halos.DEC[samplei]

         Rc = self.xc*Halos.R500[samplei]
         thetac2 = ( RtoD * Rc * (1.0+Halos.Z[samplei]) / Halos.DIS[samplei] )**2
         (pixi,pixj) = Map.pixInsideArea(RA,DEC,10.0*np.sqrt(thetac2))
         pixValue = np.zeros(len(pixi))
         summ = 0.0
         sigmaRA = Map.dRA/2.0
         sigmaDEC= Map.dDEC/2.0

         print "Sample Number : ", len(pixi)
         for i in range(len(pixi)):
            (sRA,sDEC) = Map.pix2ang(pixi[i],pixj[i]) 
            stheta2  = ((sRA - RA)**2 + (sDEC-DEC)**2)
            pixValue[i] = ( 1.0 + stheta2/thetac2 )**self.power
         summ = np.sum(pixValue[:])
         pixValue = Halos.F[samplei] * pixValue / summ / Map.MPa.effArea
         for i in range(len(pixValue[:])): 
            if (pixi[i] < 0 or pixj[i] < 0 
               or pixi[i] > Map.nside-1 
               or pixj[i] > Map.nside-1): pass
            else: Map.MAP[pixi[i],pixj[i]] += pixValue[i] 

