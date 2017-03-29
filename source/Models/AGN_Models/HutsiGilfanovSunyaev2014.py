def flux_cal(x,y,z,z_redshift,h_0,lnlobs): 
  from XCat_Objects import fluxfac
  #! --- convert comoving distance to Gpc for assumed h0
  r2     = ( x**2 + y**2 + z**2 ) / (1000.0*h_0)**2
  dlGpc2 = r2 * (1.0+z_redshift)**2
  fx = fluxfac * exp(lnlobs) / dlGpc2
  return fx

def flux_calculator(r2,Z,lobs): 
  from XCat_Objects import pi, Mpc2cm 
  fac = ((1./Mpc2cm)*(1./Mpc2cm)/(4.0*pi))
  #! --- convert comoving distance to Mpc for assumed h0
  dlMpc2 = r2 * (1.0+Z)**2
  fx = fac * lobs / dlMpc2
  return fx

def f_duty(z):  return ( 0.47**2*z/((z-0.83)**2+0.47**2) + 0.11 )
def Gamma_z(z): return ( 0.84 - 0.21*z )

def HutsiGilfanovSunyaev2014(AGNs,Halos):

  from XCat_Objects import fluxfac, ln10
  import random as rd
  import numpy  as np

  n_Halos = Halos.number_of_halos

  L0 = 1e41
  M0 = 27.3

  IP = Halos.InputParam

  AGNID = 1

  for i in range(n_Halos):

    rand_var = rd.random()
    #Following Hutsi Gilfanov Sunyaev 2014 (Model 1 f_duty <= 1)
    if ( rand_var < f_duty(Halos.Z_red[i]) ):

       AGNs.AGNID.append(AGNID)
       AGNs.HALOID.append(Halos.ID[i])
       AGNs.Z.append(Halos.Z_red[i])
       AGNs.RA.append(Halos.RA[i])
       AGNs.DEC.append(Halos.DEC[i])
       AGNs.XCatRA.append(Halos.XCatRA[i])
       AGNs.XCatDEC.append(Halos.XCatDEC[i])

       # (eq. 1) for Lx[2-10]
       Lx = np.power(Halos.M500[i]/M0,Gamma_z(Halos.Z_red[i])) * L0
       # Appling K-correction
       LobsLrest = 1. #IP.xray_band.K_corr(Halos.Z_red[i],lgT[i])
       Lobs = LobsLrest * Lx
       dis2 = (Halos.pd[i]/IP.h_0)**2
       flux = flux_calculator(dis2,Halos.Z_red[i],Lobs)

       AGNs.LxAGN.append(Lx)
       AGNs.fluxAGN.append(flux)

       AGNID += 1

       print "AGN flux = ", flux

  print "AGNs are created successfully!" 

