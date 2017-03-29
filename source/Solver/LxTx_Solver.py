import numpy as np


# Evolutionary Factor
def Evolution_factor(Cosmology, zred):
   Ez = np.sqrt(Cosmology.Omega_M*(1.0+zred)**3 + Cosmology.Omega_DE)
   return Ez


# h(z)
def Hubz(Cosmology, zred):
   Ez = Evolution_factor(Cosmology, zred)
   return Ez # Cosmology.h_0*Ez


class LxTx_Solver:

    def __init__(self):
        pass

    def solve(self, Halos):

        import numpy.random as npr
        from Objects import RtoD

        n_data = Halos.number_of_halos
        IP = Halos.InputParam

        # --- set up Lx , Tx scaling relations
        LxS = IP.LxScaling # Lxm_Parameters
        TxS = IP.TxScaling # Txm_Parameters

        # calcualte evolution factor
        E = Evolution_factor(IP, Halos.Z_red)

        # calculate mean-log-Lx + scatter
        Ep = Evolution_factor(IP, LxS.z_p)
        Mp = LxS.M_p * IP.h_0
           
        Halos.lgLx = np.log(LxS.Norm) + LxS.M_slope * np.log(Halos.M500/Mp) +\
                     LxS.E_slope * np.log(E/Ep) + npr.normal(0.0, LxS.sig, n_data)
        Halos.lgLx /= np.log(10.0)
        Halos.lgLx += - 44.

        # calculate mean-log-Tx + scatter
        Ep = Evolution_factor(IP, TxS.z_p)
        Mp = TxS.M_p * IP.h_0

        Halos.lgT = np.log(TxS.Norm) + TxS.M_slope*np.log(Halos.M500/Mp) +\
                    TxS.E_slope*np.log(E/Ep) + npr.normal(0.0,TxS.sig,n_data)
        Halos.lgT /= np.log(10.0)

        # calculate core radius, flux, and beta
        for i in range(n_data):

            fx = IP.xray_band.Luminosity2FluxWithCube(Halos.Z_red[i], Halos.lgT[i], Halos.lgLx[i])
            Halos.lgFx[i] = np.log10(fx+1e-40)

            Rcbar = RtoD * (IP._param['xc_bar']*Halos.R500[i]) * (1.0 + Halos.Z_red[i]) / Halos.pd[i]
            Halos.Rc[i] = npr.lognormal(np.log(Rcbar), IP._param['xc_sig'])

            betaBarC = np.log(IP._param['SB_beta_bar']) +\
                       IP._param["xc_Beta_r"] * IP._param['SB_beta_sig'] / IP._param['xc_bar'] *\
                       (np.log(Halos.Rc[i]) - np.log(Rcbar))

            sigBeta = np.sqrt(1.0 - IP._param["xc_Beta_r"]**2) * IP._param['SB_beta_sig']
            Halos.beta[i] = npr.lognormal(betaBarC, sigBeta)

        print "Fluxes, Tempratures, and Luminosities are assigned successfully!"

        """
        dis2 = (Halos.pd/IP.h_0)**2
        for i in range(n_data):

            fx = IP.xray_band.Luminosity2FluxWithKcorrection( Halos.Z_red[i],Halos.lgT[i],Halos.lgLx[i],dis2[i]) * 6.2415 * 1e8 * 1387.71
        """