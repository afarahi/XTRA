from Objects import fluxfac, ln10, ergs2keV
import numpy as np
import json
import sys


def fluxUnitCanvFac(EnergyUnitStr = 'ergs'):

    from Objects import ergs2keV

    if (EnergyUnitStr == 'ergs' or EnergyUnitStr == 'erg' or EnergyUnitStr == 'ERGS' or EnergyUnitStr == 'ERG'):
        return 1, 'ergs/s/cm^2'
    elif (EnergyUnitStr == 'keV' or EnergyUnitStr == 'kev' or EnergyUnitStr == 'KeV' or EnergyUnitStr == 'KEV'):
        return ergs2keV, 'keV/s/cm^2'
    else:
        print('WARNING: Unit %s was not implemented!'%EnergyUnitStr)
        print('It uses the default energy unit which is is ergs.')
        return 1, 'ergs/s/cm^2'


class Xray_Band_Class:

    def __init__(self, option=1):

        with open('parameters/Input_Parameters.json') as fp:
            EnergyUnitStr = json.load(fp)['Energy_Unit']
        self.fluxFac, self.fluxFacStr = fluxUnitCanvFac(EnergyUnitStr=EnergyUnitStr)

        self.Conv_fac = 1.0  # 0.62
        self.Str = '[0.5-2.0]'
        self.opt_num = 1

        if (option == 1):
            # K-Correction function
            self.K_corr = self.K_corr_1

        elif (option == 2):
            self.Conv_fac = 1.
            self.Str = '[0.1-2.4]'
            self.opt_num = 2
            # K-Correction function
            self.K_corr = self.K_corr_2

        elif (option == 3):
            import astropy.io.fits as pyfits
            fdir = './parameters/Models/L_to_Flux/'
            with open('parameters/Input_Parameters.json') as fp:
                _parm = json.load(fp)
                fname = _parm['L_to_Flux_model_file_name']
                self.Str = _parm['Flux_energy_band']

            x = pyfits.open(fdir + fname)[1].data
            self.lgTCUBE = np.log10(x[0]['T'])
            self.ZCUBE = x[0]['Z']
            # self.NHCUBE= x[0]['NH']
            # self.L_CUBE = x[0]['L_CUBE'][:][:][:]
            self.L_CUBE = x[0]['CR_PERL_TIMESD2'][:][:]
            del x

        else:
            print("ERROR: Input XRay band option does not exist!")
            print("Please change the JSON file and try again.")
            sys.exit(2)

    # --- K-Correction --> K_corr = Lobs / Lrest
    # For option = 1 (Nord et al. 2008 --- arXiv:0706.2189)
    # This fit is accurate to a few percent within z=2. (T in keV)
    def K_corr_1(self, Z, lgT):
        K1 = -0.209 + lgT * (1.18 - 0.39 * lgT)
        K2 = -0.098 + lgT * (-0.092 + 0.085 * lgT)
        return 1. + Z * (K1 + K2 * Z)

    # For option = 2 (Stanek et al. 2005 --- arXiv:astro-ph/0602324)
    # This is for band [0.1-2.4] keV (equation 4 in paper). (T in keV)
    # --- set K_corr = sqrt(1+(1+log10(Tx/5))*z)
    def K_corr_2(self, Z, lgT):
        return np.sqrt(1. + (1. + lgT - 0.69897) * Z)

    def fluxCalculator(self, r2, z_redshift, lglobs):
        # convert comoving distance to Mpc for assumed h0
        dlMpc2 = r2 * (1.0 + z_redshift) ** 2
        fx = fluxfac * (10 ** lglobs) / dlMpc2
        return fx

    def Luminosity2FluxWithKcorrection(self, Z, lgT, lgLx, comdis2):

        import astropy.io.fits as pyfits

        fdir = './parameters/Models/L_to_Flux/'
        x = pyfits.open(fdir + 'L_0.5-2_L_bol_factor.fits')[1].data

        lgTCUBE = np.log10(x[0]['T'])
        fac = x[0]['LBAND_LBOL']
        del x
        iT = int(float(len(lgTCUBE) - 1) * (lgT - lgTCUBE[0]) / (lgTCUBE[-1] - lgTCUBE[0]))

        T2 = lgTCUBE[iT + 1]
        T1 = lgTCUBE[iT]
        fac1 = fac[iT]
        fac2 = fac[iT + 1]
        fac = np.log10((fac1 * (T2 - lgT) + fac2 * (lgT - T1)) / (T2 - T1))
        lgLx = fac + lgLx

        # --- set Lx[0.5-2.0] = 0.62 Lx[0.1-2.4]
        # --- set Lx[0.5-2.0] = Lx[0.1-2.4] * conv_fac
        # --- set K_corr = Lx_obs / Lx_rest
        LobsLrest = self.K_corr(Z, lgT)
        lglobs = lgLx + np.log10(self.Conv_fac * LobsLrest)
        return self.fluxCalculator(comdis2, Z, lglobs)

    def LCUBE_INDEX(self, Z=0.1, lgT=0.0, iNH=0):

        iZ = int(float(len(self.ZCUBE) - 1) * (Z - self.ZCUBE[0]) / (self.ZCUBE[-1] - self.ZCUBE[0]))
        ilgT = int(float(len(self.lgTCUBE) - 1) * (lgT - self.lgTCUBE[0]) / (self.lgTCUBE[-1] - self.lgTCUBE[0]))
        return iZ, ilgT, iNH

    def Luminosity2FluxWithCube(self, Z, lgT, lgLx):

        iZ, ilgT, iNH = self.LCUBE_INDEX(Z=Z, lgT=lgT)

        Z2 = self.ZCUBE[iZ + 1]
        Z1 = self.ZCUBE[iZ]

        T2 = self.lgTCUBE[ilgT + 1]
        T1 = self.lgTCUBE[ilgT]

        L11 = self.L_CUBE[iZ][ilgT]
        L21 = self.L_CUBE[iZ + 1][ilgT]
        L12 = self.L_CUBE[iZ][ilgT + 1]
        L22 = self.L_CUBE[iZ + 1][ilgT + 1]

        CRPERLx = (L11 * (Z2 - Z) * (T2 - lgT) + L21 * (Z - Z1) * (T2 - lgT) +
                   L12 * (Z2 - Z) * (lgT - T1) + L22 * (Z - Z1) * (lgT - T1)) / ((Z2 - Z1) * (T2 - T1))

        return CRPERLx * np.power(10.0, lgLx)