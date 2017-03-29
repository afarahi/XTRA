import json
import sys


class Input_Parameters:

    def __init__(self):

        from Objects import Xray_Band_Class
        from Objects import Temprature_scaling, Luminocity_scaling

        with open('parameters/Cosmological_Parameters.json') as fp:
            self._param = json.load(fp)

        # Cosmology
        self.h_0 = self._param['Hubble_Parameter']
        self.Omega_DE = self._param['Omega_DE']
        self.Omega_M = self._param['Omega_M']
        self.Omega_b = self._param['Omega_b']
        self.Omega_R = self._param['Omega_R']
        self.Omega_k = self._param['Omega_k']
        self.sigma_8 = self._param['sigma_8']
        self.w = self._param['w']
        self.ns = self._param['ns']

        with open('parameters/Input_Parameters.json') as fp:
            self._param = json.load(fp)

        # Models
        self.Lxm_mod = self._param['Lxm_model']
        self.Txm_mod = self._param['Txm_model']
        xray_band_option = 3  # self._param['XRay_band_Option']

        self.xray_band = Xray_Band_Class(option=xray_band_option)

        try:
            self.LxScaling = Luminocity_scaling(self.Lxm_mod)
        except KeyError:
            print("ERROR: Input Lxm Parameter method does not exist!")
            print("Please modify JSON file and try again.")
            sys.exit(2)

        try:
            self.TxScaling = Temprature_scaling(self.Txm_mod)
        except KeyError:
            print("ERROR: Input Txm Parameter method does not exist!")
            print("Please modify JSON file and try again.")
            sys.exit(2)

        print("COSMOLOGY and INPUTS are assigned SUCCESSFULLY.")
