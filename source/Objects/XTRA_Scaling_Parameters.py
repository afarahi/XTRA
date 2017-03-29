import os.path
import json


class Temprature_scaling:

    def __init__(self, label):

        fname = './parameters/Models/Txm/' + label + '_parameters.json'

        if os.path.isfile(fname) == False:
            print("Error: %s does not exists it uses Tx scaling default parameters."%s)
            exit(1)
            # fname = './parameters/Models/Txm/default_parameters.xml'

        with open(fname) as fp:
            _param = json.load(fp)

        # Parameters
        self.Norm = _param['a']
        self.M_slope = _param['M_slope']
        self.E_slope = _param['E_slope']

        self.M_p = _param['M_p']
        self.z_p = _param['z_p']

        self.sig = _param['sig']


class Luminocity_scaling:

    def __init__(self, label):

        fname = './parameters/Models/Lxm/' + label + '_parameters.json'

        if os.path.isfile(fname) == False:
            print("ERROR: %s does not exists it uses Lx scaling default parameters." % s)
            exit(1)
            # fname = './parameters/Models/Lxm/default_parameters.xml'

        with open(fname) as fp:
            _param = json.load(fp)

        # Parameters
        self.Norm = _param['a']
        self.M_slope = _param['M_slope']
        self.E_slope = _param['E_slope']

        self.M_p = _param['M_p']
        self.z_p = _param['z_p']

        self.sig = _param['sig']








