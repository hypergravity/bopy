# -*- coding: utf-8 -*-
from __future__ import print_function

"""

Author
------
Bo Zhang

Email
-----
bozhang@nao.cas.cn

Created on
----------
- Sat Jul 16 22:30:00 2016

Modifications
-------------
- Sat Jul 16 22:30:00 2016    re-format code

Aims
----
- to generate config files of STARLIGHT

"""


import numpy as np

__extra_comments__ = '''
This grid file is written using bopy.starlight
(BOPY, Bo Zhang's python package, https://github.com/hypergravity/bopy)
'''


class StarlightGrid(object):
    """ StarlightGrid class is to represent the Grid file object for STARLIGHT
    """
    # specified order of meta data
    meta_order = ['num_of_fits_to_run',
                  'base_dir',
                  'obs_dir',
                  'mask_dir',
                  'out_dir',
                  'phone_number',
                  'llow_SN',
                  'lupp_SN',
                  'Olsyn_ini',
                  'Olsyn_fin',
                  'Odlsyn',
                  'fscale_chi2',
                  'fit_fxk',
                  'IsErrSpecAvailable',
                  'IsFlagSpecAvailable']
    # default values for StarlightGrid instance
    meta = dict(num_of_fits_to_run=2,
                base_dir='/pool/projects/starlight/STARLIGHTv04/BasesDir/',
                obs_dir='/pool/projects/starlight/STARLIGHTv04/',
                mask_dir='/pool/projects/starlight/STARLIGHTv04/',
                out_dir='/pool/projects/starlight/STARLIGHTv04/',
                phone_number='-2007200',
                llow_SN=4730.0,
                lupp_SN=4780.0,
                Olsyn_ini=3400.0,
                Olsyn_fin=8900.0,
                Odlsyn=1.0,
                fit_fxk='FIT',
                fscale_chi2=1.0,
                IsErrSpecAvailable=1,
                IsFlagSpecAvailable=1)
    # specified order of arq
    arq_order = ['arq_obs',
                 'arq_config',
                 'arq_base',
                 'arq_masks',
                 'red_law_option',
                 'v0_start',
                 'vd_start',
                 'arq_out']
    # default arq data
    arq_obs = []
    arq_config = []
    arq_base = []
    arq_masks = []
    red_law_option = []
    v0_start = []
    vd_start = []
    arq_out = []
    # extra comments
    extra = __extra_comments__

    def __init__(self, **kwargs):
        """ initialize instance using arq data """
        for key in kwargs.keys():
            self.__setattr__(key, kwargs[key])

    def set_meta(self, **kwargs):
        """ set meta data """
        for key in kwargs.keys():
            self.meta[key] = kwargs[key]

    def is_arq_validate(self):
        """ check if the arq data length validate """
        try:
            n_obs = len(self.arq_obs)
            assert n_obs == len(self.arq_config) \
                or np.isscalar(self.arq_config)
            assert n_obs == len(self.arq_base) \
                or np.isscalar(self.arq_base)
            assert n_obs == len(self.arq_masks) \
                or np.isscalar(self.arq_masks)
            assert n_obs == len(self.red_law_option) \
                or np.isscalar(self.red_law_option)
            assert n_obs == len(self.v0_start) \
                or np.isscalar(self.v0_start)
            assert n_obs == len(self.vd_start) \
                or np.isscalar(self.vd_start)
            assert n_obs == len(self.arq_out)
            return True
        except AssertionError:
            return False

    def pprint_meta(self):
        """ print meta data """
        for key in self.meta_order:
            print("%s: %s" % (key, self.meta[key]))

    def _print_arq_(self, arq_key):
        """ print single arq field data """
        assert arq_key in self.arq_order

        arq_val = self.__getattribute__(arq_key)

        if np.isscalar(arq_val) or len(arq_val) <= 3:
            print(arq_val)
        else:
            print('[%s, %s, ... %s]' % (arq_val[0], arq_val[1], arq_val[-1]))

    def pprint_arq(self):
        """ print all arq data"""
        for key in self.arq_order:
            print("%s:" % key),
            self._print_arq_(key)

    def pprint(self):
        """ print a summary of the instance """
        print("")
        print('StarlightGrid summary:')
        print('############### [meta] ###############')
        self.pprint_meta()
        print('############### [arq]  ###############')
        self.pprint_arq()
        print(self.extra)
        print('######################################')

    def _meta_to_string(self, val_width=50):
        """ convert meta data to string list """
        fmt_str = "%%%ds" % (-val_width)
        return [(fmt_str + "[%s]\n") % (self.meta[key], key)
                for key in self.meta_order]

    def _arq_to_string(self, sep='   '):
        """ convert arq data to string list """
        #  1. arq data to list
        n_obs = len(self.arq_obs)
        for key in self.arq_order:
            val = self.__getattribute__(key)
            if np.isscalar(val):
                self.__setattr__(key, [val for i in xrange(n_obs)])
            else:
                assert len(val) == n_obs

        #  2. to string
        str_list = []
        for i in xrange(n_obs):
            arq_data_i = ["%s" % self.__getattribute__(key)[i]
                          for key in self.arq_order]
            str_list.append(sep.join(arq_data_i) + "\n")
        return str_list

    def write(self, filepath, meta_val_width=50, sep='   '):
        f = open(filepath, "w+")
        f.writelines(self._meta_to_string(meta_val_width))
        f.write('\n')
        f.writelines(self._arq_to_string('   '))
        f.write('\n')
        f.write(self.extra)
        f.close()


def _test_starlight_grid():
    sg = StarlightGrid(arq_obs=['0414.51901.393.cxt',
                                '0784.52327.478.cxt'],
                       arq_config='StCv04.C99.config',
                       arq_base='Base.BC03.N',
                       arq_masks=['Mask.0414.51901.393.cxt.sc1.CRAP.gm.BN',
                                  'Mask.0784.52327.478.cxt.sc2.CRAP.gm.BN'],
                       red_law_option='CCM',
                       v0_start=0.0,
                       vd_start=150.0,
                       arq_out=['0414.51901.393.cxt.sc4.C99.im.CCM.BN',
                                '0784.52327.478.cxt.sc4.C99.im.CCM.BN']
                       )
    sg.pprint_meta()
    sg.pprint_arq()
    sg.pprint()
    for s in sg._meta_to_string():
        print(s)
    for s in sg._arq_to_string(',:'):
        print(s)
    sg.write('/pool/projects/starlight/STARLIGHTv04/grid_example1.in_')


def starlight_config():
    pass


def starlight_base():
    pass


def starlight_mask():
    pass


if __name__ =='__main__':
    _test_starlight_grid()