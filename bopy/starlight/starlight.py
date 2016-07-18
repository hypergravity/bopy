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
- Sat Jul 16 22:30:00 2016    framework
- Sun Jul 17 23:00:00 2016    StarlightGrid
- Mon Jul 18 09:27:00 2016

Aims
----
- to generate config files of STARLIGHT

"""


import numpy as np

__extra_comments__ = '''
# This grid file is written using bopy.starlight
# BOPY, Bo Zhang's python package, https://github.com/hypergravity/bopy
'''


# ################
# Starlight Grid #
# ################


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


# ################
# Starlight Base #
# ################


class StarlightBase(object):
    """ StarlightBase class is to represent the Base file object for STARLIGHT
    """
    # specified order of meta data
    meta_order = ['n_base']
    # default values for StarlightBase instance
    meta = dict(n_base=45)
    # specified order of arq data
    arq_order = ['spec_file',
                 'age',
                 'z',
                 'code',
                 'mstar',
                 'yav',
                 'afe']
    # default values for bases
    spec_file = []
    age = []
    z = []
    code = []
    mstar = []
    yav = []
    afe = []
    # extra comments
    extra = __extra_comments__

    def __init__(self, **kwargs):
        """ initialize instance using arq data """
        for key in kwargs.keys():
            self.__setattr__(key, kwargs[key])
        # count spec_file
        self.meta['n_base'] = len(self.spec_file)

    def _meta_to_string(self, val_width=50):
        """ convert meta data to string list """
        fmt_str = "%%%ds" % (-val_width)
        return [(fmt_str + "[%s]\n") % (self.meta[key], key)
                for key in self.meta_order]

    def _arq_to_string(self, sep='   '):
        """ convert arq data to string list """
        #  1. arq data to list
        n_obs = len(self.spec_file)
        for key in self.arq_order:
            val = self.__getattribute__(key)
            if np.isscalar(val):
                self.__setattr__(key, [val for i in xrange(n_obs)])
            else:
                assert len(val) == n_obs

        # 2. to string
        str_list = []
        for i in xrange(n_obs):
            arq_data_i = ["%s" % self.__getattribute__(key)[i]
                          for key in self.arq_order]
            str_list.append(sep.join(arq_data_i) + "\n")
        return str_list

    def write(self, filepath, meta_val_width=50, sep='   '):
        f = open(filepath, "w+")
        f.writelines(self._meta_to_string(meta_val_width))
        # f.write('\n')
        f.writelines(self._arq_to_string('   '))
        f.write('\n')
        f.write(self.extra)
        f.close()


def _test_starlight_base():
    sg = StarlightBase(spec_file=['bc2003_hr_m42_chab_ssp_020.spec',
                                  'bc2003_hr_m42_chab_ssp_045.spec'],
                       age=[0.00100e9,
                            0.00316e9],
                       z=[0.00400,
                          0.00400],
                       code=['age020_m42',
                             'age020_m42'],
                       mstar=[1.0000,
                              0.9999],
                       yav=[0,
                            0],
                       afe=[0.0,
                            0.0]
                       )
    sg.write('/pool/projects/starlight/STARLIGHTv04/Base.BC03.N_')


_config_StCv04C99_ = dict(
    # Normalization lambdas
    l_norm=4020.0,
    llow_norm=4010.0,
    lupp_norm=4060.0,
    # Parameter Limits
    AV_low=-1.0,
    AV_upp=4.0,
    YAV_low=-0.0001,
    YAV_upp=0.0001,
    fn_low=0.7,
    fn_upp=1.3,
    v0_low=-500.0,
    v0_upp=500.0,
    vd_low=0.0,
    vd_upp=500.0,
    # Clipping options & Weight-Control-Filter
    clip_method_option='NSIGMA', # NOCLIP/NSIGMA/RELRES/ABSRES
    sig_clip_threshold=3.0,
    wei_nsig_threshold=2.0,
    # Miscellaneous
    dl_cushion=50.0,
    f_cut=0.001,
    N_int_Gauss=31,
    i_verbose=1,
    i_verbose_anneal=1,
    Is1stLineHeader=0,
    i_FastBC03_FLAG=1,
    i_FitPowerLaw=0,
    alpha_PowerLaw=-0.5,
    i_SkipExistingOutFiles=0,
    # Markov Chains technical parameters
    N_chains=7,
    xinit_max=0.50,
    i_UpdateEps=0,
    i_UpdateAlpha=2,
    Falpha=2.0,
    i_MoveOneParOnly=1,
    i_UpdateAVYAVStepSeparately=1,
    i_HelpParWithMaxR=1,
    prob_jRmax=0.2,
    i_HelpPopVectorMove2Average=1,
    prob_HelpPopVectorMove2Average=0.4,
    i_HelpAVYAVMove2Average=1,
    prob_HelpAVYAVMove2Average=0.4,
    NRC_AV_Default=10,
    # First Fits (FF) technical parameters
    Temp_ini_FF=1.0e2,
    Temp_fin_FF=1.0,
    fini_eps_FF=1.0e1,
    ffin_eps_FF=1.0e2,
    R_ini_FF=1.3,
    R_fin_FF=1.2,
    IsGRTestHard_FF=0,
    N_loops_FF=10,
    i_RestartChains_FF=1,
    fN_sim_FF=1.0e1,
    fNmax_steps_FF=1.0e4,
    eff_IDEAL_FF=0.23,
    # GR R-threshold & Method for Burn-In loop
    R_Burn_in=1.2,
    IsGRTestHard_BurnIn=1,
    # EX0s technical parameters
    EXOs_PopVector_option='MIN', # MIN/AVE
    EXOs_method_option='CUMUL', # CUMUL/SMALL
    EXOs_Threshold=0.02,
    Temp_ini_EX0s=1.0,
    Temp_fin_EX0s=1.0e-3,
    fini_eps_EX0s=1.0e2,
    ffin_eps_EX0s=1.0e3,
    R_ini_EX0s=1.2,
    R_fin_EX0s=1.0,
    IsGRTestHard_EX0s=1,
    N_loops_EX0s=10,
    i_RestartChains_EX0s=1,
    fN_sim_EX0s=1.0e2,
    fNmax_steps_EX0s=1.0e3,
    eff_IDEAL_EX0s=0.50,
    IsScaleNstepsInEX0sFits=1,
    IsNoKin4LargeBaseInEX0sFits=0,
    frac_NoKin4LargeBaseInEX0sFits=0.0,
    fEX0_MinBaseSize=0.1)
_config_StCv04C11_ = dict(
    # Normalization lambdas
    l_norm=4020.0,
    llow_norm=4010.0,
    lupp_norm=4060.0,
    # Parameter Limits
    AV_low=-1.0,
    AV_upp=4.0,
    YAV_low=-0.0001,
    YAV_upp=0.0001,
    fn_low=0.7,
    fn_upp=1.3,
    v0_low=-500.0,
    v0_upp=500.0,
    vd_low=0.0,
    vd_upp=500.0,
    # Clipping options & Weight-Control-Filter
    clip_method_option='NSIGMA', # NOCLIP/NSIGMA/RELRES/ABSRES
    sig_clip_threshold=3.0,
    wei_nsig_threshold=2.0,
    # Miscellaneous
    dl_cushion=50.0,
    f_cut=0.001,
    N_int_Gauss=31,
    i_verbose=1,
    i_verbose_anneal=0,
    Is1stLineHeader=0,
    i_FastBC03_FLAG=1,
    i_FitPowerLaw=0,
    alpha_PowerLaw=-0.5,
    i_SkipExistingOutFiles=0,
    # Markov Chains technical parameters
    N_chains=7,
    xinit_max=0.50,
    i_UpdateEps=0,
    i_UpdateAlpha=2,
    Falpha=2.0,
    i_MoveOneParOnly=1,
    i_UpdateAVYAVStepSeparately=1,
    i_HelpParWithMaxR=1,
    prob_jRmax=0.2,
    i_HelpPopVectorMove2Average=1,
    prob_HelpPopVectorMove2Average=0.4,
    i_HelpAVYAVMove2Average=1,
    prob_HelpAVYAVMove2Average=0.4,
    NRC_AV_Default=10,
    # First Fits (FF) technical parameters
    Temp_ini_FF=1.0e2,
    Temp_fin_FF=1.0,
    fini_eps_FF=1.0e1,
    ffin_eps_FF=1.0e2,
    R_ini_FF=1.3,
    R_fin_FF=1.3,
    IsGRTestHard_FF=0,
    N_loops_FF=3,
    i_RestartChains_FF=1,
    fN_sim_FF=1.0e1,
    fNmax_steps_FF=1.0e4,
    eff_IDEAL_FF=0.23,
    # GR R-threshold & Method for Burn-In loop
    R_Burn_in=1.2,
    IsGRTestHard_BurnIn=0,
    # EX0s technical parameters
    EXOs_PopVector_option='MIN', # MIN/AVE
    EXOs_method_option='CUMUL', # CUMUL/SMALL
    EXOs_Threshold=0.02,
    Temp_ini_EX0s=1.0,
    Temp_fin_EX0s=1.0e-3,
    fini_eps_EX0s=1.0e2,
    ffin_eps_EX0s=1.0e3,
    R_ini_EX0s=1.2,
    R_fin_EX0s=1.0,
    IsGRTestHard_EX0s=1,
    N_loops_EX0s=5,
    i_RestartChains_EX0s=1,
    fN_sim_EX0s=1.0e2,
    fNmax_steps_EX0s=1.0e3,
    eff_IDEAL_EX0s=0.50,
    IsScaleNstepsInEX0sFits=1,
    IsNoKin4LargeBaseInEX0sFits=1,
    frac_NoKin4LargeBaseInEX0sFits=0.0,
    fEX0_MinBaseSize=0.1)

_config_comments_ = (
    ('# Configuration parameters for StarlightChains_v04.for - Cid@Lagoa - 18/Feb/2007 #\n'
     '#\n#\n# Normalization lambdas\n#\n'),
    '#\n#\n# Parameter Limits\n#\n',
    '#\n#\n# Clipping options & Weight-Control-Filter\n#\n',
    '#\n#\n# Miscellaneous\n#\n',
    '#\n#\n# Markov Chains technical parameters\n#\n',
    '#\n#\n# First Fits (FF) technical parameters\n#\n',
    '#\n#\n# GR R-threshold & Method for Burn-In loop\n#\n',
    '#\n#\n# EX0s technical parameters\n#\n',
    (
        '\n\nCid@Lagoa - 18/February/2007'
        '\n\nOBS: A SLOW config!'
        '\n\n\n\n'
        'Technical parameters you may want to play with to obtain FAST, MEDIUM'
        '\n& SLOW fits:\n\n'
        '--------------------------------\n'
        '|  FAST  |    MEDIUM  |  LONG  |\n'
        '--------------------------------\n'
        '|  5     |    7       |  12    | [N_chains]\n'
        '|  3     |    5       |  10    | [N_loops_FF & *_EX0s]\n'
        '|  1.3   |    1.2     |  1.1   | [R_ini_FF & R_fin_FF & *_EX0s]\n'
        '|  0     |   0 or 1   |  1     | [IsGRTestHard_FF & *_BurnIn & *_EX0s]\n'
        '--------------------------------\n')
    )
_config_comments_insert_index_ = (
    0, 4, 15, 19, 30, 45, 58, 61, 65)


# ##################
# Starlight Config #
# ##################


class StarlightConfig(object):
    """ StarlightConfig class is to represent the Config file object for STARLIGHT
    """
    # specified order of meta data
    meta_order = [
        # Normalization lambdas
        'l_norm',
        'llow_norm',
        'lupp_norm',
        # Parameter Limits
        'AV_low',
        'AV_upp',
        'YAV_low',
        'YAV_upp',
        'fn_low',
        'fn_upp',
        'v0_low',
        'v0_upp',
        'vd_low',
        'vd_upp',
        # Clipping options & Weight-Control-Filter
        'clip_method_option',
        'sig_clip_threshold',
        'wei_nsig_threshold',
        # Miscellaneous
        'dl_cushion',
        'f_cut',
        'N_int_Gauss',
        'i_verbose',
        'i_verbose_anneal',
        'Is1stLineHeader',
        'i_FastBC03_FLAG',
        'i_FitPowerLaw',
        'alpha_PowerLaw',
        'i_SkipExistingOutFiles',
        # Markov Chains technical parameters
        'N_chains',
        'xinit_max',
        'i_UpdateEps',
        'i_UpdateAlpha',
        'Falpha',
        'i_MoveOneParOnly',
        'i_UpdateAVYAVStepSeparately',
        'i_HelpParWithMaxR',
        'prob_jRmax',
        'i_HelpPopVectorMove2Average',
        'prob_HelpPopVectorMove2Average',
        'i_HelpAVYAVMove2Average',
        'prob_HelpAVYAVMove2Average',
        'NRC_AV_Default',
        # First Fits (FF) technical parameters
        'Temp_ini_FF',
        'Temp_fin_FF',
        'fini_eps_FF',
        'ffin_eps_FF',
        'R_ini_FF',
        'R_fin_FF',
        'IsGRTestHard_FF',
        'N_loops_FF',
        'i_RestartChains_FF',
        'fN_sim_FF',
        'fNmax_steps_FF',
        'eff_IDEAL_FF',
        # GR R-threshold & Method for Burn-In loop
        'R_Burn_in',
        'IsGRTestHard_BurnIn',
        # EX0s technical parameters
        'EXOs_PopVector_option',
        'EXOs_method_option',
        'EXOs_Threshold',
        'Temp_ini_EX0s',
        'Temp_fin_EX0s',
        'fini_eps_EX0s',
        'ffin_eps_EX0s',
        'R_ini_EX0s',
        'R_fin_EX0s',
        'IsGRTestHard_EX0s',
        'N_loops_EX0s',
        'i_RestartChains_EX0s',
        'fN_sim_EX0s',
        'fNmax_steps_EX0s',
        'eff_IDEAL_EX0s',
        'IsScaleNstepsInEX0sFits',
        'IsNoKin4LargeBaseInEX0sFits',
        'frac_NoKin4LargeBaseInEX0sFits',
        'fEX0_MinBaseSize']
    # default values for StarlightConfig instance (StCv04.C99.config)
    meta = _config_StCv04C99_
    # extra comments
    extra = __extra_comments__
    # necessary comments
    config_comments = _config_comments_
    config_comments_insert_index = _config_comments_insert_index_

    def __init__(self, **kwargs):
        """ initialize StarlightConfig instance with StCv04.C99.config """
        self.set_meta(**kwargs)

    def set_meta(self, **kwargs):
        """ set meta data """
        for key in kwargs.keys():
            self.meta[key] = kwargs[key]

    def set_quick(self, template='StCv04.C99.config'):
        if template == 'StCv04.C99.config':
            self.set_meta(**_config_StCv04C99_)
        elif template == 'StCv04.C11.config':
            self.set_meta(**_config_StCv04C11_)
        else:
            raise(ValueError('@Cham: your option should be one of {''StCv04.C99.config'', ''StCv04.C11.config''}'))

    def _meta_to_string(self, val_width=50):
        """ convert meta data to string list """
        fmt_str = "%%%ds" % (-val_width)
        return [(fmt_str + "[%s]\n") % (self.meta[key], key)
                for key in self.meta_order]

    def write(self, filepath, meta_val_width=50):
        f = open(filepath, "w+")
        str_list = self._meta_to_string(meta_val_width)

        # insert necessary comments
        for i in range(len(self.config_comments_insert_index)-1):
            str_list.insert(
                self.config_comments_insert_index[i], self.config_comments[i])
        str_list.append(self.config_comments[-1])

        # for str in str_list:
        #     print(str)

        f.writelines(str_list)
        f.write('\n')
        f.write(self.extra)
        f.close()


def _test_starlight_config():
    sc = StarlightConfig()
    sc.set_quick('StCv04.C99.config')
    sc.meta
    sc.set_quick('StCv04.C11.config')
    sc.meta
    sc.set_meta(Is1stLineHeader=1)
    sc.write('/pool/projects/starlight/STARLIGHTv04/StCv04.C99.config_')


# ################
# Starlight Mask #
# ################


class StarlightMask(object):
    """ StarlightMask class is to represent the Mask file object for STARLIGHT
    """
    pass


def _test_starlight_mask():
    pass


if __name__ =='__main__':
    # _test_starlight_grid()
    # _test_starlight_base()
    # _test_starlight_config()
    _test_starlight_mask()